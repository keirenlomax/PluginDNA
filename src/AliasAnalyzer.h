#pragma once

#include "Analyzer.h"
#include "MeasurementResult.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

class AliasAnalyzer final : public Analyzer
{
public:
    AliasAnalyzer(double sampleRateIn, double durationIn, std::vector<juce::String> paramsIn)
        : sampleRate(sampleRateIn), duration(durationIn), params(std::move(paramsIn))
    {
        result.analyserName = "AliasDNA";
    }

    static const std::vector<double>& testFrequencies()
    {
        static const std::vector<double> frequencies {
            500.0, 1000.0, 2000.0, 5000.0, 8000.0,
            10000.0, 12000.0, 15000.0, 18000.0, 20000.0
        };
        return frequencies;
    }

    void processBlock(const BlockContext& ctx) override
    {
        auto& run = runs[ctx.runId];
        if (!run.ready)
        {
            run.ready = true;
            run.inputGainDb = ctx.inputGainDb;
            run.values = ctx.paramNamedValues;
            run.output.reserve(static_cast<std::size_t>(duration * sampleRate));
        }

        for (int i = 0; i < ctx.numSamples; ++i)
            run.output.push_back(ctx.outL != nullptr ? ctx.outL[i] : 0.0f);
    }

    void finish(const juce::File&) override { build(); }
    const MeasurementDataset& getResult() const noexcept { return result; }
    MeasurementDataset takeResult() noexcept { return std::move(result); }

private:
    struct Run
    {
        bool ready = false;
        float inputGainDb = 0.0f;
        std::map<juce::String, float> values;
        std::vector<float> output;
    };

    static double amplitudeToDb(double value)
    {
        return 20.0 * std::log10(std::max(value, 1.0e-15));
    }

    static double foldToNyquist(double frequency, double sampleRateIn)
    {
        double folded = std::fmod(std::abs(frequency), sampleRateIn);
        if (folded > sampleRateIn * 0.5)
            folded = sampleRateIn - folded;
        return std::abs(folded);
    }

    void build()
    {
        result.clear();
        result.analyserName = "AliasDNA";
        result.columns = {
            "runId", "fundamentalHz", "nyquistHz", "fundamentalMagnitude",
            "legalHarmonicPowerRatio", "aliasPowerRatio", "aliasPowerDb",
            "aliasToHarmonicRatio", "predictedFoldCount", "detectedFoldCount",
            "strongestAliasFrequencyHz", "strongestAliasRatio", "ladderIntegrity",
            "cleanAtMinus90Db", "inputGainDb"
        };
        for (const auto& p : params)
            result.columns.push_back(p.toStdString());

        const auto& frequencies = testFrequencies();
        const int count = static_cast<int>(frequencies.size());
        const int64_t segmentSamples = std::max<int64_t>(1, static_cast<int64_t>(duration * sampleRate / count));
        constexpr int fftOrder = 13;
        constexpr int fftSize = 1 << fftOrder;
        const double binHz = sampleRate / static_cast<double>(fftSize);
        juce::dsp::FFT fft(fftOrder);

        for (const auto& [runId, run] : runs)
        {
            for (int index = 0; index < count; ++index)
            {
                const double fundamental = frequencies[static_cast<std::size_t>(index)];
                if (fundamental >= sampleRate * 0.49)
                    continue;

                const int64_t segmentStart = static_cast<int64_t>(index) * segmentSamples;
                const int64_t segmentEnd = std::min<int64_t>(segmentStart + segmentSamples,
                                                              static_cast<int64_t>(run.output.size()));
                if (segmentEnd - segmentStart < fftSize)
                    continue;

                // Use the final FFT window of each segment, after plugin state has settled.
                const int64_t windowStart = segmentEnd - fftSize;
                std::vector<std::complex<float>> spectrum(static_cast<std::size_t>(fftSize));
                for (int i = 0; i < fftSize; ++i)
                {
                    const double window = 0.5 * (1.0 - std::cos(2.0 * juce::MathConstants<double>::pi
                                                               * i / static_cast<double>(fftSize - 1)));
                    spectrum[static_cast<std::size_t>(i)] = {
                        static_cast<float>(run.output[static_cast<std::size_t>(windowStart + i)] * window), 0.0f
                    };
                }
                fft.perform(spectrum.data(), spectrum.data(), false);

                auto bandMagnitude = [&](double frequency)
                {
                    const int centre = static_cast<int>(std::round(frequency / binHz));
                    double power = 0.0;
                    for (int offset = -2; offset <= 2; ++offset)
                    {
                        const int bin = centre + offset;
                        if (bin <= 0 || bin >= fftSize / 2)
                            continue;
                        const double magnitude = std::abs(spectrum[static_cast<std::size_t>(bin)]);
                        power += magnitude * magnitude;
                    }
                    return std::sqrt(power);
                };

                const double fundamentalMagnitude = bandMagnitude(fundamental);
                if (fundamentalMagnitude <= 1.0e-12)
                    continue;

                std::vector<double> legalFrequencies;
                legalFrequencies.push_back(fundamental);
                double legalPower = 0.0;
                for (int harmonic = 2; harmonic <= 10; ++harmonic)
                {
                    const double frequency = fundamental * harmonic;
                    if (frequency < sampleRate * 0.5)
                    {
                        legalFrequencies.push_back(frequency);
                        const double ratio = bandMagnitude(frequency) / fundamentalMagnitude;
                        legalPower += ratio * ratio;
                    }
                }

                double aliasPower = 0.0;
                double strongestAlias = 0.0;
                double strongestAliasFrequency = 0.0;
                int predictedFoldCount = 0;
                int detectedFoldCount = 0;

                for (int harmonic = 2; harmonic <= 10; ++harmonic)
                {
                    const double rawFrequency = fundamental * harmonic;
                    if (rawFrequency < sampleRate * 0.5)
                        continue;

                    ++predictedFoldCount;
                    const double folded = foldToNyquist(rawFrequency, sampleRate);
                    if (folded <= binHz * 3.0 || folded >= sampleRate * 0.5 - binHz * 3.0)
                        continue;

                    const bool overlapsLegal = std::any_of(legalFrequencies.begin(), legalFrequencies.end(),
                        [&](double legal) { return std::abs(legal - folded) <= binHz * 4.0; });
                    if (overlapsLegal)
                        continue; // Cannot separate alias energy from a legitimate harmonic at this bin.

                    const double ratio = bandMagnitude(folded) / fundamentalMagnitude;
                    aliasPower += ratio * ratio;
                    if (ratio >= std::pow(10.0, -90.0 / 20.0))
                        ++detectedFoldCount;
                    if (ratio > strongestAlias)
                    {
                        strongestAlias = ratio;
                        strongestAliasFrequency = folded;
                    }
                }

                const double aliasRatio = std::sqrt(aliasPower);
                const double legalRatio = std::sqrt(legalPower);
                const double aliasToHarmonic = legalRatio > 1.0e-15 ? aliasRatio / legalRatio : 0.0;
                const double integrity = 1.0 / (1.0 + aliasToHarmonic);

                std::vector<double> row {
                    static_cast<double>(runId), fundamental, sampleRate * 0.5, fundamentalMagnitude,
                    legalRatio, aliasRatio, amplitudeToDb(aliasRatio), aliasToHarmonic,
                    static_cast<double>(predictedFoldCount), static_cast<double>(detectedFoldCount),
                    strongestAliasFrequency, strongestAlias, integrity,
                    aliasRatio < std::pow(10.0, -90.0 / 20.0) ? 1.0 : 0.0,
                    run.inputGainDb
                };
                for (const auto& p : params)
                {
                    const auto it = run.values.find(p);
                    row.push_back(it != run.values.end() ? it->second : 0.0);
                }
                result.rows.push_back(std::move(row));
            }
        }
    }

    double sampleRate = 96000.0;
    double duration = 5.0;
    std::vector<juce::String> params;
    std::map<int, Run> runs;
    MeasurementDataset result;
};

inline std::unique_ptr<Analyzer> createAliasAnalyzer(double sampleRate, double duration,
                                                      const std::vector<juce::String>& params)
{
    return std::make_unique<AliasAnalyzer>(sampleRate, duration, params);
}
