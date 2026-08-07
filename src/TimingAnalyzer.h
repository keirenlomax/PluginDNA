#pragma once

#include "Analyzer.h"
#include "MeasurementResult.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>

class TimingAnalyzer final : public Analyzer
{
public:
    TimingAnalyzer(double sampleRateIn, int impulseSampleIn,
                   std::vector<juce::String> parameterNamesIn,
                   int captureSamplesIn = 32768)
        : sampleRate(sampleRateIn), impulseSample(impulseSampleIn),
          parameterNames(std::move(parameterNamesIn)), captureSamples(captureSamplesIn)
    {
        result.analyserName = "TimingDNA";
    }

    void processBlock(const BlockContext& ctx) override
    {
        auto& run = runs[ctx.runId];
        if (!run.initialised)
        {
            run.initialised = true;
            run.inputGainDb = ctx.inputGainDb;
            run.paramValues = ctx.paramNamedValues;
            run.input.assign(static_cast<std::size_t>(captureSamples), 0.0f);
            run.output.assign(static_cast<std::size_t>(captureSamples), 0.0f);
        }

        for (int i = 0; i < ctx.numSamples; ++i)
        {
            const auto absolute = ctx.firstSample + i;
            if (absolute < 0 || absolute >= captureSamples)
                continue;
            const auto index = static_cast<std::size_t>(absolute);
            run.input[index] = ctx.inL != nullptr ? ctx.inL[i] : 0.0f;
            run.output[index] = ctx.outL != nullptr ? ctx.outL[i] : 0.0f;
        }
    }

    void finish(const juce::File&) override { buildResult(); }

    const MeasurementDataset& getResult() const noexcept { return result; }
    MeasurementDataset takeResult() noexcept { return std::move(result); }

private:
    struct RunState
    {
        bool initialised = false;
        float inputGainDb = 0.0f;
        std::map<juce::String, float> paramValues;
        std::vector<float> input;
        std::vector<float> output;
    };

    static double amplitudeToDb(double amplitude)
    {
        return 20.0 * std::log10(std::max(amplitude, 1.0e-12));
    }

    static double wrapDegrees(double degrees)
    {
        while (degrees > 180.0) degrees -= 360.0;
        while (degrees <= -180.0) degrees += 360.0;
        return degrees;
    }

    void appendParameters(std::vector<double>& row, const RunState& run) const
    {
        for (const auto& name : parameterNames)
        {
            const auto it = run.paramValues.find(name);
            row.push_back(it != run.paramValues.end() ? static_cast<double>(it->second) : 0.0);
        }
        row.push_back(static_cast<double>(run.inputGainDb));
    }

    void buildResult()
    {
        result.clear();
        result.analyserName = "TimingDNA";
        result.columns = {
            "runId", "recordType", "frequencyHz", "magnitudeDb", "phaseDegrees",
            "latencySamples", "latencyMs", "preRingingDb", "postRingingDb",
            "settlingMs", "overshootDb", "undershootDb", "riseTimeSamples",
            "riseTimeMs", "fallTimeSamples", "fallTimeMs", "impulseWidthSamples",
            "impulseWidthMs", "energyCentroidMs", "preRingingEnergyDb",
            "postRingingEnergyDb", "earlyEnergyFraction", "lateEnergyFraction"
        };
        for (const auto& name : parameterNames)
            result.columns.push_back(name.toStdString());
        result.columns.push_back("inputGainDb");

        constexpr int fftOrder = 15;
        constexpr int fftSize = 1 << fftOrder;
        juce::dsp::FFT fft(fftOrder);

        for (const auto& [runId, run] : runs)
        {
            if (run.input.size() < static_cast<std::size_t>(fftSize)
                || run.output.size() < static_cast<std::size_t>(fftSize))
                continue;

            int peakIndex = impulseSample;
            double peak = 0.0;
            for (int i = std::max(0, impulseSample); i < fftSize; ++i)
            {
                const double value = std::abs(static_cast<double>(run.output[static_cast<std::size_t>(i)]));
                if (value > peak) { peak = value; peakIndex = i; }
            }

            const double threshold = peak * 0.001; // -60 dB from main peak.
            double pre = 0.0;
            for (int i = impulseSample; i < peakIndex; ++i)
                pre = std::max(pre, std::abs(static_cast<double>(run.output[static_cast<std::size_t>(i)])));

            double post = 0.0;
            int lastAbove = peakIndex;
            for (int i = peakIndex + 1; i < fftSize; ++i)
            {
                const double value = std::abs(static_cast<double>(run.output[static_cast<std::size_t>(i)]));
                post = std::max(post, value);
                if (value > threshold) lastAbove = i;
            }

            const double inputPeak = std::abs(static_cast<double>(run.input[static_cast<std::size_t>(impulseSample)]));
            const double latencySamples = static_cast<double>(peakIndex - impulseSample);
            const double latencyMs = 1000.0 * latencySamples / sampleRate;
            const double preDb = peak > 0.0 ? amplitudeToDb(pre / peak) : -240.0;
            const double postDb = peak > 0.0 ? amplitudeToDb(post / peak) : -240.0;
            const double settlingMs = 1000.0 * static_cast<double>(std::max(0, lastAbove - peakIndex)) / sampleRate;
            const double overshootDb = inputPeak > 0.0 ? amplitudeToDb(peak / inputPeak) : -240.0;

            int rise10 = peakIndex, rise90 = peakIndex, fall90 = peakIndex, fall10 = peakIndex;
            const double ten = peak * 0.10, ninety = peak * 0.90, half = peak * 0.50;
            for (int i = std::max(0, impulseSample); i <= peakIndex; ++i)
            {
                const double a = std::abs(static_cast<double>(run.output[static_cast<std::size_t>(i)]));
                if (a >= ten && rise10 == peakIndex) rise10 = i;
                if (a >= ninety) { rise90 = i; break; }
            }
            bool found90 = false;
            for (int i = peakIndex; i < fftSize; ++i)
            {
                const double a = std::abs(static_cast<double>(run.output[static_cast<std::size_t>(i)]));
                if (!found90 && a <= ninety) { fall90 = i; found90 = true; }
                if (found90 && a <= ten) { fall10 = i; break; }
            }
            int widthStart = peakIndex, widthEnd = peakIndex;
            for (int i = peakIndex; i >= 0; --i)
                if (std::abs(static_cast<double>(run.output[static_cast<std::size_t>(i)])) >= half) widthStart = i; else break;
            for (int i = peakIndex; i < fftSize; ++i)
                if (std::abs(static_cast<double>(run.output[static_cast<std::size_t>(i)])) >= half) widthEnd = i; else break;

            double negativeTrough = 0.0, totalEnergy = 0.0, weightedEnergy = 0.0;
            double preEnergy = 0.0, postEnergy = 0.0, earlyEnergy = 0.0, lateEnergy = 0.0;
            const int earlyEnd = std::min(fftSize, peakIndex + static_cast<int>(0.001 * sampleRate));
            for (int i = 0; i < fftSize; ++i)
            {
                const double value = static_cast<double>(run.output[static_cast<std::size_t>(i)]);
                negativeTrough = std::min(negativeTrough, value);
                const double e = value * value;
                totalEnergy += e;
                weightedEnergy += (i - impulseSample) * e;
                if (i < peakIndex) preEnergy += e; else if (i > peakIndex) postEnergy += e;
                if (i >= peakIndex && i < earlyEnd) earlyEnergy += e;
                if (i >= earlyEnd) lateEnergy += e;
            }
            const double riseSamples = std::max(0, rise90 - rise10);
            const double fallSamples = std::max(0, fall10 - fall90);
            const double widthSamples = std::max(0, widthEnd - widthStart + 1);
            const double energyCentroidMs = totalEnergy > 0.0
                ? 1000.0 * (weightedEnergy / totalEnergy) / sampleRate : 0.0;
            const double undershootDb = peak > 0.0 ? amplitudeToDb(std::abs(negativeTrough) / peak) : -240.0;
            const double preEnergyDb = totalEnergy > 0.0 ? 10.0 * std::log10(std::max(preEnergy / totalEnergy, 1.0e-24)) : -240.0;
            const double postEnergyDb = totalEnergy > 0.0 ? 10.0 * std::log10(std::max(postEnergy / totalEnergy, 1.0e-24)) : -240.0;
            const double earlyFraction = totalEnergy > 0.0 ? earlyEnergy / totalEnergy : 0.0;
            const double lateFraction = totalEnergy > 0.0 ? lateEnergy / totalEnergy : 0.0;

            std::vector<double> summary {
                static_cast<double>(runId), 2.0, 0.0, 0.0, 0.0,
                latencySamples, latencyMs, preDb, postDb, settlingMs, overshootDb,
                undershootDb, riseSamples, 1000.0 * riseSamples / sampleRate,
                fallSamples, 1000.0 * fallSamples / sampleRate,
                widthSamples, 1000.0 * widthSamples / sampleRate, energyCentroidMs,
                preEnergyDb, postEnergyDb, earlyFraction, lateFraction
            };
            appendParameters(summary, run);
            result.rows.push_back(std::move(summary));

            std::vector<std::complex<float>> inputTime(fftSize), outputTime(fftSize);
            std::vector<std::complex<float>> inputFreq(fftSize), outputFreq(fftSize);
            for (int i = 0; i < fftSize; ++i)
            {
                inputTime[static_cast<std::size_t>(i)] = { run.input[static_cast<std::size_t>(i)], 0.0f };
                outputTime[static_cast<std::size_t>(i)] = { run.output[static_cast<std::size_t>(i)], 0.0f };
            }
            fft.perform(inputTime.data(), inputFreq.data(), false);
            fft.perform(outputTime.data(), outputFreq.data(), false);

            const double binHz = sampleRate / static_cast<double>(fftSize);
            for (int bin = 1; bin < fftSize / 2; ++bin)
            {
                const double frequency = bin * binHz;
                if (frequency < 20.0 || frequency > std::min(20000.0, sampleRate * 0.48))
                    continue;
                const auto in = inputFreq[static_cast<std::size_t>(bin)];
                const auto out = outputFreq[static_cast<std::size_t>(bin)];
                if (std::abs(in) < 1.0e-12f)
                    continue;
                const std::complex<double> transfer = std::complex<double>(out.real(), out.imag())
                                                    / std::complex<double>(in.real(), in.imag());
                const double magDb = amplitudeToDb(std::abs(transfer));
                const double phaseDeg = wrapDegrees(std::arg(transfer) * 180.0
                                                    / juce::MathConstants<double>::pi);
                std::vector<double> row {
                    static_cast<double>(runId), 1.0, frequency, magDb, phaseDeg,
                    latencySamples, latencyMs, preDb, postDb, settlingMs, overshootDb,
                    undershootDb, riseSamples, 1000.0 * riseSamples / sampleRate,
                    fallSamples, 1000.0 * fallSamples / sampleRate,
                    widthSamples, 1000.0 * widthSamples / sampleRate, energyCentroidMs,
                    preEnergyDb, postEnergyDb, earlyFraction, lateFraction
                };
                appendParameters(row, run);
                result.rows.push_back(std::move(row));
            }
        }
    }

    double sampleRate = 96000.0;
    int impulseSample = 8192;
    std::vector<juce::String> parameterNames;
    int captureSamples = 32768;
    std::map<int, RunState> runs;
    MeasurementDataset result;
};

inline std::unique_ptr<Analyzer> createTimingAnalyzer(double sampleRate, int impulseSample,
                                                       const std::vector<juce::String>& parameterNames)
{
    return std::make_unique<TimingAnalyzer>(sampleRate, impulseSample, parameterNames);
}
