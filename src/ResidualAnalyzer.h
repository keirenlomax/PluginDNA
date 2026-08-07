#pragma once

#include "Analyzer.h"
#include "MeasurementResult.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>

class ResidualAnalyzer final : public Analyzer
{
public:
    ResidualAnalyzer(double sampleRateIn, std::vector<juce::String> parameterNamesIn,
                     int captureSamplesIn = 65536)
        : sampleRate(sampleRateIn), parameterNames(std::move(parameterNamesIn)),
          captureSamples(captureSamplesIn)
    {
        result.analyserName = "ResidualDNA";
    }

    void processBlock(const BlockContext& ctx) override
    {
        auto& run = runs[ctx.runId];
        if (!run.initialised)
        {
            run.initialised = true;
            run.inputGainDb = ctx.inputGainDb;
            run.paramValues = ctx.paramNamedValues;
            run.input.reserve((std::size_t)captureSamples);
            run.output.reserve((std::size_t)captureSamples);
        }
        if ((int)run.input.size() >= captureSamples) return;
        const int remaining = captureSamples - (int)run.input.size();
        const int count = std::min(remaining, ctx.numSamples);
        for (int i = 0; i < count; ++i)
        {
            run.input.push_back(ctx.inL ? ctx.inL[i] : 0.0f);
            run.output.push_back(ctx.outL ? ctx.outL[i] : 0.0f);
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
        std::vector<float> input, output;
    };

    static double toDb(double x) { return 20.0 * std::log10(std::max(x, 1.0e-12)); }

    int bestLag(const std::vector<float>& x, const std::vector<float>& y) const
    {
        int winner = 0; double best = -1.0;
        const int n = (int)std::min(x.size(), y.size());
        const int start = std::min(4096, n / 8);
        const int end = std::min(n, start + 32768);
        for (int lag = -128; lag <= 128; ++lag)
        {
            double xy = 0.0, xx = 0.0, yy = 0.0;
            for (int i = start; i < end; ++i)
            {
                const int j = i + lag;
                if (j < 0 || j >= n) continue;
                const double a = x[(std::size_t)i], b = y[(std::size_t)j];
                xy += a * b; xx += a * a; yy += b * b;
            }
            const double corr = std::abs(xy) / std::sqrt(std::max(xx * yy, 1.0e-30));
            if (corr > best) { best = corr; winner = lag; }
        }
        return winner;
    }

    struct ResidualSpectrum
    {
        std::array<double, 4> bandDb {-240.0, -240.0, -240.0, -240.0};
        double centroidHz = 0.0;
        double spreadHz = 0.0;
        double flatness = 0.0;
        double tonality = 0.0;
        double entropyBits = 0.0;
        double entropyNormalized = 0.0;
        double effectiveBandwidthHz = 0.0;
        double dominantFrequencyHz = 0.0;
        double dominantFrequencyDb = -240.0;
    };

    ResidualSpectrum analyseResidualSpectrum(const std::vector<double>& residual) const
    {
        const int fftSize = 65536;
        std::vector<std::complex<float>> time((std::size_t)fftSize), freq((std::size_t)fftSize);
        const int n = std::min((int)residual.size(), fftSize);
        for (int i = 0; i < n; ++i)
        {
            const double w = 0.5 - 0.5 * std::cos(2.0 * juce::MathConstants<double>::pi * i / (fftSize - 1));
            time[(std::size_t)i] = {(float)(residual[(std::size_t)i] * w), 0.0f};
        }
        juce::dsp::FFT fft(16);
        fft.perform(time.data(), freq.data(), false);
        ResidualSpectrum spectrum;
        std::array<double, 4> energy {0,0,0,0};
        std::array<int, 4> count {0,0,0,0};
        double totalPower = 0.0, weightedHz = 0.0, weightedHz2 = 0.0;
        double logMagnitudeSum = 0.0, magnitudeSum = 0.0;
        std::vector<double> binPowers;
        binPowers.reserve(static_cast<std::size_t>(fftSize / 2));
        int spectralBins = 0;
        double strongestPower = 0.0;
        for (int k = 1; k < fftSize / 2; ++k)
        {
            const double hz = k * sampleRate / fftSize;
            int band = hz < 200.0 ? 0 : hz < 4000.0 ? 1 : hz < 20000.0 ? 2 : 3;
            const double power = std::norm(freq[(std::size_t)k]);
            const double magnitude = std::sqrt(power);
            energy[(std::size_t)band] += power;
            count[(std::size_t)band]++;
            totalPower += power;
            binPowers.push_back(power);
            weightedHz += hz * power;
            weightedHz2 += hz * hz * power;
            logMagnitudeSum += std::log(std::max(magnitude, 1.0e-30));
            magnitudeSum += magnitude;
            ++spectralBins;
            if (power > strongestPower)
            {
                strongestPower = power;
                spectrum.dominantFrequencyHz = hz;
            }
        }
        for (int b = 0; b < 4; ++b)
            spectrum.bandDb[(std::size_t)b] = count[(std::size_t)b] > 0
                ? toDb(std::sqrt(energy[(std::size_t)b] / count[(std::size_t)b]) / fftSize)
                : -240.0;
        if (totalPower > 0.0)
        {
            spectrum.centroidHz = weightedHz / totalPower;
            spectrum.spreadHz = std::sqrt(std::max(0.0, weightedHz2 / totalPower
                                                       - spectrum.centroidHz * spectrum.centroidHz));
            spectrum.effectiveBandwidthHz = 2.0 * spectrum.spreadHz;
            for (const double power : binPowers)
            {
                const double probability = power / totalPower;
                if (probability > 0.0)
                    spectrum.entropyBits -= probability * std::log2(probability);
            }
            if (binPowers.size() > 1)
                spectrum.entropyNormalized = juce::jlimit(0.0, 1.0,
                    spectrum.entropyBits / std::log2(static_cast<double>(binPowers.size())));
        }
        if (spectralBins > 0 && magnitudeSum > 0.0)
        {
            const double geometric = std::exp(logMagnitudeSum / spectralBins);
            const double arithmetic = magnitudeSum / spectralBins;
            spectrum.flatness = juce::jlimit(0.0, 1.0, geometric / std::max(arithmetic, 1.0e-30));
            spectrum.tonality = 1.0 - spectrum.flatness;
        }
        spectrum.dominantFrequencyDb = strongestPower > 0.0
            ? toDb(std::sqrt(strongestPower) / fftSize) : -240.0;
        return spectrum;
    }

    void buildResult()
    {
        result.columns = {"runId","lagSamples","gainDb","residualDbRelativeToOutput",
                          "lowResidualDb","midResidualDb","highResidualDb","ultrasonicResidualDb",
                          "residualPeakDbRelativeToOutput","residualCrestDb","residualMeanAbs",
                          "residualZeroCrossingRate","residualKurtosis","residualSparsity",
                          "residualSpectralCentroidHz","residualSpectralSpreadHz","residualEffectiveBandwidthHz",
                          "residualSpectralEntropyBits","residualSpectralEntropyNormalized","residualSpectralFlatness",
                          "residualTonality","residualDominantFrequencyHz","residualDominantFrequencyDb",
                          "inputGainDb"};
        for (const auto& p : parameterNames) result.columns.push_back(p.toStdString());
        for (const auto& [runId, run] : runs)
        {
            const int n = (int)std::min(run.input.size(), run.output.size());
            if (n < 4096) continue;
            const int lag = bestLag(run.input, run.output);
            double xy = 0.0, xx = 0.0, yy = 0.0;
            for (int i = 0; i < n; ++i)
            {
                const int j = i + lag; if (j < 0 || j >= n) continue;
                const double x = run.input[(std::size_t)i], y = run.output[(std::size_t)j];
                xy += x*y; xx += x*x; yy += y*y;
            }
            const double gain = xx > 0.0 ? xy / xx : 1.0;
            std::vector<double> residual; residual.reserve((std::size_t)n);
            double rr = 0.0, absoluteSum = 0.0, fourthMoment = 0.0, residualPeak = 0.0;
            int zeroCrossings = 0; int used = 0; double previousResidual = 0.0;
            for (int i = 0; i < n; ++i)
            {
                const int j = i + lag; if (j < 0 || j >= n) continue;
                const double r = run.output[(std::size_t)j] - gain * run.input[(std::size_t)i];
                residual.push_back(r); rr += r*r; absoluteSum += std::abs(r);
                fourthMoment += r*r*r*r; residualPeak = std::max(residualPeak, std::abs(r));
                if (used > 0 && ((r >= 0.0) != (previousResidual >= 0.0))) ++zeroCrossings;
                previousResidual = r; used++;
            }
            const double residualRms = used ? std::sqrt(rr / used) : 0.0;
            const double outputRms = used ? std::sqrt(yy / used) : 0.0;
            const auto spectrum = analyseResidualSpectrum(residual);
            const double meanAbs = used ? absoluteSum / used : 0.0;
            const double crestDb = residualRms > 0.0 ? toDb(residualPeak / residualRms) : 0.0;
            const double zeroCrossingRate = used > 1 ? static_cast<double>(zeroCrossings) / (used - 1) : 0.0;
            const double kurtosis = residualRms > 0.0 && used > 0
                ? (fourthMoment / used) / std::pow(residualRms, 4.0) : 0.0;
            const double sparsity = residualPeak > 0.0 ? 1.0 - meanAbs / residualPeak : 0.0;
            std::vector<double> row { (double)runId, (double)lag, toDb(std::abs(gain)),
                toDb(residualRms / std::max(outputRms, 1.0e-12)),
                spectrum.bandDb[0], spectrum.bandDb[1], spectrum.bandDb[2], spectrum.bandDb[3],
                toDb(residualPeak / std::max(outputRms, 1.0e-12)), crestDb, meanAbs,
                zeroCrossingRate, kurtosis, sparsity, spectrum.centroidHz, spectrum.spreadHz,
                spectrum.effectiveBandwidthHz, spectrum.entropyBits, spectrum.entropyNormalized,
                spectrum.flatness, spectrum.tonality, spectrum.dominantFrequencyHz,
                spectrum.dominantFrequencyDb, run.inputGainDb };
            for (const auto& p : parameterNames)
            {
                auto it = run.paramValues.find(p); row.push_back(it != run.paramValues.end() ? it->second : 0.0);
            }
            result.rows.push_back(std::move(row));
        }
    }

    double sampleRate;
    std::vector<juce::String> parameterNames;
    int captureSamples;
    std::map<int, RunState> runs;
    MeasurementDataset result;
};

inline std::unique_ptr<Analyzer> createResidualAnalyzer(double sampleRate,
    const std::vector<juce::String>& parameterNames)
{
    return std::make_unique<ResidualAnalyzer>(sampleRate, parameterNames);
}
