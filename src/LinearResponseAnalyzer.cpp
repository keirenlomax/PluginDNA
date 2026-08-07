#include "LinearResponseAnalyzer.h"
#include "JuceHeader.h"
#include <cmath>
#include <fstream>
#include <iostream>

LinearResponseAnalyzer::LinearResponseAnalyzer(const juce::File& outDir, int fftSize,
                                               const std::vector<juce::String>& paramNames,
                                               const juce::String& signalType,
                                               double sweepStartHzIn,
                                               double sweepEndHzIn,
                                               double sweepDurationSecondsIn)
    : fftSize(fftSize), paramNames(paramNames), outputDir(outDir), signalType(signalType),
      sweepStartHz(sweepStartHzIn), sweepEndHz(sweepEndHzIn),
      sweepDurationSeconds(sweepDurationSecondsIn) {}

LinearResponseAnalyzer::~LinearResponseAnalyzer() {}

void LinearResponseAnalyzer::applyHannWindow(std::vector<float>& buffer) {
    const int N = (int)buffer.size();
    for (int i = 0; i < N; ++i) {
        float window = 0.5f * (1.0f - std::cos(2.0f * juce::MathConstants<float>::pi * (float)i / (float)(N - 1)));
        buffer[i] *= window;
    }
}

void LinearResponseAnalyzer::processFFTWindow(RunSpectrum& spectrum) {
    if ((int)spectrum.inBuffer.size() < fftSize || (int)spectrum.outBuffer.size() < fftSize)
        return;

    // Apply window
    applyHannWindow(spectrum.inBuffer);
    applyHannWindow(spectrum.outBuffer);

    // Transform input and plugin output independently.  The previous code used
    // outFFT as the destination for the input transform, so the plugin output
    // was never transformed and the resulting transfer magnitudes were invalid.
    juce::dsp::FFT fft((int)std::log2(fftSize));
    std::vector<std::complex<float>> inTime(fftSize);
    std::vector<std::complex<float>> outTime(fftSize);
    std::vector<std::complex<float>> inFFT(fftSize);
    std::vector<std::complex<float>> outFFT(fftSize);

    for (int i = 0; i < fftSize; ++i) {
        inTime[i] = std::complex<float>(spectrum.inBuffer[i], 0.0f);
        outTime[i] = std::complex<float>(spectrum.outBuffer[i], 0.0f);
    }

    fft.perform(inTime.data(), inFFT.data(), false);
    fft.perform(outTime.data(), outFFT.data(), false);

    // Accumulate magnitude squared
    const int numBins = fftSize / 2;
    if ((int)spectrum.sumInMagSq.size() < numBins) {
        spectrum.sumInMagSq.resize(numBins, 0.0);
        spectrum.sumOutMagSq.resize(numBins, 0.0);
    }

    for (int k = 0; k < numBins; ++k) {
        float inMag = std::abs(inFFT[k]);
        float outMag = std::abs(outFFT[k]);
        spectrum.sumInMagSq[k] += (double)(inMag * inMag);
        spectrum.sumOutMagSq[k] += (double)(outMag * outMag);
    }

    spectrum.numAverages++;

    // Clear buffers for next window
    spectrum.inBuffer.clear();
    spectrum.outBuffer.clear();
}


void LinearResponseAnalyzer::processHarmonicFingerprintSamples(RunSpectrum& spectrum,
                                                                const BlockContext& ctx)
{
    if (!signalType.equalsIgnoreCase("sweep") || sweepStartHz <= 0.0
        || sweepEndHz <= sweepStartHz || sweepDurationSeconds <= 0.0)
        return;

    const int bandCount = RunSpectrum::harmonicBands;
    const int maxOrder = RunSpectrum::maxHarmonicOrder;

    if (spectrum.harmonicCounts.empty())
    {
        spectrum.inputFundSin.assign((std::size_t)bandCount, 0.0);
        spectrum.inputFundCos.assign((std::size_t)bandCount, 0.0);
        spectrum.outputHarmSin.assign((std::size_t)maxOrder,
                                      std::vector<double>((std::size_t)bandCount, 0.0));
        spectrum.outputHarmCos.assign((std::size_t)maxOrder,
                                      std::vector<double>((std::size_t)bandCount, 0.0));
        spectrum.harmonicCounts.assign((std::size_t)bandCount, 0);
        spectrum.sweepPhase = 0.0;
    }

    const double totalSamples = sweepDurationSeconds * ctx.sampleRate;
    const double logStart = std::log(sweepStartHz);
    const double logRange = std::log(sweepEndHz / sweepStartHz);
    const double twoPi = 2.0 * juce::MathConstants<double>::pi;

    for (int i = 0; i < ctx.numSamples; ++i)
    {
        const double absoluteSample = static_cast<double>(ctx.firstSample + i);
        if (absoluteSample >= totalSamples)
            break;

        const double position = juce::jlimit(0.0, 1.0, absoluteSample / totalSamples);
        const double frequency = std::exp(logStart + position * logRange);
        const int band = juce::jlimit(0, bandCount - 1,
            static_cast<int>(std::floor(position * static_cast<double>(bandCount))));

        const double input = static_cast<double>(ctx.inL[i]);
        const double output = static_cast<double>(ctx.outL[i]);

        spectrum.inputFundSin[(std::size_t)band] += input * std::sin(spectrum.sweepPhase);
        spectrum.inputFundCos[(std::size_t)band] += input * std::cos(spectrum.sweepPhase);

        for (int order = 1; order <= maxOrder; ++order)
        {
            if (frequency * static_cast<double>(order) >= 0.48 * ctx.sampleRate)
                continue;

            const double harmonicPhase = static_cast<double>(order) * spectrum.sweepPhase;
            spectrum.outputHarmSin[(std::size_t)(order - 1)][(std::size_t)band]
                += output * std::sin(harmonicPhase);
            spectrum.outputHarmCos[(std::size_t)(order - 1)][(std::size_t)band]
                += output * std::cos(harmonicPhase);
        }

        spectrum.harmonicCounts[(std::size_t)band]++;
        spectrum.sweepPhase += twoPi * frequency / ctx.sampleRate;
        if (spectrum.sweepPhase >= twoPi)
            spectrum.sweepPhase = std::fmod(spectrum.sweepPhase, twoPi);
    }
}

void LinearResponseAnalyzer::processBlock(const BlockContext& ctx) {
    auto& spectrum = perRunSpectra[ctx.runId];

    // Initialize on first block
    if (spectrum.sumInMagSq.empty()) {
        spectrum.paramValues = ctx.paramNamedValues;
        spectrum.inputGainDb = ctx.inputGainDb;
        spectrum.sampleRate = ctx.sampleRate;
    }

    processHarmonicFingerprintSamples(spectrum, ctx);

    // Accumulate samples
    for (int i = 0; i < ctx.numSamples; ++i) {
        spectrum.inBuffer.push_back(ctx.inL[i]);
        spectrum.outBuffer.push_back(ctx.outL[i]);

        // Process FFT window when we have enough samples
        if ((int)spectrum.inBuffer.size() >= fftSize) {
            processFFTWindow(spectrum);
        }
    }
}

void LinearResponseAnalyzer::buildResult()
{
    result.clear();
    result.analyserName = "LinearResponse";

    result.columns.push_back("runId");
    result.columns.push_back("freqHz");
    result.columns.push_back("magDb");

    for (const auto& paramName : paramNames)
        result.columns.push_back(paramName.toStdString());

    result.columns.push_back("inputGainDb");

    for (const auto& [runId, spectrum] : perRunSpectra)
    {
        if (spectrum.numAverages == 0)
            continue;

        const int numBins = fftSize / 2;
        const double binHz = spectrum.sampleRate / static_cast<double>(fftSize);

        for (int bin = 0; bin < numBins; ++bin)
        {
            const double magIn = std::sqrt(
                spectrum.sumInMagSq[static_cast<std::size_t>(bin)]
                / static_cast<double>(spectrum.numAverages));

            const double magOut = std::sqrt(
                spectrum.sumOutMagSq[static_cast<std::size_t>(bin)]
                / static_cast<double>(spectrum.numAverages));

            if (magIn <= 0.0)
                continue;

            const double transferMagnitude = magOut / magIn;
            const double magDb = 20.0 * std::log10(
                std::max(transferMagnitude, 1.0e-10));

            const double freqHz = static_cast<double>(bin) * binHz;

            std::vector<double> row;
            row.reserve(result.columns.size());

            row.push_back(static_cast<double>(runId));
            row.push_back(freqHz);
            row.push_back(magDb);

            for (const auto& paramName : paramNames)
            {
                double value = 0.0;
                const auto valueIt = spectrum.paramValues.find(paramName);

                if (valueIt != spectrum.paramValues.end())
                    value = static_cast<double>(valueIt->second);

                row.push_back(value);
            }

            row.push_back(static_cast<double>(spectrum.inputGainDb));
            result.rows.push_back(std::move(row));
        }
    }
}


void LinearResponseAnalyzer::buildHarmonicFingerprintResult()
{
    harmonicFingerprintResult.clear();
    harmonicFingerprintResult.analyserName = "HarmonicFingerprint";
    harmonicFingerprintResult.columns = { "runId", "fundamentalHz", "harmonicOrder",
                                          "magDbRelativeToInput" };
    for (const auto& paramName : paramNames)
        harmonicFingerprintResult.columns.push_back(paramName.toStdString());
    harmonicFingerprintResult.columns.push_back("inputGainDb");

    if (!signalType.equalsIgnoreCase("sweep"))
        return;

    const int bandCount = RunSpectrum::harmonicBands;
    const int maxOrder = RunSpectrum::maxHarmonicOrder;
    const double logRange = std::log(sweepEndHz / sweepStartHz);

    for (const auto& [runId, spectrum] : perRunSpectra)
    {
        if ((int)spectrum.harmonicCounts.size() != bandCount)
            continue;

        for (int band = 0; band < bandCount; ++band)
        {
            const auto count = spectrum.harmonicCounts[(std::size_t)band];
            if (count < 32)
                continue;

            const double centrePosition = (static_cast<double>(band) + 0.5)
                                        / static_cast<double>(bandCount);
            const double fundamentalHz = sweepStartHz * std::exp(centrePosition * logRange);
            const double inputAmplitude = (2.0 / static_cast<double>(count))
                * std::hypot(spectrum.inputFundSin[(std::size_t)band],
                             spectrum.inputFundCos[(std::size_t)band]);
            if (inputAmplitude <= 1.0e-12)
                continue;

            for (int order = 1; order <= maxOrder; ++order)
            {
                if (fundamentalHz * static_cast<double>(order) >= 0.48 * spectrum.sampleRate)
                    continue;

                const double outputAmplitude = (2.0 / static_cast<double>(count))
                    * std::hypot(spectrum.outputHarmSin[(std::size_t)(order - 1)][(std::size_t)band],
                                 spectrum.outputHarmCos[(std::size_t)(order - 1)][(std::size_t)band]);
                const double relative = outputAmplitude / inputAmplitude;
                const double magDb = 20.0 * std::log10(std::max(relative, 1.0e-12));

                std::vector<double> row;
                row.reserve(harmonicFingerprintResult.columns.size());
                row.push_back(static_cast<double>(runId));
                row.push_back(fundamentalHz);
                row.push_back(static_cast<double>(order));
                row.push_back(magDb);
                for (const auto& paramName : paramNames)
                {
                    const auto it = spectrum.paramValues.find(paramName);
                    row.push_back(it != spectrum.paramValues.end()
                        ? static_cast<double>(it->second) : 0.0);
                }
                row.push_back(static_cast<double>(spectrum.inputGainDb));
                harmonicFingerprintResult.rows.push_back(std::move(row));
            }
        }
    }
}

void LinearResponseAnalyzer::finish(const juce::File& outDir)
{
    juce::ignoreUnused(outDir);
    buildResult();
    buildHarmonicFingerprintResult();
}

std::unique_ptr<Analyzer> createLinearResponseAnalyzer(const juce::File& outDir, int fftSize,
                                                       const std::vector<juce::String>& paramNames,
                                                       const juce::String& signalType,
                                                       double sweepStartHz,
                                                       double sweepEndHz,
                                                       double sweepDurationSeconds) {
    return std::make_unique<LinearResponseAnalyzer>(outDir, fftSize, paramNames, signalType,
                                                    sweepStartHz, sweepEndHz,
                                                    sweepDurationSeconds);
}
