#include "ThdAnalyzer.h"
#include "JuceHeader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

ThdAnalyzer::ThdAnalyzer(
    const juce::File& outDir,
    int fftSize,
    double fundamentalFreq,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType)
    : fftSize(fftSize),
      fundamentalFreq(fundamentalFreq),
      paramNames(paramNames),
      outputDir(outDir),
      signalType(signalType)
{
}

ThdAnalyzer::~ThdAnalyzer() = default;

void ThdAnalyzer::applyHannWindow(std::vector<float>& buffer)
{
    const int N = static_cast<int>(buffer.size());

    for (int i = 0; i < N; ++i)
    {
        const float window = 0.5f *
            (1.0f - std::cos(
                2.0f * juce::MathConstants<float>::pi *
                static_cast<float>(i) /
                static_cast<float>(N - 1)));

        buffer[static_cast<std::size_t>(i)] *= window;
    }
}

ThdAnalyzer::HarmonicMeasurement ThdAnalyzer::computeHarmonics(
    const std::vector<std::complex<float>>& fftResult,
    double sampleRate,
    int64_t centreSample)
{
    HarmonicMeasurement measurement;
    measurement.centreSample = centreSample;

    const double binHz = sampleRate / static_cast<double>(fftSize);
    const int fundamentalBin = static_cast<int>(
        std::round(fundamentalFreq / binHz));

    if (fundamentalBin <= 0 || fundamentalBin >= fftSize / 2)
        return measurement;

    // A Hann-windowed sine spreads across neighbouring bins unless the test
    // frequency is perfectly coherent with the FFT. Integrating a small bin
    // neighbourhood gives stable harmonic amplitudes without changing the
    // existing test signal.
    auto bandMagnitude = [&fftResult, this](int centreBin)
    {
        double power = 0.0;
        for (int offset = -2; offset <= 2; ++offset)
        {
            const int bin = centreBin + offset;
            if (bin <= 0 || bin >= fftSize / 2)
                continue;

            const double magnitude = std::abs(
                fftResult[static_cast<std::size_t>(bin)]);
            power += magnitude * magnitude;
        }
        return std::sqrt(power);
    };

    const double fundamentalMagnitude = bandMagnitude(fundamentalBin);
    if (fundamentalMagnitude <= 1.0e-12)
        return measurement;

    double harmonicPowerSum = 0.0;
    for (int harmonic = 2; harmonic <= 10; ++harmonic)
    {
        const int harmonicBin = harmonic * fundamentalBin;
        if (harmonicBin >= fftSize / 2)
            break;

        const double ratio = bandMagnitude(harmonicBin) / fundamentalMagnitude;
        measurement.harmonicRatios[static_cast<std::size_t>(harmonic - 2)] = ratio;
        harmonicPowerSum += ratio * ratio;
    }

    measurement.thd = std::sqrt(harmonicPowerSum);

    // Classify energy that does not belong to the expected H2-H10 ladder.
    // This catches foldback/alias products and other out-of-ladder components
    // while the FFT buffer is still available.
    std::vector<bool> expected(static_cast<std::size_t>(fftSize / 2), false);
    auto markBand = [&](int centreBin)
    {
        for (int offset = -3; offset <= 3; ++offset)
        {
            const int bin = centreBin + offset;
            if (bin > 0 && bin < fftSize / 2) expected[static_cast<std::size_t>(bin)] = true;
        }
    };
    markBand(fundamentalBin);
    for (int harmonic = 2; harmonic <= 10; ++harmonic)
    {
        const double rawFrequency = harmonic * fundamentalFreq;
        const int directBin = harmonic * fundamentalBin;
        if (directBin < fftSize / 2) markBand(directBin);

        if (rawFrequency >= sampleRate * 0.5)
        {
            double folded = std::fmod(rawFrequency, sampleRate);
            if (folded > sampleRate * 0.5) folded = sampleRate - folded;
            const int foldedBin = static_cast<int>(std::round(folded / binHz));
            if (foldedBin > 0 && foldedBin < fftSize / 2)
            {
                measurement.aliasPowerRatio += std::pow(bandMagnitude(foldedBin) / fundamentalMagnitude, 2.0);
                ++measurement.foldedHarmonicCount;
                markBand(foldedBin);
            }
        }
    }

    double unexpectedPower = 0.0;
    double strongestUnexpected = 0.0;
    int strongestBin = 0;
    for (int bin = 1; bin < fftSize / 2; ++bin)
    {
        if (expected[static_cast<std::size_t>(bin)]) continue;
        const double magnitude = std::abs(fftResult[static_cast<std::size_t>(bin)]);
        const double ratio = magnitude / fundamentalMagnitude;
        unexpectedPower += ratio * ratio;
        if (ratio > strongestUnexpected) { strongestUnexpected = ratio; strongestBin = bin; }
    }
    measurement.unexpectedPowerRatio = std::sqrt(unexpectedPower);
    measurement.aliasPowerRatio = std::sqrt(measurement.aliasPowerRatio);
    measurement.aliasToHarmonicRatio = measurement.thd > 0.0
        ? measurement.unexpectedPowerRatio / measurement.thd : 0.0;
    measurement.strongestUnexpectedFrequencyHz = strongestBin * binHz;
    measurement.strongestUnexpectedRatio = strongestUnexpected;
    measurement.ladderIntegrity = 1.0 / (1.0 + measurement.aliasToHarmonicRatio);
    return measurement;
}

void ThdAnalyzer::processFFTWindow(
    RunThdData& data,
    int64_t centreSample)
{
    if (static_cast<int>(data.buffer.size()) < fftSize)
        return;

    applyHannWindow(data.buffer);

    juce::dsp::FFT fft(static_cast<int>(std::log2(fftSize)));
    std::vector<std::complex<float>> fftResult(
        static_cast<std::size_t>(fftSize));

    for (int i = 0; i < fftSize; ++i)
    {
        fftResult[static_cast<std::size_t>(i)] =
            std::complex<float>(
                data.buffer[static_cast<std::size_t>(i)],
                0.0f);
    }

    fft.perform(fftResult.data(), fftResult.data(), false);

    data.thdResults.push_back(
        computeHarmonics(fftResult, data.sampleRate, centreSample));
    data.buffer.clear();
}

void ThdAnalyzer::processBlock(const BlockContext& ctx)
{
    auto& data = perRunData[ctx.runId];

    if (data.buffer.empty() && data.thdResults.empty())
    {
        data.paramValues = ctx.paramNamedValues;
        data.inputGainDb = ctx.inputGainDb;
        data.sampleRate = ctx.sampleRate;
    }

    for (int i = 0; i < ctx.numSamples; ++i)
    {
        data.buffer.push_back(ctx.outL[i]);

        if (static_cast<int>(data.buffer.size()) >= fftSize)
        {
            const int64_t centreSample =
                ctx.firstSample + i - fftSize / 2;

            processFFTWindow(data, centreSample);
        }
    }
}

void ThdAnalyzer::buildResult()
{
    result.clear();
    result.analyserName = "Thd";

    result.columns.push_back("runId");
    result.columns.push_back("centreSample");
    result.columns.push_back("thd");
    for (int harmonic = 2; harmonic <= 10; ++harmonic)
        result.columns.push_back("h" + std::to_string(harmonic));
    result.columns.push_back("unexpectedPowerRatio");
    result.columns.push_back("aliasPowerRatio");
    result.columns.push_back("aliasToHarmonicRatio");
    result.columns.push_back("strongestUnexpectedFrequencyHz");
    result.columns.push_back("strongestUnexpectedRatio");
    result.columns.push_back("ladderIntegrity");
    result.columns.push_back("foldedHarmonicCount");

    for (const auto& paramName : paramNames)
        result.columns.push_back(paramName.toStdString());

    result.columns.push_back("inputGainDb");

    for (const auto& [runId, data] : perRunData)
    {
        for (const auto& measurement : data.thdResults)
        {
            std::vector<double> row;
            row.reserve(result.columns.size());

            row.push_back(static_cast<double>(runId));
            row.push_back(static_cast<double>(measurement.centreSample));
            row.push_back(measurement.thd);
            for (double ratio : measurement.harmonicRatios)
                row.push_back(ratio);
            row.push_back(measurement.unexpectedPowerRatio);
            row.push_back(measurement.aliasPowerRatio);
            row.push_back(measurement.aliasToHarmonicRatio);
            row.push_back(measurement.strongestUnexpectedFrequencyHz);
            row.push_back(measurement.strongestUnexpectedRatio);
            row.push_back(measurement.ladderIntegrity);
            row.push_back(static_cast<double>(measurement.foldedHarmonicCount));

            for (const auto& paramName : paramNames)
            {
                double value = 0.0;
                const auto valueIt = data.paramValues.find(paramName);

                if (valueIt != data.paramValues.end())
                    value = static_cast<double>(valueIt->second);

                row.push_back(value);
            }

            row.push_back(static_cast<double>(data.inputGainDb));
            result.rows.push_back(std::move(row));
        }
    }
}

void ThdAnalyzer::finish(const juce::File& outDir)
{
    juce::ignoreUnused(outDir);
    buildResult();
}

std::unique_ptr<Analyzer> createThdAnalyzer(
    const juce::File& outDir,
    int fftSize,
    double fundamentalFreq,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType)
{
    return std::make_unique<ThdAnalyzer>(
        outDir,
        fftSize,
        fundamentalFreq,
        paramNames,
        signalType);
}
