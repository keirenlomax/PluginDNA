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

double ThdAnalyzer::computeTHD(
    const std::vector<std::complex<float>>& fftResult,
    double sampleRate)
{
    const double binHz = sampleRate / static_cast<double>(fftSize);
    const int fundamentalBin = static_cast<int>(
        std::round(fundamentalFreq / binHz));

    if (fundamentalBin <= 0 || fundamentalBin >= fftSize / 2)
        return 0.0;

    const float fundamentalMagnitude =
        std::abs(fftResult[static_cast<std::size_t>(fundamentalBin)]);

    const double fundamentalPower =
        static_cast<double>(fundamentalMagnitude * fundamentalMagnitude);

    if (fundamentalPower <= 0.0)
        return 0.0;

    double harmonicPowerSum = 0.0;
    const int maxHarmonic = std::min(
        10,
        (fftSize / 2) / fundamentalBin);

    for (int harmonic = 2; harmonic <= maxHarmonic; ++harmonic)
    {
        const int harmonicBin = harmonic * fundamentalBin;

        if (harmonicBin >= fftSize / 2)
            break;

        const float harmonicMagnitude =
            std::abs(fftResult[static_cast<std::size_t>(harmonicBin)]);

        harmonicPowerSum += static_cast<double>(
            harmonicMagnitude * harmonicMagnitude);
    }

    return std::sqrt(harmonicPowerSum / fundamentalPower);
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

    const double thd = computeTHD(fftResult, data.sampleRate);
    data.thdResults.push_back({ centreSample, thd });
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

    for (const auto& paramName : paramNames)
        result.columns.push_back(paramName.toStdString());

    result.columns.push_back("inputGainDb");

    for (const auto& [runId, data] : perRunData)
    {
        for (const auto& [centreSample, thd] : data.thdResults)
        {
            std::vector<double> row;
            row.reserve(result.columns.size());

            row.push_back(static_cast<double>(runId));
            row.push_back(static_cast<double>(centreSample));
            row.push_back(thd);

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
    buildResult();

    // Keep the existing CSV export, now written from the in-memory dataset.
    const juce::String filename =
        "grid_thd_" + signalType.toLowerCase() + ".csv";

    const juce::File csvFile = outDir.getChildFile(filename);
    std::ofstream out(csvFile.getFullPathName().toStdString());

    if (!out.is_open())
    {
        std::cerr
            << "Failed to open "
            << filename.toStdString()
            << " for writing"
            << std::endl;

        return;
    }

    for (std::size_t column = 0;
         column < result.columns.size();
         ++column)
    {
        if (column > 0)
            out << ',';

        out << result.columns[column];
    }

    out << '\n';

    for (const auto& row : result.rows)
    {
        for (std::size_t column = 0; column < row.size(); ++column)
        {
            if (column > 0)
                out << ',';

            out << row[column];
        }

        out << '\n';
    }

    std::cout
        << "Thd in-memory result: "
        << result.getRowCount()
        << " rows"
        << std::endl;
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
