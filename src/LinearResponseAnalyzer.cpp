#include "LinearResponseAnalyzer.h"
#include "JuceHeader.h"
#include <cmath>
#include <fstream>
#include <iostream>

LinearResponseAnalyzer::LinearResponseAnalyzer(const juce::File& outDir, int fftSize,
                                               const std::vector<juce::String>& paramNames,
                                               const juce::String& signalType)
    : fftSize(fftSize), paramNames(paramNames), outputDir(outDir), signalType(signalType) {}

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

    // Perform FFT
    juce::dsp::FFT fft((int)std::log2(fftSize));
    std::vector<std::complex<float>> inFFT(fftSize);
    std::vector<std::complex<float>> outFFT(fftSize);

    // Copy to complex buffers
    for (int i = 0; i < fftSize; ++i) {
        inFFT[i] = std::complex<float>(spectrum.inBuffer[i], 0.0f);
        outFFT[i] = std::complex<float>(spectrum.outBuffer[i], 0.0f);
    }

    // Perform FFT
    fft.perform(inFFT.data(), outFFT.data(), false);

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

void LinearResponseAnalyzer::processBlock(const BlockContext& ctx) {
    auto& spectrum = perRunSpectra[ctx.runId];

    // Initialize on first block
    if (spectrum.sumInMagSq.empty()) {
        spectrum.paramValues = ctx.paramNamedValues;
        spectrum.inputGainDb = ctx.inputGainDb;
        spectrum.sampleRate = ctx.sampleRate;
    }

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

void LinearResponseAnalyzer::finish(const juce::File& outDir)
{
    buildResult();

    // Keep the existing CSV export, now written from the in-memory dataset.
    const juce::String filename =
        "grid_linear_response_" + signalType.toLowerCase() + ".csv";

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
        << "LinearResponse in-memory result: "
        << result.getRowCount()
        << " rows"
        << std::endl;
}

std::unique_ptr<Analyzer> createLinearResponseAnalyzer(const juce::File& outDir, int fftSize,
                                                       const std::vector<juce::String>& paramNames,
                                                       const juce::String& signalType) {
    return std::make_unique<LinearResponseAnalyzer>(outDir, fftSize, paramNames, signalType);
}
