#pragma once

#include "Analyzer.h"
#include "JuceHeader.h"
#include "MeasurementResult.h"

#include <complex>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>
#include <utility>

struct LinearResponseAnalyzer : public Analyzer
{
    LinearResponseAnalyzer(
        const juce::File& outDir,
        int fftSize,
        const std::vector<juce::String>& paramNames,
        const juce::String& signalType,
        double sweepStartHz,
        double sweepEndHz,
        double sweepDurationSeconds);

    ~LinearResponseAnalyzer() override;

    void processBlock(const BlockContext& ctx) override;
    void finish(const juce::File& outDir) override;

    /**
        Returns the most recently completed in-memory linear-response dataset.

        This is populated during finish() and retained in memory for optional export.
    */
    const MeasurementDataset& getResult() const noexcept
    {
        return result;
    }

    MeasurementDataset takeResult() noexcept
    {
        return std::move(result);
    }

    const MeasurementDataset& getHarmonicFingerprintResult() const noexcept
    {
        return harmonicFingerprintResult;
    }

    MeasurementDataset takeHarmonicFingerprintResult() noexcept
    {
        return std::move(harmonicFingerprintResult);
    }

private:
    struct RunSpectrum
    {
        std::vector<double> sumInMagSq;
        std::vector<double> sumOutMagSq;
        int numAverages = 0;
        std::vector<float> inBuffer;
        std::vector<float> outBuffer;
        std::map<juce::String, float> paramValues;
        float inputGainDb = 0.0f;
        double sampleRate = 48000.0;

        static constexpr int harmonicBands = 48;
        static constexpr int maxHarmonicOrder = 7;
        std::vector<double> inputFundSin;
        std::vector<double> inputFundCos;
        std::vector<std::vector<double>> outputHarmSin;
        std::vector<std::vector<double>> outputHarmCos;
        std::vector<std::int64_t> harmonicCounts;
        double sweepPhase = 0.0;
    };

    void buildResult();
    void buildHarmonicFingerprintResult();
    void processFFTWindow(RunSpectrum& spectrum);
    void processHarmonicFingerprintSamples(RunSpectrum& spectrum, const BlockContext& ctx);
    void applyHannWindow(std::vector<float>& buffer);

    std::map<int, RunSpectrum> perRunSpectra;
    int fftSize;
    std::vector<juce::String> paramNames;
    juce::File outputDir;
    juce::String signalType;
    double sweepStartHz = 20.0;
    double sweepEndHz = 20000.0;
    double sweepDurationSeconds = 5.0;

    MeasurementDataset result;
    MeasurementDataset harmonicFingerprintResult;
};

std::unique_ptr<Analyzer> createLinearResponseAnalyzer(
    const juce::File& outDir,
    int fftSize,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType,
    double sweepStartHz,
    double sweepEndHz,
    double sweepDurationSeconds);
