#pragma once

#include "Analyzer.h"
#include "JuceHeader.h"
#include "MeasurementResult.h"

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include <utility>

struct ThdAnalyzer : public Analyzer
{
    ThdAnalyzer(
        const juce::File& outDir,
        int fftSize,
        double fundamentalFreq,
        const std::vector<juce::String>& paramNames,
        const juce::String& signalType);

    ~ThdAnalyzer() override;

    void processBlock(const BlockContext& ctx) override;
    void finish(const juce::File& outDir) override;

    /**
        Returns the most recently completed in-memory THD dataset.

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

private:
    struct HarmonicMeasurement
    {
        int64_t centreSample = 0;
        double thd = 0.0;
        std::array<double, 9> harmonicRatios {};
    };

    struct RunThdData
    {
        std::vector<float> buffer;
        std::vector<HarmonicMeasurement> thdResults;
        std::map<juce::String, float> paramValues;
        float inputGainDb = 0.0f;
        double sampleRate = 48000.0;
    };

    void buildResult();
    void processFFTWindow(RunThdData& data, int64_t centreSample);
    void applyHannWindow(std::vector<float>& buffer);
    HarmonicMeasurement computeHarmonics(
        const std::vector<std::complex<float>>& fftResult,
        double sampleRate,
        int64_t centreSample);

    std::map<int, RunThdData> perRunData;
    int fftSize;
    double fundamentalFreq;
    std::vector<juce::String> paramNames;
    juce::File outputDir;
    juce::String signalType;

    MeasurementDataset result;
};

std::unique_ptr<Analyzer> createThdAnalyzer(
    const juce::File& outDir,
    int fftSize,
    double fundamentalFreq,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType);
