#pragma once

#include "Analyzer.h"
#include "JuceHeader.h"
#include "MeasurementResult.h"

#include <complex>
#include <map>
#include <memory>
#include <vector>

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

        This is populated during finish(), alongside the existing CSV export.
    */
    const MeasurementDataset& getResult() const noexcept
    {
        return result;
    }

private:
    struct RunThdData
    {
        std::vector<float> buffer;
        std::vector<std::pair<int64_t, double>> thdResults;
        std::map<juce::String, float> paramValues;
        float inputGainDb = 0.0f;
        double sampleRate = 48000.0;
    };

    void buildResult();
    void processFFTWindow(RunThdData& data, int64_t centreSample);
    void applyHannWindow(std::vector<float>& buffer);
    double computeTHD(
        const std::vector<std::complex<float>>& fftResult,
        double sampleRate);

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
