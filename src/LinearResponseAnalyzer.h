#pragma once

#include "Analyzer.h"
#include "JuceHeader.h"
#include "MeasurementResult.h"

#include <complex>
#include <map>
#include <memory>
#include <vector>

struct LinearResponseAnalyzer : public Analyzer
{
    LinearResponseAnalyzer(
        const juce::File& outDir,
        int fftSize,
        const std::vector<juce::String>& paramNames,
        const juce::String& signalType);

    ~LinearResponseAnalyzer() override;

    void processBlock(const BlockContext& ctx) override;
    void finish(const juce::File& outDir) override;

    /**
        Returns the most recently completed in-memory linear-response dataset.

        This is populated during finish(), alongside the existing CSV export.
    */
    const MeasurementDataset& getResult() const noexcept
    {
        return result;
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
    };

    void buildResult();
    void processFFTWindow(RunSpectrum& spectrum);
    void applyHannWindow(std::vector<float>& buffer);

    std::map<int, RunSpectrum> perRunSpectra;
    int fftSize;
    std::vector<juce::String> paramNames;
    juce::File outputDir;
    juce::String signalType;

    MeasurementDataset result;
};

std::unique_ptr<Analyzer> createLinearResponseAnalyzer(
    const juce::File& outDir,
    int fftSize,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType);
