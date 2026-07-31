#pragma once

#include "Analyzer.h"
#include "JuceHeader.h"
#include "MeasurementResult.h"

#include <map>
#include <memory>
#include <vector>
#include <utility>

struct RunStats
{
    double sumSqInL = 0.0;
    double sumSqInR = 0.0;
    double sumSqOutL = 0.0;
    double sumSqOutR = 0.0;

    float peakInL = 0.0f;
    float peakInR = 0.0f;
    float peakOutL = 0.0f;
    float peakOutR = 0.0f;

    int64_t sampleCount = 0;
};

struct RmsPeakAnalyzer : public Analyzer
{
    RmsPeakAnalyzer(const juce::File& outDir,
                    const std::vector<juce::String>& paramNames,
                    const juce::String& signalType);

    ~RmsPeakAnalyzer() override;

    void processBlock(const BlockContext& ctx) override;
    void finish(const juce::File& outDir) override;

    /**
        Returns the most recently completed in-memory RMS/peak dataset.

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
    void buildResult();

    std::map<int, RunStats> perRunStats;

    // runId -> parameter name -> normalised parameter value
    std::map<int, std::map<juce::String, float>> runParamValues;

    // runId -> input gain in dB
    std::map<int, float> runInputGainDb;

    std::vector<juce::String> paramNames;
    juce::File outputDir;
    juce::String signalType;

    MeasurementDataset result;
};

std::unique_ptr<Analyzer> createRmsPeakAnalyzer(
    const juce::File& outDir,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType);
