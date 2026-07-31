#pragma once

#include "Analyzer.h"
#include "JuceHeader.h"
#include "MeasurementResult.h"

#include <map>
#include <memory>
#include <vector>

struct TransferCurveAnalyzer : public Analyzer
{
    TransferCurveAnalyzer(
        const juce::File& outDir,
        int numBins,
        const std::vector<juce::String>& paramNames,
        const juce::String& signalType);

    ~TransferCurveAnalyzer() override;

    void processBlock(const BlockContext& ctx) override;
    void finish(const juce::File& outDir) override;

    /**
        Returns the most recently completed in-memory transfer-curve dataset.

        This is populated during finish(), alongside the existing CSV export.
    */
    const MeasurementDataset& getResult() const noexcept
    {
        return result;
    }

private:
    struct BinData
    {
        double sumY = 0.0;
        int count = 0;
    };

    struct RunBinData
    {
        std::vector<BinData> bins;
        std::map<juce::String, float> paramValues;
        float inputGainDb = 0.0f;
    };

    void buildResult();

    int getBinIndex(float x) const;
    float getBinCenter(int binIndex) const;

    std::map<int, RunBinData> perRunBins;
    int numBins;
    std::vector<juce::String> paramNames;
    juce::File outputDir;
    juce::String signalType;

    MeasurementDataset result;
};

std::unique_ptr<Analyzer> createTransferCurveAnalyzer(
    const juce::File& outDir,
    int numBins,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType);
