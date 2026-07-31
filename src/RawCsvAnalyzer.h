#pragma once

#include "Analyzer.h"
#include "JuceHeader.h"
#include "MeasurementResult.h"

#include <fstream>
#include <memory>
#include <vector>

struct RawCsvAnalyzer : public Analyzer
{
    RawCsvAnalyzer(const juce::File& outDir, const juce::String& signalType);
    ~RawCsvAnalyzer() override;

    void processBlock(const BlockContext& ctx) override;
    void finish(const juce::File& outDir) override;

    /**
        Returns the most recently completed in-memory raw-sample dataset.

        This is populated during finish(), alongside the existing CSV export.
    */
    const MeasurementDataset& getResult() const noexcept
    {
        return result;
    }

private:
    void buildResult();

    std::unique_ptr<std::ofstream> csvFile;
    bool headerWritten = false;
    juce::String signalType;

    std::vector<std::vector<double>> capturedRows;
    bool capturedInputRight = false;
    bool capturedOutputRight = false;

    MeasurementDataset result;
};

std::unique_ptr<Analyzer> createRawCsvAnalyzer(
    const juce::File& outDir,
    const juce::String& signalType);
