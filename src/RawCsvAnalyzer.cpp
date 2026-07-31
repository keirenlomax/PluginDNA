#include "RawCsvAnalyzer.h"

#include <iostream>
#include <utility>

RawCsvAnalyzer::RawCsvAnalyzer(
    const juce::File& outDir,
    const juce::String& signalType)
    : signalType(signalType)
{
    juce::ignoreUnused(outDir);
}

RawCsvAnalyzer::~RawCsvAnalyzer() = default;

void RawCsvAnalyzer::processBlock(const BlockContext& ctx)
{
    const bool hasInputRight = ctx.inR != nullptr;
    const bool hasOutputRight = ctx.outR != nullptr;

    capturedInputRight = capturedInputRight || hasInputRight;
    capturedOutputRight = capturedOutputRight || hasOutputRight;

    for (int i = 0; i < ctx.numSamples; ++i)
    {
        const int64_t sampleIndex = ctx.firstSample + i;
        const double timeSec =
            static_cast<double>(sampleIndex) / ctx.sampleRate;

        std::vector<double> row;
        row.reserve(7);

        row.push_back(static_cast<double>(ctx.runId));
        row.push_back(static_cast<double>(sampleIndex));
        row.push_back(timeSec);
        row.push_back(static_cast<double>(ctx.inL[i]));

        if (hasInputRight)
            row.push_back(static_cast<double>(ctx.inR[i]));

        row.push_back(static_cast<double>(ctx.outL[i]));

        if (hasOutputRight)
            row.push_back(static_cast<double>(ctx.outR[i]));

        capturedRows.push_back(std::move(row));
    }
}

void RawCsvAnalyzer::buildResult()
{
    result.clear();
    result.analyserName = "RawCsv";

    result.columns.push_back("runId");
    result.columns.push_back("sample");
    result.columns.push_back("time_sec");
    result.columns.push_back("inL");

    if (capturedInputRight)
        result.columns.push_back("inR");

    result.columns.push_back("outL");

    if (capturedOutputRight)
        result.columns.push_back("outR");

    result.rows = std::move(capturedRows);
}

void RawCsvAnalyzer::finish(const juce::File& outDir)
{
    juce::ignoreUnused(outDir);
    buildResult();
}

std::unique_ptr<Analyzer> createRawCsvAnalyzer(
    const juce::File& outDir,
    const juce::String& signalType)
{
    return std::make_unique<RawCsvAnalyzer>(outDir, signalType);
}
