#include "RmsPeakAnalyzer.h"

#include <cmath>
#include <fstream>
#include <iostream>

RmsPeakAnalyzer::RmsPeakAnalyzer(
    const juce::File& outDir,
    const std::vector<juce::String>& parameterNames,
    const juce::String& measuredSignalType)
    : paramNames(parameterNames),
      outputDir(outDir),
      signalType(measuredSignalType)
{
}

RmsPeakAnalyzer::~RmsPeakAnalyzer() = default;

void RmsPeakAnalyzer::processBlock(const BlockContext& ctx)
{
    auto& stats = perRunStats[ctx.runId];

    // Store metadata on the first block of each run.
    if (runParamValues.find(ctx.runId) == runParamValues.end())
    {
        runParamValues[ctx.runId] = ctx.paramNamedValues;
        runInputGainDb[ctx.runId] = ctx.inputGainDb;
    }

    for (int i = 0; i < ctx.numSamples; ++i)
    {
        const float inL = ctx.inL[i];

        stats.sumSqInL += static_cast<double>(inL * inL);
        stats.peakInL = std::max(stats.peakInL, std::abs(inL));

        if (ctx.inR != nullptr)
        {
            const float inR = ctx.inR[i];

            stats.sumSqInR += static_cast<double>(inR * inR);
            stats.peakInR = std::max(stats.peakInR, std::abs(inR));
        }

        const float outL = ctx.outL[i];

        stats.sumSqOutL += static_cast<double>(outL * outL);
        stats.peakOutL = std::max(stats.peakOutL, std::abs(outL));

        if (ctx.outR != nullptr)
        {
            const float outR = ctx.outR[i];

            stats.sumSqOutR += static_cast<double>(outR * outR);
            stats.peakOutR = std::max(stats.peakOutR, std::abs(outR));
        }

        ++stats.sampleCount;
    }
}

void RmsPeakAnalyzer::buildResult()
{
    result.clear();
    result.analyserName = "RmsPeak";

    result.columns.push_back("runId");

    for (const auto& paramName : paramNames)
        result.columns.push_back(paramName.toStdString());

    result.columns.push_back("inputGainDb");

    result.columns.push_back("rmsInL");
    result.columns.push_back("rmsInR");
    result.columns.push_back("rmsOutL");
    result.columns.push_back("rmsOutR");

    result.columns.push_back("peakInL");
    result.columns.push_back("peakInR");
    result.columns.push_back("peakOutL");
    result.columns.push_back("peakOutR");

    for (const auto& [runId, stats] : perRunStats)
    {
        std::vector<double> row;

        row.reserve(result.columns.size());
        row.push_back(static_cast<double>(runId));

        const auto paramIt = runParamValues.find(runId);

        for (const auto& paramName : paramNames)
        {
            double value = 0.0;

            if (paramIt != runParamValues.end())
            {
                const auto valueIt = paramIt->second.find(paramName);

                if (valueIt != paramIt->second.end())
                    value = static_cast<double>(valueIt->second);
            }

            row.push_back(value);
        }

        double inputGainDb = 0.0;

        const auto gainIt = runInputGainDb.find(runId);

        if (gainIt != runInputGainDb.end())
            inputGainDb = static_cast<double>(gainIt->second);

        row.push_back(inputGainDb);

        const double rmsInL =
            stats.sampleCount > 0
                ? std::sqrt(stats.sumSqInL / stats.sampleCount)
                : 0.0;

        const double rmsInR =
            stats.sampleCount > 0
                ? std::sqrt(stats.sumSqInR / stats.sampleCount)
                : 0.0;

        const double rmsOutL =
            stats.sampleCount > 0
                ? std::sqrt(stats.sumSqOutL / stats.sampleCount)
                : 0.0;

        const double rmsOutR =
            stats.sampleCount > 0
                ? std::sqrt(stats.sumSqOutR / stats.sampleCount)
                : 0.0;

        row.push_back(rmsInL);
        row.push_back(rmsInR);
        row.push_back(rmsOutL);
        row.push_back(rmsOutR);

        row.push_back(static_cast<double>(stats.peakInL));
        row.push_back(static_cast<double>(stats.peakInR));
        row.push_back(static_cast<double>(stats.peakOutL));
        row.push_back(static_cast<double>(stats.peakOutR));

        result.rows.push_back(std::move(row));
    }
}

void RmsPeakAnalyzer::finish(const juce::File& outDir)
{
    // First create the PluginDNA in-memory result.
    buildResult();

    // Keep the existing CSV export.
    const juce::String filename =
        "grid_rms_peak_" + signalType.toLowerCase() + ".csv";

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

    // Write the in-memory dataset to CSV.
    for (std::size_t column = 0;
         column < result.columns.size();
         ++column)
    {
        if (column > 0)
            out << ",";

        out << result.columns[column];
    }

    out << "\n";

    for (const auto& row : result.rows)
    {
        for (std::size_t column = 0;
             column < row.size();
             ++column)
        {
            if (column > 0)
                out << ",";

            out << row[column];
        }

        out << "\n";
    }

    std::cout
        << "RmsPeak in-memory result: "
        << result.getRowCount()
        << " rows"
        << std::endl;
}

std::unique_ptr<Analyzer> createRmsPeakAnalyzer(
    const juce::File& outDir,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType)
{
    return std::make_unique<RmsPeakAnalyzer>(
        outDir,
        paramNames,
        signalType);
}
