#include "TransferCurveAnalyzer.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

TransferCurveAnalyzer::TransferCurveAnalyzer(
    const juce::File& outDir,
    int requestedNumBins,
    const std::vector<juce::String>& measuredParamNames,
    const juce::String& measuredSignalType)
    : numBins(requestedNumBins),
      paramNames(measuredParamNames),
      outputDir(outDir),
      signalType(measuredSignalType)
{
}

TransferCurveAnalyzer::~TransferCurveAnalyzer() = default;

int TransferCurveAnalyzer::getBinIndex(float x) const
{
    // Map x from [-1, 1] to [0, numBins - 1].
    const float normalized = (x + 1.0f) * 0.5f;
    const int bin = static_cast<int>(normalized * numBins);

    return std::clamp(bin, 0, numBins - 1);
}

float TransferCurveAnalyzer::getBinCenter(int binIndex) const
{
    // Inverse of getBinIndex().
    const float normalized =
        (static_cast<float>(binIndex) + 0.5f)
        / static_cast<float>(numBins);

    return normalized * 2.0f - 1.0f;
}

void TransferCurveAnalyzer::processBlock(const BlockContext& ctx)
{
    auto& runData = perRunBins[ctx.runId];

    // Initialise bins and run metadata on the first block for this run.
    if (runData.bins.empty())
    {
        runData.bins.resize(static_cast<std::size_t>(numBins));
        runData.paramValues = ctx.paramNamedValues;
        runData.inputGainDb = ctx.inputGainDb;
    }

    // Accumulate the input-to-output mapping.
    for (int sample = 0; sample < ctx.numSamples; ++sample)
    {
        const float x = ctx.inL[sample];
        const float y = ctx.outL[sample];
        const int binIndex = getBinIndex(x);

        auto& bin = runData.bins[static_cast<std::size_t>(binIndex)];
        bin.sumY += static_cast<double>(y);
        ++bin.count;
    }
}

void TransferCurveAnalyzer::buildResult()
{
    result.clear();
    result.analyserName = "TransferCurve";

    result.columns.push_back("runId");
    result.columns.push_back("binIndex");
    result.columns.push_back("x");
    result.columns.push_back("meanY");
    result.columns.push_back("count");

    for (const auto& paramName : paramNames)
        result.columns.push_back(paramName.toStdString());

    result.columns.push_back("inputGainDb");

    for (const auto& [runId, runData] : perRunBins)
    {
        for (int binIndex = 0; binIndex < numBins; ++binIndex)
        {
            const auto& bin =
                runData.bins[static_cast<std::size_t>(binIndex)];

            if (bin.count == 0)
                continue;

            std::vector<double> row;
            row.reserve(result.columns.size());

            row.push_back(static_cast<double>(runId));
            row.push_back(static_cast<double>(binIndex));
            row.push_back(static_cast<double>(getBinCenter(binIndex)));
            row.push_back(bin.sumY / static_cast<double>(bin.count));
            row.push_back(static_cast<double>(bin.count));

            for (const auto& paramName : paramNames)
            {
                double value = 0.0;
                const auto valueIt = runData.paramValues.find(paramName);

                if (valueIt != runData.paramValues.end())
                    value = static_cast<double>(valueIt->second);

                row.push_back(value);
            }

            row.push_back(static_cast<double>(runData.inputGainDb));
            result.rows.push_back(std::move(row));
        }
    }
}

void TransferCurveAnalyzer::finish(const juce::File& outDir)
{
    buildResult();

    // Keep the existing CSV export, now written from the in-memory dataset.
    const juce::String filename =
        "grid_transfer_curves_" + signalType.toLowerCase() + ".csv";

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
        << "TransferCurve in-memory result: "
        << result.getRowCount()
        << " rows"
        << std::endl;
}

std::unique_ptr<Analyzer> createTransferCurveAnalyzer(
    const juce::File& outDir,
    int numBins,
    const std::vector<juce::String>& paramNames,
    const juce::String& signalType)
{
    return std::make_unique<TransferCurveAnalyzer>(
        outDir,
        numBins,
        paramNames,
        signalType);
}
