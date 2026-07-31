#include "RawCsvAnalyzer.h"

#include <iostream>
#include <utility>

RawCsvAnalyzer::RawCsvAnalyzer(
    const juce::File& outDir,
    const juce::String& signalType)
    : signalType(signalType)
{
    const juce::String filename =
        "raw_" + signalType.toLowerCase() + ".csv";

    const juce::File outputFile = outDir.getChildFile(filename);

    csvFile = std::make_unique<std::ofstream>(
        outputFile.getFullPathName().toStdString());

    if (!csvFile->is_open())
    {
        std::cerr
            << "Failed to open "
            << filename.toStdString()
            << " for writing"
            << std::endl;

        csvFile.reset();
    }
}

RawCsvAnalyzer::~RawCsvAnalyzer()
{
    if (csvFile && csvFile->is_open())
        csvFile->close();
}

void RawCsvAnalyzer::processBlock(const BlockContext& ctx)
{
    const bool hasInputRight = ctx.inR != nullptr;
    const bool hasOutputRight = ctx.outR != nullptr;

    capturedInputRight = capturedInputRight || hasInputRight;
    capturedOutputRight = capturedOutputRight || hasOutputRight;

    if (csvFile && !headerWritten)
    {
        *csvFile << "runId,sample,time_sec,inL";

        if (hasInputRight)
            *csvFile << ",inR";

        *csvFile << ",outL";

        if (hasOutputRight)
            *csvFile << ",outR";

        *csvFile << "\n";
        headerWritten = true;
    }

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

        if (csvFile)
        {
            *csvFile
                << ctx.runId << ","
                << sampleIndex << ","
                << timeSec << ","
                << ctx.inL[i];

            if (hasInputRight)
                *csvFile << "," << ctx.inR[i];

            *csvFile << "," << ctx.outL[i];

            if (hasOutputRight)
                *csvFile << "," << ctx.outR[i];

            *csvFile << "\n";
        }
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

    if (csvFile && csvFile->is_open())
    {
        csvFile->close();
        csvFile.reset();
    }
}

std::unique_ptr<Analyzer> createRawCsvAnalyzer(
    const juce::File& outDir,
    const juce::String& signalType)
{
    return std::make_unique<RawCsvAnalyzer>(outDir, signalType);
}
