#pragma once

#include <cstdint>
#include <string>
#include <vector>

/**
    One named set of measurements produced by an analyser.

    This deliberately resembles the existing CSV structure:
    - columns describe the values
    - rows contain the measured numeric data

    It gives PluginDNA access to measurements in memory while allowing
    the existing CSV exporters to continue working.
*/
struct MeasurementDataset
{
    std::string analyserName;
    std::vector<std::string> columns;
    std::vector<std::vector<double>> rows;

    bool isEmpty() const noexcept
    {
        return rows.empty();
    }

    std::size_t getRowCount() const noexcept
    {
        return rows.size();
    }

    void clear()
    {
        columns.clear();
        rows.clear();
    }
};

/**
    Complete measurement output for one plugin analysis session.

    Later this will also hold:
    - plugin format and version
    - parameter state
    - sample-rate and block-size information
    - DNA fingerprints
    - optional notes and tags
*/
struct MeasurementResult
{
    std::string pluginName;
    std::string pluginManufacturer;
    std::string pluginVersion;
    std::string pluginFormat;

    double sampleRate = 0.0;
    int blockSize = 0;

    std::int64_t createdAtUnixTime = 0;

    std::vector<MeasurementDataset> datasets;

    bool isEmpty() const noexcept
    {
        return datasets.empty();
    }

    void addDataset(MeasurementDataset dataset)
    {
        datasets.push_back(std::move(dataset));
    }

    const MeasurementDataset* findDataset(
        const std::string& analyserNameToFind) const noexcept
    {
        for (const auto& dataset : datasets)
        {
            if (dataset.analyserName == analyserNameToFind)
                return &dataset;
        }

        return nullptr;
    }

    void clear()
    {
        pluginName.clear();
        pluginManufacturer.clear();
        pluginVersion.clear();
        pluginFormat.clear();

        sampleRate = 0.0;
        blockSize = 0;
        createdAtUnixTime = 0;

        datasets.clear();
    }
};
