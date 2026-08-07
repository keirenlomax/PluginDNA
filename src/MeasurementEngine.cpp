#include "MeasurementEngine.h"
#include "BucketSpec.h"
#include "LinearResponseAnalyzer.h"
#include "InteractionAnalyzer.h"
#include "TimingAnalyzer.h"
#include "ResidualAnalyzer.h"
#include "BoundaryAnalyzer.h"
#include "StereoAnalyzer.h"
#include "SummingAnalyzer.h"
#include "AliasAnalyzer.h"
#include "PluginLoader.h"
#include "RawCsvAnalyzer.h"
#include "RmsPeakAnalyzer.h"
#include "ThdAnalyzer.h"
#include "TransferCurveAnalyzer.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>

std::vector<RunConfig> buildRunGrid(const Config& config, const std::vector<juce::String>& paramNames) {
    std::vector<RunConfig> runs;

    // Convert ParameterBucketConfig to BucketSpec and generate values
    std::vector<std::pair<juce::String, std::vector<float>>> paramValueLists;

    for (const auto& bucketConfig : config.parameterBuckets) {
        BucketSpec spec;
        spec.paramName = bucketConfig.paramName;
        spec.strategy = BucketSpec::strategyFromString(bucketConfig.strategy);
        spec.min = bucketConfig.min;
        spec.max = bucketConfig.max;
        spec.numBuckets = bucketConfig.numBuckets;
        spec.values = bucketConfig.values;

        auto values = spec.generateValues();
        if (bucketConfig.includePluginDefault && !bucketConfig.strategy.equalsIgnoreCase("Enumerated"))
        {
            const auto defaultValue = juce::jlimit(0.0f, 1.0f, bucketConfig.pluginDefaultValue);
            const auto alreadyIncluded = std::any_of(values.begin(), values.end(), [defaultValue](float value)
            {
                return std::abs(value - defaultValue) < 1.0e-5f;
            });
            if (!alreadyIncluded)
                values.push_back(defaultValue);
            std::sort(values.begin(), values.end());
            values.erase(std::unique(values.begin(), values.end(), [](float a, float b)
            {
                return std::abs(a - b) < 1.0e-5f;
            }), values.end());
        }
        paramValueLists.push_back({bucketConfig.paramName, values});
    }

    // Build Cartesian product of parameter values and input gain buckets
    int runId = 0;

    // Helper function to generate combinations recursively
    std::function<void(int, std::map<juce::String, float>)> generateCombinations;
    generateCombinations = [&](int paramIndex, std::map<juce::String, float> currentParams) {
        if (paramIndex >= (int)paramValueLists.size()) {
            // All parameters set, now combine with input gain buckets
            for (float inputGainDb : config.inputGainBucketsDb) {
                RunConfig run;
                run.runId = runId++;
                run.paramValues = currentParams;
                run.inputGainDb = inputGainDb;
                runs.push_back(run);
            }
            return;
        }

        // Try all values for current parameter
        const auto& [paramName, values] = paramValueLists[paramIndex];
        for (float value : values) {
            auto newParams = currentParams;
            newParams[paramName] = value;
            generateCombinations(paramIndex + 1, newParams);
        }
    };

    generateCombinations(0, {});
    return runs;
}

std::vector<std::unique_ptr<Analyzer>> createAnalyzers(const Config& config, const juce::File& outDir,
                                                       const std::vector<juce::String>& paramNames) {
    std::vector<std::unique_ptr<Analyzer>> analyzers;

    if (config.signalType.equalsIgnoreCase("interaction"))
    {
        analyzers.push_back(createInteractionAnalyzer(config.sampleRate, config.seconds, paramNames));
        return analyzers;
    }
    if (config.signalType.equalsIgnoreCase("timing"))
    {
        analyzers.push_back(createTimingAnalyzer(config.sampleRate, 8192, paramNames));
        return analyzers;
    }
    if (config.signalType.equalsIgnoreCase("residual"))
    {
        analyzers.push_back(createResidualAnalyzer(config.sampleRate, paramNames));
        return analyzers;
    }
    if (config.signalType.equalsIgnoreCase("boundary"))
    {
        analyzers.push_back(createBoundaryAnalyzer(config.sampleRate, config.seconds, paramNames));
        return analyzers;
    }
    if (config.signalType.equalsIgnoreCase("stereo"))
    {
        analyzers.push_back(createStereoAnalyzer(static_cast<int64_t>(config.seconds * config.sampleRate / 4.0), paramNames));
        return analyzers;
    }
    if (config.signalType.equalsIgnoreCase("summing"))
    {
        analyzers.push_back(createSummingAnalyzer(static_cast<int64_t>(config.seconds * config.sampleRate / 3.0), paramNames));
        return analyzers;
    }
    if (config.signalType.equalsIgnoreCase("alias"))
    {
        analyzers.push_back(createAliasAnalyzer(config.sampleRate, config.seconds, paramNames));
        return analyzers;
    }

    for (const auto& analyzerName : config.analyzers) {
        if (analyzerName.equalsIgnoreCase("RawCsv")) {
            analyzers.push_back(createRawCsvAnalyzer(outDir, config.signalType));
        } else if (analyzerName.equalsIgnoreCase("RmsPeak")) {
            analyzers.push_back(createRmsPeakAnalyzer(outDir, paramNames, config.signalType));
        } else if (analyzerName.equalsIgnoreCase("TransferCurve")) {
            analyzers.push_back(createTransferCurveAnalyzer(outDir, 512, paramNames, config.signalType));
        } else if (analyzerName.equalsIgnoreCase("LinearResponse")) {
            if (config.signalType.equalsIgnoreCase("noise") || config.signalType.equalsIgnoreCase("sweep")) {
                analyzers.push_back(createLinearResponseAnalyzer(outDir, 4096, paramNames, config.signalType,
                                                                config.sweepStartHz, config.sweepEndHz,
                                                                config.seconds));
            } else {
                std::cerr << "Warning: LinearResponse analyzer requires noise or sweep signal type" << std::endl;
            }
        } else if (analyzerName.equalsIgnoreCase("Thd")) {
            if (config.signalType.equalsIgnoreCase("sine")) {
                analyzers.push_back(
                    createThdAnalyzer(outDir, 2048, config.sineFrequency, paramNames, config.signalType));
            } else {
                std::cerr << "Warning: Thd analyzer requires sine signal type" << std::endl;
            }
        } else {
            std::cerr << "Warning: Unknown analyzer: " << analyzerName << std::endl;
        }
    }

    return analyzers;
}

void runMeasurementGrid(juce::AudioPluginInstance& plugin, double sampleRate, int blockSize, int64_t totalSamples,
                        const std::vector<RunConfig>& runs, const std::vector<std::unique_ptr<Analyzer>>& analyzers,
                        const Config& config, const juce::File& outDir, std::function<void(int)> progressCallback,
                        const std::vector<juce::AudioPluginInstance*>& serialPlugins) {
    auto paramMap = buildParameterMap(plugin, false); // Plugin 1 parameters
    std::vector<std::map<juce::String, juce::AudioProcessorParameter*>> serialParamMaps;
    serialParamMaps.reserve(serialPlugins.size());
    for (auto* serialPlugin : serialPlugins)
        serialParamMaps.push_back(serialPlugin != nullptr ? buildParameterMap(*serialPlugin, false)
                                                          : std::map<juce::String, juce::AudioProcessorParameter*>{});

    // Build parameter name list in order
    std::vector<juce::String> paramNames;
    for (const auto& bucket : config.parameterBuckets) {
        paramNames.push_back(bucket.paramName);
    }

    juce::AudioBuffer<float> inputBuffer(2, blockSize);
    juce::AudioBuffer<float> outputBuffer(2, blockSize);
    juce::MidiBuffer midiBuffer;

    int runCount = 0;
    for (const auto& run : runs) {
        runCount++;
        if (progressCallback) {
            progressCallback(run.runId);
        }
        // Set plugin parameters
        for (const auto& [paramName, value] : run.paramValues) {
            if (paramName.startsWith("P1::")) {
                setParameterValue(plugin, paramMap, paramName.substring(4), value);
            } else if (paramName.startsWith("P") && paramName.contains("::")) {
                const int delimiter = paramName.indexOf("::");
                const int stageNumber = paramName.substring(1, delimiter).getIntValue();
                const int serialIndex = stageNumber - 2;
                if (serialIndex >= 0 && serialIndex < (int)serialPlugins.size() && serialPlugins[(size_t)serialIndex] != nullptr)
                    setParameterValue(*serialPlugins[(size_t)serialIndex], serialParamMaps[(size_t)serialIndex],
                                      paramName.substring(delimiter + 2), value);
            } else {
                setParameterValue(plugin, paramMap, paramName, value);
            }
        }

        // Convert input gain from dB to linear amplitude
        float inputGainLinear = std::pow(10.0f, run.inputGainDb / 20.0f);

        // Create signal generator
        std::unique_ptr<SineGenerator> sineGen;
        std::unique_ptr<NoiseGenerator> noiseGen;
        std::unique_ptr<SweepGenerator> sweepGen;
        std::unique_ptr<InteractionGenerator> interactionGen;
        std::unique_ptr<DeltaGenerator> deltaGen;
        std::unique_ptr<BoundaryGenerator> boundaryGen;
        std::unique_ptr<StereoTestGenerator> stereoGen;
        std::unique_ptr<SummingTestGenerator> summingGen;
        std::unique_ptr<AliasStressGenerator> aliasGen;

        if (config.signalType.equalsIgnoreCase("sine")) {
            sineGen = std::make_unique<SineGenerator>();
            sineGen->sampleRate = sampleRate;
            sineGen->frequency = config.sineFrequency;
            sineGen->amplitude = inputGainLinear;
        } else if (config.signalType.equalsIgnoreCase("noise")) {
            noiseGen = std::make_unique<NoiseGenerator>();
            noiseGen->amplitude = inputGainLinear;
        } else if (config.signalType.equalsIgnoreCase("sweep")) {
            sweepGen = std::make_unique<SweepGenerator>();
            sweepGen->sampleRate = sampleRate;
            sweepGen->startHz = config.sweepStartHz;
            sweepGen->endHz = config.sweepEndHz;
            sweepGen->duration = config.seconds;
            sweepGen->amplitude = inputGainLinear;
            sweepGen->reset();
        } else if (config.signalType.equalsIgnoreCase("interaction")) {
            interactionGen = std::make_unique<InteractionGenerator>();
            interactionGen->sampleRate = sampleRate;
            interactionGen->duration = config.seconds;
            interactionGen->amplitude = inputGainLinear;
            interactionGen->reset();
        } else if (config.signalType.equalsIgnoreCase("timing")) {
            deltaGen = std::make_unique<DeltaGenerator>();
            deltaGen->impulseSample = 8192;
            deltaGen->amplitude = inputGainLinear;
            deltaGen->reset();
        } else if (config.signalType.equalsIgnoreCase("residual")) {
            noiseGen = std::make_unique<NoiseGenerator>();
            noiseGen->amplitude = inputGainLinear;
        } else if (config.signalType.equalsIgnoreCase("boundary")) {
            boundaryGen = std::make_unique<BoundaryGenerator>();
            boundaryGen->sampleRate = sampleRate;
            boundaryGen->duration = config.seconds;
            boundaryGen->amplitude = inputGainLinear;
            boundaryGen->reset();
        } else if (config.signalType.equalsIgnoreCase("stereo")) {
            stereoGen = std::make_unique<StereoTestGenerator>();
            stereoGen->sampleRate = sampleRate; stereoGen->duration = config.seconds; stereoGen->amplitude = inputGainLinear; stereoGen->reset();
        } else if (config.signalType.equalsIgnoreCase("summing")) {
            summingGen = std::make_unique<SummingTestGenerator>();
            summingGen->sampleRate = sampleRate; summingGen->duration = config.seconds; summingGen->amplitude = inputGainLinear; summingGen->reset();
        } else if (config.signalType.equalsIgnoreCase("alias")) {
            aliasGen = std::make_unique<AliasStressGenerator>();
            aliasGen->sampleRate = sampleRate; aliasGen->duration = config.seconds; aliasGen->amplitude = inputGainLinear; aliasGen->reset();
        }

        // Process samples
        int64_t currentSample = 0;
        while (currentSample < totalSamples) {
            int numThisBlock = (int)std::min((int64_t)blockSize, totalSamples - currentSample);
            // Clear buffers
            inputBuffer.clear();
            outputBuffer.clear();

            // Fill input with test signal
            if (sineGen) {
                sineGen->fillBlock(inputBuffer, numThisBlock);
            } else if (noiseGen) {
                noiseGen->fillBlock(inputBuffer, numThisBlock);
            } else if (sweepGen) {
                sweepGen->fillBlock(inputBuffer, numThisBlock);
            } else if (interactionGen) {
                interactionGen->fillBlock(inputBuffer, numThisBlock);
            } else if (deltaGen) {
                deltaGen->fillBlock(inputBuffer, numThisBlock);
            } else if (boundaryGen) {
                boundaryGen->fillBlock(inputBuffer, numThisBlock);
            } else if (stereoGen) {
                stereoGen->fillBlock(inputBuffer, numThisBlock);
            } else if (summingGen) {
                summingGen->fillBlock(inputBuffer, numThisBlock);
            } else if (aliasGen) {
                aliasGen->fillBlock(inputBuffer, numThisBlock);
            }

            // Copy input to output buffer (processBlock works in-place)
            outputBuffer.makeCopyOf(inputBuffer);

            // Process through plugin (modifies outputBuffer in-place)
            plugin.processBlock(outputBuffer, midiBuffer);
            for (auto* serialPlugin : serialPlugins)
                if (serialPlugin != nullptr)
                    serialPlugin->processBlock(outputBuffer, midiBuffer);

            // Build BlockContext
            BlockContext ctx;
            ctx.firstSample = currentSample;
            ctx.sampleRate = sampleRate;
            ctx.numSamples = numThisBlock;
            ctx.inL = inputBuffer.getReadPointer(0);
            ctx.inR = inputBuffer.getNumChannels() > 1 ? inputBuffer.getReadPointer(1) : nullptr;
            ctx.outL = outputBuffer.getReadPointer(0);
            ctx.outR = outputBuffer.getNumChannels() > 1 ? outputBuffer.getReadPointer(1) : nullptr;
            ctx.runId = run.runId;
            ctx.paramNamedValues = run.paramValues;
            ctx.inputGainDb = run.inputGainDb;

            // Build params vector in fixed order
            for (const auto& paramName : paramNames) {
                float value = 0.0f;
                auto it = run.paramValues.find(paramName);
                if (it != run.paramValues.end())
                    value = it->second;
                ctx.params.push_back(value);
            }

            // Process through analyzers
            for (auto& analyzer : analyzers) {
                analyzer->processBlock(ctx);
            }

            currentSample += numThisBlock;
        }
    }

}
