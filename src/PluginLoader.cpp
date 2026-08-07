#include "PluginLoader.h"
#include <iostream>

std::unique_ptr<juce::AudioPluginInstance> loadPluginInstance(const juce::File& pluginFile, double sampleRate,
                                                              int blockSize, juce::String& errorMessageOut) {
    errorMessageOut.clear();

    // VST3 and Audio Unit plug-ins may be files or macOS bundles.
    if (!pluginFile.exists()) {
        errorMessageOut = "Plugin file does not exist: " + pluginFile.getFullPathName();
        std::cerr << errorMessageOut << std::endl;
        return nullptr;
    }

    juce::AudioPluginFormatManager formatManager;
    #if JUCE_PLUGINHOST_VST3
    formatManager.addFormat(std::make_unique<juce::VST3PluginFormat>());
    #endif
    #if JUCE_PLUGINHOST_AU && JUCE_MAC
    formatManager.addFormat(std::make_unique<juce::AudioUnitPluginFormat>());
    #endif

    juce::String errorMessage;
    juce::PluginDescription description;

    // For VST3 bundles on macOS, we need to check the bundle structure
    juce::String pluginPath = pluginFile.getFullPathName();

    // If it's a .vst3 bundle, ensure we're pointing to the bundle directory
    if (pluginFile.isDirectory() && pluginPath.endsWithIgnoreCase(".vst3")) {
        // VST3 bundle is correct - use as-is
    } else if (pluginFile.isDirectory()) {
        // Might be pointing inside the bundle, try to find the .vst3 parent
        juce::File current = pluginFile;
        while (!current.isRoot() && !current.getFileName().endsWithIgnoreCase(".vst3")) {
            current = current.getParentDirectory();
        }
        if (current.getFileName().endsWithIgnoreCase(".vst3")) {
            pluginPath = current.getFullPathName();
        }
    }

    const int numFormats = formatManager.getNumFormats();

    // Try to find the plugin by scanning
    bool foundFormat = false;

    for (int i = 0; i < numFormats; ++i) {
        auto* format = formatManager.getFormat(i);
        juce::String formatName = format->getName();

        if (format->fileMightContainThisPluginType(pluginPath)) {
            foundFormat = true;
            juce::OwnedArray<juce::PluginDescription> found;
            format->findAllTypesForFile(found, pluginPath);
            if (found.size() > 0) {
                description = *found[0];
                std::cerr << "Found plugin: " << description.name << " (" << formatName << ")" << std::endl;
                break;
            } else {
                std::cerr << "Format " << formatName << " recognized file but found no plugins" << std::endl;
            }
        }
    }

    if (!foundFormat) {
        errorMessageOut = "No supported plugin format recognized for: " + pluginPath +
                          "\nSupported formats in this build: VST3 and Audio Unit.";
        std::cerr << errorMessageOut << std::endl;
        std::cerr << "Available formats: ";
        for (int i = 0; i < numFormats; ++i) {
            std::cerr << formatManager.getFormat(i)->getName();
            if (i < numFormats - 1)
                std::cerr << ", ";
        }
        std::cerr << std::endl;
        return nullptr;
    }

    if (description.name.isEmpty()) {
        errorMessageOut = "Failed to get plugin description for: " + pluginPath +
                          "\nThe file may be corrupted or not a valid plugin.";
        std::cerr << errorMessageOut << std::endl;
        return nullptr;
    }

    auto instance = formatManager.createPluginInstance(description, sampleRate, blockSize, errorMessage);

    if (instance == nullptr) {
        errorMessageOut = "Failed to create plugin instance: " + errorMessage;
        std::cerr << errorMessageOut << std::endl;
        return nullptr;
    }

    instance->prepareToPlay(sampleRate, blockSize);

    return instance;
}

std::map<juce::String, juce::AudioProcessorParameter*> buildParameterMap(juce::AudioPluginInstance& plugin,
                                                                         bool uiOnly) {
    std::map<juce::String, juce::AudioProcessorParameter*> paramMap;

    for (auto* param : plugin.getParameters()) {
        if (param == nullptr)
            continue;

        // Filter out non-UI parameters if requested
        if (uiOnly) {
            // Skip meta parameters (internal automation)
            if (param->isMetaParameter())
                continue;

            // Skip non-automatable parameters (typically internal)
            if (!param->isAutomatable())
                continue;

            // Skip MIDI CC parameters (they're not UI-exposed)
            juce::String paramName = param->getName(512);
            if (paramName.containsIgnoreCase("MIDI CC") || paramName.containsIgnoreCase("midi cc"))
                continue;
        }

        juce::String paramName = param->getName(512);
        paramName = paramName.trim().toLowerCase();
        paramMap[paramName] = param;
    }

    return paramMap;
}

void setParameterValue(juce::AudioPluginInstance& plugin,
                       const std::map<juce::String, juce::AudioProcessorParameter*>& paramMap,
                       const juce::String& paramName, float normalizedValue) {
    juce::String searchName = paramName.trim().toLowerCase();
    auto it = paramMap.find(searchName);

    if (it == paramMap.end()) {
        std::cerr << "Warning: Parameter not found: " << paramName << std::endl;
        return;
    }

    it->second->setValueNotifyingHost(normalizedValue);
}
