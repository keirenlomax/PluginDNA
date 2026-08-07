#include "ParameterConfigComponent.h"

ParameterConfigComponent::ParameterConfigComponent(const juce::String& paramName, juce::AudioProcessorParameter* parameter) : paramName(paramName) {
    nameLabel.setText(paramName, juce::dontSendNotification);
    nameLabel.setFont(juce::Font(16.0f, juce::Font::bold));
    addAndMakeVisible(nameLabel);

    strategyLabel.setText("Strategy:", juce::dontSendNotification);
    addAndMakeVisible(strategyLabel);

    strategyCombo.addItem("Linear", 1);
    strategyCombo.addItem("ExplicitValues", 2);
    strategyCombo.addItem("Log", 3);
    strategyCombo.addItem("EdgeAndCenter", 4);
    strategyCombo.addItem("Enumerated", 5);
    strategyCombo.setSelectedId(1);
    strategyCombo.addListener(this);
    addAndMakeVisible(strategyCombo);

    minLabel.setText("Min:", juce::dontSendNotification);
    addAndMakeVisible(minLabel);
    minEditor.setText("0.0", juce::dontSendNotification);
    minEditor.addListener(this);
    addAndMakeVisible(minEditor);

    maxLabel.setText("Max:", juce::dontSendNotification);
    addAndMakeVisible(maxLabel);
    maxEditor.setText("1.0", juce::dontSendNotification);
    maxEditor.addListener(this);
    addAndMakeVisible(maxEditor);

    numBucketsLabel.setText("Buckets:", juce::dontSendNotification);
    addAndMakeVisible(numBucketsLabel);
    numBucketsEditor.setText("5", juce::dontSendNotification);
    numBucketsEditor.addListener(this);
    addAndMakeVisible(numBucketsEditor);

    valuesLabel.setText("Values (comma-separated):", juce::dontSendNotification);
    addAndMakeVisible(valuesLabel);
    valuesEditor.setText("0.0, 1.0", juce::dontSendNotification);
    valuesEditor.setMultiLine(true);
    valuesEditor.setScrollbarsShown(true);
    valuesEditor.setReturnKeyStartsNewLine(true);
    valuesEditor.addListener(this);
    addAndMakeVisible(valuesEditor);

    currentConfig.paramName = paramName;
    currentConfig.strategy = "Linear";
    currentConfig.min = 0.0f;
    currentConfig.max = 1.0f;
    currentConfig.numBuckets = 5;
    currentConfig.values = {0.0f, 1.0f};

    detectEnumeratedStates(parameter);
    updateUI();
}


namespace
{
    bool looksNumericParameterText(juce::String text)
    {
        text = text.trim().toLowerCase().removeCharacters(" ");
        const juce::StringArray suffixes { "%", "db", "hz", "khz", "mhz", "ms", "samples", "sample", "x" };
        for (const auto& suffix : suffixes)
            if (text.endsWith(suffix)) { text = text.dropLastCharacters(suffix.length()); break; }
        if (text.isEmpty()) return false;
        return text.containsOnly("+-0123456789.e");
    }
}

void ParameterConfigComponent::detectEnumeratedStates(juce::AudioProcessorParameter* parameter)
{
    if (parameter == nullptr) return;

    std::vector<float> stateValues;
    std::vector<juce::String> stateLabels;
    const int steps = parameter->getNumSteps();

    if ((parameter->isDiscrete() || steps < juce::AudioProcessor::getDefaultNumParameterSteps())
        && steps >= 2 && steps <= 64)
    {
        for (int i = 0; i < steps; ++i)
        {
            const float value = steps == 1 ? 0.0f : static_cast<float>(i) / static_cast<float>(steps - 1);
            auto label = parameter->getText(value, 128).trim();
            if (label.isEmpty()) label = "State " + juce::String(i + 1);
            stateValues.push_back(value);
            stateLabels.push_back(label);
        }
    }
    else
    {
        // Some plug-ins expose selectors as a continuous 0..1 parameter. Sample the
        // displayed text and identify contiguous named regions (e.g. Dark/Rock/Lush).
        constexpr int samples = 1025;
        juce::String previous;
        int regionStart = 0;
        struct Region { int start = 0; int end = 0; juce::String label; };
        std::vector<Region> regions;
        for (int i = 0; i < samples; ++i)
        {
            const float value = static_cast<float>(i) / static_cast<float>(samples - 1);
            auto label = parameter->getText(value, 128).trim();
            if (i == 0) { previous = label; regionStart = 0; continue; }
            if (label != previous)
            {
                regions.push_back({ regionStart, i - 1, previous });
                regionStart = i; previous = label;
            }
        }
        regions.push_back({ regionStart, samples - 1, previous });

        bool valid = regions.size() >= 2 && regions.size() <= 64;
        for (const auto& region : regions)
            valid = valid && !region.label.isEmpty() && !looksNumericParameterText(region.label);

        if (valid)
        {
            for (const auto& region : regions)
            {
                const float midpoint = 0.5f * static_cast<float>(region.start + region.end)
                                     / static_cast<float>(samples - 1);
                stateValues.push_back(midpoint);
                stateLabels.push_back(region.label);
            }
        }
    }

    if (stateValues.size() >= 2)
    {
        currentConfig.strategy = "Enumerated";
        currentConfig.values = std::move(stateValues);
        currentConfig.valueLabels = std::move(stateLabels);
        currentConfig.numBuckets = static_cast<int>(currentConfig.values.size());
        currentConfig.includePluginDefault = false; // each mode already has one representative value
        strategyCombo.setSelectedId(5, juce::dontSendNotification);

        juce::String list;
        for (size_t i = 0; i < currentConfig.valueLabels.size(); ++i)
        {
            if (i > 0) list += "\n";
            list += juce::String(static_cast<int>(i) + 1) + ". " + currentConfig.valueLabels[i]
                 + "  [" + juce::String(currentConfig.values[i], 4) + "]";
        }
        valuesEditor.setMultiLine(true);
        valuesEditor.setText(list, juce::dontSendNotification);
    }
}

ParameterConfigComponent::~ParameterConfigComponent() {}

void ParameterConfigComponent::paint(juce::Graphics& g) {
    g.fillAll(juce::Colours::darkgrey);
    g.setColour(juce::Colours::lightgrey);
    g.drawRect(getLocalBounds(), 1);
}

void ParameterConfigComponent::resized() {
    auto bounds = getLocalBounds().reduced(10);
    constexpr int rowHeight = 30;

    nameLabel.setBounds(bounds.removeFromTop(rowHeight));
    bounds.removeFromTop(5);

    auto strategyRow = bounds.removeFromTop(rowHeight);
    strategyLabel.setBounds(strategyRow.removeFromLeft(100));
    strategyCombo.setBounds(strategyRow.removeFromLeft(180));
    bounds.removeFromTop(6);

    const bool enumerated = currentConfig.strategy == "Enumerated";

    // Mode selectors are the important content for enumerated parameters, so
    // show them immediately under Strategy rather than below irrelevant range controls.
    if (enumerated)
    {
        valuesLabel.setBounds(bounds.removeFromTop(26));
        valuesEditor.setBounds(bounds.removeFromTop(180));
        return;
    }

    auto minMaxRow = bounds.removeFromTop(rowHeight);
    minLabel.setBounds(minMaxRow.removeFromLeft(50));
    minEditor.setBounds(minMaxRow.removeFromLeft(80));
    minMaxRow.removeFromLeft(20);
    maxLabel.setBounds(minMaxRow.removeFromLeft(50));
    maxEditor.setBounds(minMaxRow.removeFromLeft(80));
    bounds.removeFromTop(5);

    auto bucketsRow = bounds.removeFromTop(rowHeight);
    numBucketsLabel.setBounds(bucketsRow.removeFromLeft(100));
    numBucketsEditor.setBounds(bucketsRow.removeFromLeft(80));
    bounds.removeFromTop(5);

    valuesLabel.setBounds(bounds.removeFromTop(rowHeight));
    valuesEditor.setBounds(bounds.removeFromTop(42));
}

void ParameterConfigComponent::comboBoxChanged(juce::ComboBox* comboBox) {
    if (comboBox == &strategyCombo) {
        updateUI();
    }
}

void ParameterConfigComponent::textEditorTextChanged(juce::TextEditor& editor) {
    // Update config when text changes
    if (&editor == &minEditor) {
        currentConfig.min = (float)editor.getText().getDoubleValue();
    } else if (&editor == &maxEditor) {
        currentConfig.max = (float)editor.getText().getDoubleValue();
    } else if (&editor == &numBucketsEditor) {
        currentConfig.numBuckets = editor.getText().getIntValue();
    } else if (&editor == &valuesEditor) {
        juce::StringArray tokens;
        tokens.addTokens(editor.getText(), ",", "");
        currentConfig.values.clear();
        for (const auto& token : tokens) {
            float val = (float)token.trim().getDoubleValue();
            currentConfig.values.push_back(val);
        }
    }
}

void ParameterConfigComponent::updateUI() {
    int strategyId = strategyCombo.getSelectedId();

    bool showMinMax = (strategyId == 1 || strategyId == 3 || strategyId == 4); // Linear, Log, EdgeAndCenter
    bool showBuckets = (strategyId == 1 || strategyId == 3);                   // Linear, Log
    bool showValues = (strategyId == 2 || strategyId == 5);                    // ExplicitValues / Enumerated

    minLabel.setVisible(showMinMax);
    minEditor.setVisible(showMinMax);
    maxLabel.setVisible(showMinMax);
    maxEditor.setVisible(showMinMax);
    numBucketsLabel.setVisible(showBuckets);
    numBucketsEditor.setVisible(showBuckets);
    valuesLabel.setVisible(showValues);
    valuesEditor.setVisible(showValues);

    // Update strategy string
    if (strategyId == 1)
        currentConfig.strategy = "Linear";
    else if (strategyId == 2)
        currentConfig.strategy = "ExplicitValues";
    else if (strategyId == 3)
        currentConfig.strategy = "Log";
    else if (strategyId == 4)
        currentConfig.strategy = "EdgeAndCenter";
    else if (strategyId == 5)
        currentConfig.strategy = "Enumerated";

    const auto modeCount = static_cast<int>(currentConfig.valueLabels.size());
    valuesLabel.setText(strategyId == 5
        ? juce::String("Detected Modes (") + juce::String(modeCount) + ")"
        : "Values (comma-separated):",
        juce::dontSendNotification);
    valuesEditor.setReadOnly(strategyId == 5);

    resized();
}

ParameterBucketConfig ParameterConfigComponent::getConfig() const {
    return currentConfig;
}

void ParameterConfigComponent::setConfig(const ParameterBucketConfig& config) {
    currentConfig = config;

    if (config.strategy == "Linear")
        strategyCombo.setSelectedId(1);
    else if (config.strategy == "ExplicitValues")
        strategyCombo.setSelectedId(2);
    else if (config.strategy == "Log")
        strategyCombo.setSelectedId(3);
    else if (config.strategy == "EdgeAndCenter")
        strategyCombo.setSelectedId(4);
    else if (config.strategy == "Enumerated")
        strategyCombo.setSelectedId(5);

    minEditor.setText(juce::String(config.min), juce::dontSendNotification);
    maxEditor.setText(juce::String(config.max), juce::dontSendNotification);
    numBucketsEditor.setText(juce::String(config.numBuckets), juce::dontSendNotification);

    juce::String valuesStr;
    if (config.strategy == "Enumerated" && config.valueLabels.size() == config.values.size())
    {
        for (size_t i = 0; i < config.valueLabels.size(); ++i)
        {
            if (i > 0) valuesStr += "\n";
            valuesStr += juce::String(static_cast<int>(i) + 1) + ". " + config.valueLabels[i]
                       + "  [" + juce::String(config.values[i], 4) + "]";
        }
    }
    else
    {
        for (size_t i = 0; i < config.values.size(); ++i) {
            if (i > 0) valuesStr += ", ";
            valuesStr += juce::String(config.values[i]);
        }
    }
    valuesEditor.setText(valuesStr, juce::dontSendNotification);

    updateUI();
}
