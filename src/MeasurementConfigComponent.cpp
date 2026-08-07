#include "MeasurementConfigComponent.h"

MeasurementConfigComponent::MeasurementConfigComponent() {
    signalTypeLabel.setText("Signal Type:", juce::dontSendNotification);
    addAndMakeVisible(signalTypeLabel);

    signalTypeCombo.addItem("All (Recommended)", 1);
    signalTypeCombo.addSeparator();
    signalTypeCombo.addItem("Sine", 2);
    signalTypeCombo.addItem("Noise", 3);
    signalTypeCombo.addItem("Sweep", 4);
    signalTypeCombo.setSelectedId(1);
    signalTypeCombo.addListener(this);
    addAndMakeVisible(signalTypeCombo);

    sineFreqLabel.setText("Sine Frequency (Hz):", juce::dontSendNotification);
    addAndMakeVisible(sineFreqLabel);
    sineFreqEditor.setText("1000.0", juce::dontSendNotification);
    sineFreqEditor.addListener(this);
    addAndMakeVisible(sineFreqEditor);

    sweepStartLabel.setText("Sweep Start (Hz):", juce::dontSendNotification);
    addAndMakeVisible(sweepStartLabel);
    sweepStartEditor.setText("20.0", juce::dontSendNotification);
    sweepStartEditor.addListener(this);
    addAndMakeVisible(sweepStartEditor);

    sweepEndLabel.setText("Sweep End (Hz):", juce::dontSendNotification);
    addAndMakeVisible(sweepEndLabel);
    sweepEndEditor.setText("20000.0", juce::dontSendNotification);
    sweepEndEditor.addListener(this);
    addAndMakeVisible(sweepEndEditor);

    sampleRateLabel.setText("Sample Rate:", juce::dontSendNotification);
    addAndMakeVisible(sampleRateLabel);
    sampleRateEditor.setText("96000", juce::dontSendNotification);
    sampleRateEditor.addListener(this);
    addAndMakeVisible(sampleRateEditor);

    secondsLabel.setText("Duration (seconds):", juce::dontSendNotification);
    addAndMakeVisible(secondsLabel);
    secondsEditor.setText("5.0", juce::dontSendNotification);
    secondsEditor.addListener(this);
    addAndMakeVisible(secondsEditor);

    blockSizeLabel.setText("Block Size:", juce::dontSendNotification);
    addAndMakeVisible(blockSizeLabel);
    blockSizeEditor.setText("256", juce::dontSendNotification);
    blockSizeEditor.addListener(this);
    addAndMakeVisible(blockSizeEditor);

    inputGainLabel.setText("Input Level Range (dBFS)", juce::dontSendNotification);
    addAndMakeVisible(inputGainLabel);

    inputGainStartLabel.setText("Start", juce::dontSendNotification);
    addAndMakeVisible(inputGainStartLabel);
    inputGainStartEditor.setText("-18", juce::dontSendNotification);
    inputGainStartEditor.addListener(this);
    addAndMakeVisible(inputGainStartEditor);

    inputGainEndLabel.setText("End", juce::dontSendNotification);
    addAndMakeVisible(inputGainEndLabel);
    inputGainEndEditor.setText("0", juce::dontSendNotification);
    inputGainEndEditor.addListener(this);
    addAndMakeVisible(inputGainEndEditor);

    inputGainStepLabel.setText("Step", juce::dontSendNotification);
    addAndMakeVisible(inputGainStepLabel);
    inputGainStepCombo.addItem("1 dB", 1);
    inputGainStepCombo.addItem("3 dB", 3);
    inputGainStepCombo.addItem("6 dB", 6);
    inputGainStepCombo.setSelectedId(3);
    inputGainStepCombo.addListener(this);
    addAndMakeVisible(inputGainStepCombo);

    analysisSuiteGroup.setText("Analysis Suite");
    addAndMakeVisible(analysisSuiteGroup);

    analyzersLabel.setText("Choose analysers", juce::dontSendNotification);
    addAndMakeVisible(analyzersLabel);

    allAnalyzersButton.setButtonText("All");
    noAnalyzersButton.setButtonText("None");
    for (auto* button : { &allAnalyzersButton, &noAnalyzersButton })
    {
        button->addListener(this);
        addAndMakeVisible(button);
    }

    rawCsvButton.setButtonText("Raw CSV");
    rawCsvButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(rawCsvButton);

    rmsPeakButton.setButtonText("RMS/Peak");
    rmsPeakButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(rmsPeakButton);

    transferCurveButton.setButtonText("Transfer Curve");
    transferCurveButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(transferCurveButton);

    linearResponseButton.setButtonText("Linear Response");
    linearResponseButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(linearResponseButton);

    thdButton.setButtonText("THD / Harmonics");
    thdButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(thdButton);

    interactionButton.setButtonText("Interaction / IMD");
    interactionButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(interactionButton);

    timingButton.setButtonText("Timing / Phase");
    timingButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(timingButton);

    residualButton.setButtonText("Residual");
    residualButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(residualButton);

    boundaryButton.setButtonText("Boundary");
    boundaryButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(boundaryButton);

    stereoButton.setButtonText("Stereo");
    stereoButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(stereoButton);

    summingButton.setButtonText("Summing");
    summingButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(summingButton);

    aliasButton.setButtonText("Alias Stress");
    aliasButton.setToggleState(true, juce::dontSendNotification);
    addAndMakeVisible(aliasButton);
    envelopeButton.setButtonText("Dynamic Envelope"); envelopeButton.setToggleState(true, juce::dontSendNotification); addAndMakeVisible(envelopeButton);
    hysteresisButton.setButtonText("Hysteresis / Memory"); hysteresisButton.setToggleState(true, juce::dontSendNotification); addAndMakeVisible(hysteresisButton);
    silenceButton.setButtonText("Silence / Self-Noise"); silenceButton.setToggleState(true, juce::dontSendNotification); addAndMakeVisible(silenceButton);
    truePeakButton.setButtonText("True Peak"); truePeakButton.setToggleState(true, juce::dontSendNotification); addAndMakeVisible(truePeakButton);

    analyzersLabel.setColour(juce::Label::textColourId, juce::Colours::white);
    analyzersLabel.setFont(juce::Font(juce::FontOptions(15.0f, juce::Font::bold)));
    for (auto* toggle : { &rawCsvButton, &rmsPeakButton, &transferCurveButton, &linearResponseButton,
                          &thdButton, &interactionButton, &timingButton, &residualButton,
                          &boundaryButton, &stereoButton, &summingButton, &aliasButton, &envelopeButton, &hysteresisButton, &silenceButton, &truePeakButton })
    {
        toggle->setColour(juce::ToggleButton::textColourId, juce::Colours::white);
        toggle->setColour(juce::ToggleButton::tickColourId, juce::Colours::limegreen);
        toggle->setColour(juce::ToggleButton::tickDisabledColourId, juce::Colours::grey);
        toggle->addListener(this);
    }

    updateUI();
}

MeasurementConfigComponent::~MeasurementConfigComponent() {}

void MeasurementConfigComponent::paint(juce::Graphics& g) {
    g.fillAll(juce::Colour::fromRGB(45, 48, 54));
}

void MeasurementConfigComponent::resized() {
    analysisSuiteGroup.setBounds(getLocalBounds());
    auto bounds = getLocalBounds().reduced(12, 28);
    constexpr int rowHeight = 28;
    constexpr int gap = 5;

    analyzersLabel.setBounds(bounds.removeFromTop(24));
    auto presetRow = bounds.removeFromTop(28);
    const int presetWidth = 72;
    allAnalyzersButton.setBounds(presetRow.removeFromLeft(presetWidth));
    presetRow.removeFromLeft(6);
    noAnalyzersButton.setBounds(presetRow.removeFromLeft(presetWidth));
    bounds.removeFromTop(6);

    constexpr int toggleHeight = 27;
    constexpr int toggleGap = 2;
    const int halfWidth = (bounds.getWidth() - 8) / 2;
    auto placePair = [&](juce::ToggleButton& left, juce::ToggleButton& right)
    {
        auto row = bounds.removeFromTop(toggleHeight);
        left.setBounds(row.removeFromLeft(halfWidth));
        row.removeFromLeft(8);
        right.setBounds(row);
        bounds.removeFromTop(toggleGap);
    };
    placePair(rawCsvButton, rmsPeakButton);
    placePair(transferCurveButton, linearResponseButton);
    placePair(thdButton, interactionButton);
    placePair(timingButton, residualButton);
    placePair(boundaryButton, stereoButton);
    placePair(summingButton, aliasButton);
    placePair(envelopeButton, hysteresisButton);
    placePair(silenceButton, truePeakButton);
    bounds.removeFromTop(10);

    auto placeLabelEditor = [&](juce::Label& label, juce::TextEditor& editor,
                                int labelWidth = 170, int editorWidth = 110)
    {
        auto row = bounds.removeFromTop(rowHeight);
        label.setBounds(row.removeFromLeft(labelWidth));
        row.removeFromLeft(8);
        editor.setBounds(row.removeFromLeft(editorWidth));
        bounds.removeFromTop(gap);
    };

    auto signalRow = bounds.removeFromTop(rowHeight);
    signalTypeLabel.setBounds(signalRow.removeFromLeft(120));
    signalRow.removeFromLeft(8);
    signalTypeCombo.setBounds(signalRow.removeFromLeft(190));
    bounds.removeFromTop(gap);

    if (sineFreqLabel.isVisible())
        placeLabelEditor(sineFreqLabel, sineFreqEditor);

    if (sweepStartLabel.isVisible())
    {
        auto row = bounds.removeFromTop(rowHeight);
        sweepStartLabel.setBounds(row.removeFromLeft(120));
        sweepStartEditor.setBounds(row.removeFromLeft(76));
        row.removeFromLeft(8);
        sweepEndLabel.setBounds(row.removeFromLeft(110));
        sweepEndEditor.setBounds(row.removeFromLeft(76));
        bounds.removeFromTop(gap);
    }

    auto audioRow = bounds.removeFromTop(rowHeight);
    sampleRateLabel.setBounds(audioRow.removeFromLeft(88));
    sampleRateEditor.setBounds(audioRow.removeFromLeft(78));
    audioRow.removeFromLeft(8);
    secondsLabel.setBounds(audioRow.removeFromLeft(105));
    secondsEditor.setBounds(audioRow.removeFromLeft(62));
    audioRow.removeFromLeft(8);
    blockSizeLabel.setBounds(audioRow.removeFromLeft(70));
    blockSizeEditor.setBounds(audioRow.removeFromLeft(62));
    bounds.removeFromTop(gap);

    inputGainLabel.setBounds(bounds.removeFromTop(24));
    auto gainRow = bounds.removeFromTop(rowHeight);
    inputGainStartLabel.setBounds(gainRow.removeFromLeft(42));
    inputGainStartEditor.setBounds(gainRow.removeFromLeft(58));
    gainRow.removeFromLeft(10);
    inputGainEndLabel.setBounds(gainRow.removeFromLeft(32));
    inputGainEndEditor.setBounds(gainRow.removeFromLeft(58));
    gainRow.removeFromLeft(10);
    inputGainStepLabel.setBounds(gainRow.removeFromLeft(34));
    inputGainStepCombo.setBounds(gainRow.removeFromLeft(82));
}

void MeasurementConfigComponent::comboBoxChanged(juce::ComboBox* comboBox) {
    if (comboBox == &signalTypeCombo)
        updateUI();
}

void MeasurementConfigComponent::textEditorTextChanged(juce::TextEditor& editor) {
    // Values are read when building config
}

void MeasurementConfigComponent::buttonClicked(juce::Button* button) {
    if (button == &allAnalyzersButton)
        setAllAnalyzers(true);
    else if (button == &noAnalyzersButton)
        setAllAnalyzers(false);
}

void MeasurementConfigComponent::setAllAnalyzers(bool enabled)
{
    for (auto* toggle : { &rawCsvButton, &rmsPeakButton, &transferCurveButton, &linearResponseButton,
                          &thdButton, &interactionButton, &timingButton, &residualButton,
                          &boundaryButton, &stereoButton, &summingButton, &aliasButton, &envelopeButton, &hysteresisButton, &silenceButton, &truePeakButton })
        toggle->setToggleState(enabled, juce::dontSendNotification);
}

void MeasurementConfigComponent::updateUI() {
    int signalType = signalTypeCombo.getSelectedId();
    bool showSine = (signalType == 1 || signalType == 2);
    bool showSweep = (signalType == 1 || signalType == 4);

    sineFreqLabel.setVisible(showSine);
    sineFreqEditor.setVisible(showSine);
    sweepStartLabel.setVisible(showSweep);
    sweepStartEditor.setVisible(showSweep);
    sweepEndLabel.setVisible(showSweep);
    sweepEndEditor.setVisible(showSweep);

    resized();
}

void MeasurementConfigComponent::fillConfig(Config& config) {
    int signalType = signalTypeCombo.getSelectedId();
    if (signalType == 1)
        config.signalType = "all";
    else if (signalType == 2)
        config.signalType = "sine";
    else if (signalType == 3)
        config.signalType = "noise";
    else if (signalType == 4)
        config.signalType = "sweep";

    config.sineFrequency = sineFreqEditor.getText().getDoubleValue();
    config.sweepStartHz = sweepStartEditor.getText().getDoubleValue();
    config.sweepEndHz = sweepEndEditor.getText().getDoubleValue();

    config.sampleRate = sampleRateEditor.getText().getDoubleValue();
    config.seconds = secondsEditor.getText().getDoubleValue();
    config.blockSize = blockSizeEditor.getText().getIntValue();

    // Generate one canonical input-level grid for the entire analyser suite.
    // Start/end can be edited; step is deliberately constrained to 1, 3 or 6 dB.
    double gainStart = inputGainStartEditor.getText().getDoubleValue();
    double gainEnd = inputGainEndEditor.getText().getDoubleValue();
    double gainStep = static_cast<double>(inputGainStepCombo.getSelectedId());
    if (gainStep != 1.0 && gainStep != 3.0 && gainStep != 6.0)
        gainStep = 3.0;
    if (gainStart > gainEnd)
        std::swap(gainStart, gainEnd);

    config.inputGainBucketsDb.clear();
    for (double gain = gainStart; gain <= gainEnd + 1.0e-9; gain += gainStep)
        config.inputGainBucketsDb.push_back(static_cast<float>(gain));

    // Always include the requested end point, even when the range is not an exact
    // multiple of the selected step.
    if (config.inputGainBucketsDb.empty() ||
        std::abs(static_cast<double>(config.inputGainBucketsDb.back()) - gainEnd) > 1.0e-6)
        config.inputGainBucketsDb.push_back(static_cast<float>(gainEnd));

    // Analyzers
    config.analyzers.clear();
    if (rawCsvButton.getToggleState())
        config.analyzers.push_back("RawCsv");
    if (rmsPeakButton.getToggleState())
        config.analyzers.push_back("RmsPeak");
    if (transferCurveButton.getToggleState())
        config.analyzers.push_back("TransferCurve");
    if (linearResponseButton.getToggleState())
        config.analyzers.push_back("LinearResponse");
    if (thdButton.getToggleState())
        config.analyzers.push_back("Thd");
    if (interactionButton.getToggleState())
        config.analyzers.push_back("Interaction");
    if (timingButton.getToggleState())
        config.analyzers.push_back("Timing");
    if (residualButton.getToggleState())
        config.analyzers.push_back("Residual");
    if (boundaryButton.getToggleState())
        config.analyzers.push_back("Boundary");
    if (stereoButton.getToggleState())
        config.analyzers.push_back("Stereo");
    if (summingButton.getToggleState())
        config.analyzers.push_back("Summing");
    if (aliasButton.getToggleState())
        config.analyzers.push_back("Alias");
    if (envelopeButton.getToggleState()) config.analyzers.push_back("Envelope");
    if (hysteresisButton.getToggleState()) config.analyzers.push_back("Hysteresis");
    if (silenceButton.getToggleState()) config.analyzers.push_back("Silence");
    if (truePeakButton.getToggleState()) config.analyzers.push_back("TruePeak");
}

void MeasurementConfigComponent::loadFromConfig(const Config& config) {
    if (config.signalType == "all")
        signalTypeCombo.setSelectedId(1);
    else if (config.signalType == "sine")
        signalTypeCombo.setSelectedId(2);
    else if (config.signalType == "noise")
        signalTypeCombo.setSelectedId(3);
    else if (config.signalType == "sweep")
        signalTypeCombo.setSelectedId(4);

    sineFreqEditor.setText(juce::String(config.sineFrequency), juce::dontSendNotification);
    sweepStartEditor.setText(juce::String(config.sweepStartHz), juce::dontSendNotification);
    sweepEndEditor.setText(juce::String(config.sweepEndHz), juce::dontSendNotification);

    sampleRateEditor.setText(juce::String(config.sampleRate), juce::dontSendNotification);
    secondsEditor.setText(juce::String(config.seconds), juce::dontSendNotification);
    blockSizeEditor.setText(juce::String(config.blockSize), juce::dontSendNotification);

    // Reconstruct the range controls from an existing config. Legacy/custom lists
    // still load safely; their first/last values become the visible range and the
    // closest supported step (1/3/6 dB) is selected.
    if (!config.inputGainBucketsDb.empty())
    {
        const double start = config.inputGainBucketsDb.front();
        const double end = config.inputGainBucketsDb.back();
        inputGainStartEditor.setText(juce::String(start), juce::dontSendNotification);
        inputGainEndEditor.setText(juce::String(end), juce::dontSendNotification);

        double observedStep = 3.0;
        if (config.inputGainBucketsDb.size() >= 2)
            observedStep = std::abs(static_cast<double>(config.inputGainBucketsDb[1]) -
                                    static_cast<double>(config.inputGainBucketsDb[0]));

        int selectedStep = 3;
        if (std::abs(observedStep - 1.0) <= std::abs(observedStep - 3.0) &&
            std::abs(observedStep - 1.0) <= std::abs(observedStep - 6.0))
            selectedStep = 1;
        else if (std::abs(observedStep - 6.0) < std::abs(observedStep - 3.0))
            selectedStep = 6;
        inputGainStepCombo.setSelectedId(selectedStep, juce::dontSendNotification);
    }

    auto hasAnalyzer = [&](const juce::String& name) {
        return std::find(config.analyzers.begin(), config.analyzers.end(), name) != config.analyzers.end();
    };

    rawCsvButton.setToggleState(hasAnalyzer("RawCsv"), juce::dontSendNotification);
    rmsPeakButton.setToggleState(hasAnalyzer("RmsPeak"), juce::dontSendNotification);
    transferCurveButton.setToggleState(hasAnalyzer("TransferCurve"), juce::dontSendNotification);
    linearResponseButton.setToggleState(hasAnalyzer("LinearResponse"), juce::dontSendNotification);
    thdButton.setToggleState(hasAnalyzer("Thd"), juce::dontSendNotification);
    interactionButton.setToggleState(hasAnalyzer("Interaction"), juce::dontSendNotification);
    timingButton.setToggleState(hasAnalyzer("Timing"), juce::dontSendNotification);
    residualButton.setToggleState(hasAnalyzer("Residual"), juce::dontSendNotification);
    boundaryButton.setToggleState(hasAnalyzer("Boundary"), juce::dontSendNotification);
    stereoButton.setToggleState(hasAnalyzer("Stereo"), juce::dontSendNotification);
    summingButton.setToggleState(hasAnalyzer("Summing"), juce::dontSendNotification);
    aliasButton.setToggleState(hasAnalyzer("Alias"), juce::dontSendNotification);

    updateUI();
}
