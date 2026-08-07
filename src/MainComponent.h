#pragma once

#include "Config.h"
#include "JuceHeader.h"
#include "MeasurementEngine.h"
#include "MeasurementResult.h"
#include "PluginLoader.h"
#include <memory>
#include <mutex>
#include <atomic>
#include <vector>

class ParameterConfigComponent;
class MeasurementConfigComponent;

class MainComponent : public juce::Component,
                      public juce::Button::Listener,
                      public juce::TextEditor::Listener,
                      public juce::ListBoxModel {
public:
    MainComponent();
    ~MainComponent() override;

    void paint(juce::Graphics& g) override;
    void resized() override;
    void buttonClicked(juce::Button* button) override;
    void textEditorTextChanged(juce::TextEditor& editor) override;

    // ListBoxModel
    int getNumRows() override;
    void paintListBoxItem(int rowNumber, juce::Graphics& g, int width, int height, bool rowIsSelected) override;
    void listBoxItemClicked(int row, const juce::MouseEvent& e) override;

private:
    void loadPlugin();
    void refreshAdditionalSerialStageVisibility();
    void loadAdditionalSerialStage(int index);
    void loadSerialPlugin();
    void scanSerialPluginParameters();
    void saveCurrentParameterPanelState();
    void loadParameterPanelState(bool serial);
    void openPluginEditor(juce::AudioPluginInstance* plugin, std::unique_ptr<juce::DocumentWindow>& window, const juce::String& title);
    void scanPluginParameters();
    void updateParameterList();
    void runMeasurement();
    Config buildConfigFromUI();
    void showError(const juce::String& message);
    void clearResultsSummary();
    juce::String buildSummaryText() const;
    void exportMeasurementData();
    void exportEvidencePack(const juce::File& overrideOutputDirectory = {});

    struct PendingExportDataset
    {
        juce::String filename;
        MeasurementDataset dataset;
    };

    struct AdditionalSerialStage
    {
        int stageNumber = 0;
        std::unique_ptr<juce::Label> label;
        std::unique_ptr<juce::TextEditor> pathEditor;
        std::unique_ptr<juce::TextButton> browseButton;
        std::unique_ptr<juce::TextButton> loadButton;
        std::unique_ptr<juce::TextButton> openButton;
        std::unique_ptr<juce::TextButton> removeButton;
        std::unique_ptr<juce::TextButton> parameterButton;
        std::unique_ptr<juce::TextButton> targetButton;
        std::unique_ptr<juce::Label> infoLabel;
        std::vector<juce::String> availableParameters;
        std::vector<bool> selectedParameters;
        std::map<juce::String, juce::AudioProcessorParameter*> parameterMap;
        std::vector<ParameterBucketConfig> savedBuckets;
        std::unique_ptr<juce::AudioPluginInstance> pluginInstance;
        std::unique_ptr<juce::DocumentWindow> editorWindow;
    };

    // Plugin selection
    juce::Label buildVersionLabel;
    juce::Label pluginPathLabel;
    juce::TextEditor pluginPathEditor;
    juce::TextButton browseButton;
    juce::TextButton loadPluginButton;
    juce::TextButton openPluginButton;

    // Optional second processor for serial stage analysis
    juce::Label serialPluginLabel;
    juce::TextEditor serialPluginPathEditor;
    juce::TextButton serialBrowseButton;
    juce::TextButton serialLoadButton;
    juce::TextButton serialOpenButton;
    juce::Label serialPluginInfoLabel;
    std::vector<std::unique_ptr<AdditionalSerialStage>> additionalSerialStages;
    static constexpr int maxSerialStages = 8;

    // Plugin info
    juce::Label pluginInfoLabel;

    // Parameter selection and configuration
    juce::GroupComponent parameterScanGroup;
    juce::TextButton plugin1ParameterButton;
    juce::TextButton plugin2ParameterButton;
    juce::Label parameterTargetLabel;
    juce::Label parametersLabel;
    juce::TextEditor parameterSearchEditor;
    juce::ListBox parameterListBox;
    std::vector<juce::String> availableParameters;
    std::vector<bool> selectedParameters;
    std::vector<int> filteredParameterIndices;
    bool editingSerialParameters = false;
    int editingParameterStage = 1;
    std::vector<juce::String> primaryAvailableParameters;
    std::vector<juce::String> serialAvailableParameters;
    std::vector<bool> primarySelectedParameters;
    std::vector<bool> serialSelectedParameters;
    std::map<juce::String, juce::AudioProcessorParameter*> primaryParameterMap;
    std::map<juce::String, juce::AudioProcessorParameter*> serialParameterMap;
    std::vector<ParameterBucketConfig> primarySavedBuckets;
    std::vector<ParameterBucketConfig> serialSavedBuckets;
    juce::TextButton selectAllButton;
    juce::TextButton deselectAllButton;
    juce::ToggleButton noParameterScanButton;

    // Parameter configuration panel
    juce::Viewport parameterConfigViewport;
    juce::Component parameterConfigContainer;
    std::vector<std::unique_ptr<ParameterConfigComponent>> parameterConfigComponents;

    // Measurement configuration
    std::unique_ptr<MeasurementConfigComponent> measurementConfig;

    // Output
    juce::Label outputPathLabel;
    juce::TextEditor outputPathEditor;
    juce::TextButton browseOutputButton;
    juce::TextButton exportDataButton;
    juce::TextButton exportEvidenceButton;

    // Human-readable measurement summary
    juce::GroupComponent resultsSummaryGroup;
    juce::TextButton copySummaryButton;
    juce::Label identityRoleLabel;
    juce::Label identityWhyLabel;
    juce::Label identityWatchLabel;
    juce::Label levelSummaryLabel;
    juce::Label peakSummaryLabel;
    juce::Label dynamicsSummaryLabel;
    juce::Label harmonicSummaryLabel;
    juce::Label harmonicBalanceLabel;
    juce::Label harmonicGrowthLabel;
    juce::Label toneBassLabel;
    juce::Label toneMidLabel;
    juce::Label toneTrebleLabel;
    juce::Label toneLargestLabel;
    juce::Label behaviourSummaryLabel;
    juce::Label behaviourChangesLabel;
    juce::Label operatingRangeLabel;
    juce::Label characterStartsLabel;
    juce::Label underStressLabel;
    juce::Label nonlinearityLabel;
    juce::Label curveBehaviourLabel;
    juce::Label waveformShapeLabel;
    juce::Label waveformKneeLabel;
    juce::Label waveformStyleLabel;
    juce::Label temporalResponseLabel;
    juce::Label temporalAttackLabel;
    juce::Label temporalAfterLabel;

    // Run button
    juce::TextButton runMeasurementButton;
    juce::TextButton stopMeasurementButton;
    std::atomic<bool> measurementCancelRequested { false };

    // Progress
    juce::Label progressLabel;
    juce::ProgressBar progressBar;
    double progress = 0.0;

    // Plugin instance
    std::unique_ptr<juce::AudioPluginInstance> pluginInstance;
    std::unique_ptr<juce::AudioPluginInstance> serialPluginInstance;
    std::unique_ptr<juce::DocumentWindow> pluginEditorWindow;
    std::unique_ptr<juce::DocumentWindow> serialPluginEditorWindow;
    std::map<juce::String, juce::AudioProcessorParameter*> parameterMap;

    std::vector<PendingExportDataset> pendingExportDatasets;
    Config pendingEvidenceConfig;
    juce::String pendingEvidenceSummary;
    bool hasPendingEvidencePack = false;
    std::mutex pendingExportMutex;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainComponent)
};
