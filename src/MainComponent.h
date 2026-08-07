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



    struct SpectralResponseComponent : public juce::Component
    {
        struct Trace
        {
            double inputGainDb = 0.0;
            double parameterValue = 0.0;
            juce::String label;
            std::vector<double> frequencyHz;
            std::vector<double> magnitudeDb;
        };

        void setData(juce::String parameterNameIn,
                     double displayedParameterValueIn,
                     std::vector<Trace> tracesIn)
        {
            parameterName = std::move(parameterNameIn);
            displayedParameterValue = displayedParameterValueIn;
            traces = std::move(tracesIn);
            repaint();
        }

        void clear()
        {
            parameterName.clear();
            traces.clear();
            repaint();
        }

        void paint(juce::Graphics& g) override
        {
            const auto bg = getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId);
            const auto text = getLookAndFeel().findColour(juce::Label::textColourId);
            const auto accent = getLookAndFeel().findColour(juce::TextButton::buttonOnColourId);
            const auto muted = text.withAlpha(0.55f);

            g.fillAll(bg);
            auto bounds = getLocalBounds().toFloat().reduced(10.0f);

            g.setColour(text);
            g.setFont(juce::Font(16.0f, juce::Font::bold));
            g.drawText("SPECTRAL DNA", bounds.removeFromTop(24.0f), juce::Justification::centredLeft);

            g.setFont(juce::Font(13.0f));
            g.setColour(muted);
            juce::String subtitle = traces.empty()
                ? "Run a sweep-based parameter scan to see reference vs processed response"
                : parameterName.fromFirstOccurrenceOf("::", false, false)
                    + " " + juce::String(displayedParameterValue, 2)
                    + "   ·   overlays = input level";
            g.drawText(subtitle, bounds.removeFromTop(22.0f), juce::Justification::centredLeft);

            if (traces.empty())
                return;

            bounds.removeFromTop(8.0f);
            auto legendArea = bounds.removeFromBottom(44.0f);
            auto plot = bounds.reduced(54.0f, 24.0f);
            lastPlot = plot;

            // Fixed display range: enough to reveal tonal movement without tiny changes disappearing.
            constexpr double yMin = -12.0;
            constexpr double yMax = 12.0;
            const auto xFor = [&](double hz)
            {
                const double x = (std::log10(juce::jlimit(20.0, 20000.0, hz)) - std::log10(20.0))
                               / (std::log10(20000.0) - std::log10(20.0));
                return plot.getX() + (float)x * plot.getWidth();
            };
            const auto yFor = [&](double db)
            {
                const double n = (juce::jlimit(yMin, yMax, db) - yMin) / (yMax - yMin);
                return plot.getBottom() - (float)n * plot.getHeight();
            };

            // Grid.
            const double freqs[] = {20,50,100,200,500,1000,2000,5000,10000,20000};
            const double dbs[] = {-12,-6,0,6,12};

            g.setFont(juce::Font(11.0f));
            for (double f : freqs)
            {
                const float x = xFor(f);
                g.setColour(text.withAlpha((f==1000.0)?0.16f:0.08f));
                g.drawVerticalLine((int)x, plot.getY(), plot.getBottom());
                g.setColour(muted);
                juce::String lab;
                if (f >= 1000.0) lab = juce::String(f/1000.0, f==1000.0||f==2000.0||f==5000.0?0:0) + "k";
                else lab = juce::String((int)f);
                g.drawText(lab, juce::Rectangle<float>(x-24.0f, plot.getBottom()+4.0f, 48.0f, 18.0f),
                           juce::Justification::centred);
            }
            for (double db : dbs)
            {
                const float y = yFor(db);
                g.setColour(text.withAlpha(db==0.0?0.26f:0.08f));
                g.drawHorizontalLine((int)y, plot.getX(), plot.getRight());
                g.setColour(muted);
                g.drawText((db>0?"+":"") + juce::String((int)db),
                           juce::Rectangle<float>(plot.getX()-48.0f,y-9.0f,42.0f,18.0f),
                           juce::Justification::centredRight);
            }

            // Reference = perfectly flat input-to-output response.
            g.setColour(text.withAlpha(0.75f));
            juce::Path reference;
            reference.startNewSubPath(plot.getX(), yFor(0.0));
            reference.lineTo(plot.getRight(), yFor(0.0));
            g.strokePath(reference, juce::PathStrokeType(1.8f));

            // Processed overlays: lowest input is faintest; highest is strongest.
            for (size_t ti=0; ti<traces.size(); ++ti)
            {
                const auto& t = traces[ti];
                if (t.frequencyHz.size()<2 || t.magnitudeDb.size()!=t.frequencyHz.size())
                    continue;

                const float alpha = 0.25f + 0.70f * ((float)ti / std::max(1.0f,(float)traces.size()-1.0f));
                juce::Path path;
                bool started=false;
                for (size_t i=0; i<t.frequencyHz.size(); ++i)
                {
                    const double f=t.frequencyHz[i], m=t.magnitudeDb[i];
                    if (!std::isfinite(f)||!std::isfinite(m)||f<20.0||f>20000.0) continue;
                    const float x=xFor(f), y=yFor(m);
                    if (!started){path.startNewSubPath(x,y);started=true;}
                    else path.lineTo(x,y);
                }
                g.setColour(accent.withAlpha(alpha));
                g.strokePath(path, juce::PathStrokeType(ti+1==traces.size()?2.3f:1.3f,
                                                        juce::PathStrokeType::curved,
                                                        juce::PathStrokeType::rounded));
            }

            // Labels.
            g.setColour(text);
            g.setFont(juce::Font(11.5f));
            g.drawText("reference 0 dB", juce::Rectangle<float>(plot.getRight()-100.0f,yFor(0.0)-18.0f,96.0f,16.0f),
                       juce::Justification::centredRight);

            // Legend across bottom.
            float lx=legendArea.getX();
            const float ly=legendArea.getY()+8.0f;
            g.setColour(muted);
            g.drawText("INPUT", juce::Rectangle<float>(lx,ly,42.0f,18.0f), juce::Justification::centredLeft);
            lx += 48.0f;
            for (size_t ti=0; ti<traces.size(); ++ti)
            {
                const float alpha = 0.25f + 0.70f * ((float)ti / std::max(1.0f,(float)traces.size()-1.0f));
                g.setColour(accent.withAlpha(alpha));
                g.drawLine(lx,ly+9.0f,lx+22.0f,ly+9.0f,ti+1==traces.size()?2.3f:1.3f);
                g.setColour(text.withAlpha(0.8f));
                const auto label=juce::String(traces[ti].inputGainDb,0)+" dB";
                g.drawText(label,juce::Rectangle<float>(lx+26.0f,ly,52.0f,18.0f),juce::Justification::centredLeft);
                lx += 82.0f;
            }

            g.setColour(muted);
            g.drawText("Frequency →   processed response relative to input",
                       legendArea.removeFromBottom(18.0f), juce::Justification::centredRight);
        }

        juce::String parameterName;
        double displayedParameterValue = 0.0;
        std::vector<Trace> traces;
        juce::Rectangle<float> lastPlot;
    };

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

    SpectralResponseComponent spectralResponse;

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
