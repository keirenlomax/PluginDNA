#include "JuceHeader.h"
#include "MainWindow.h"

class PluginDNAApplication : public juce::JUCEApplication {
public:
    PluginDNAApplication() {}

    const juce::String getApplicationName() override {
        return "PluginDNA";
    }
    const juce::String getApplicationVersion() override {
        return "12.1.0";
    }
    bool moreThanOneInstanceAllowed() override {
        return true;
    }

    void initialise(const juce::String& commandLine) override {
        mainWindow = std::make_unique<MainWindow>(getApplicationName());
        mainWindow->setVisible(true);
    }

    void shutdown() override {
        mainWindow = nullptr;
    }

    void systemRequestedQuit() override {
        quit();
    }

private:
    std::unique_ptr<MainWindow> mainWindow;
};

START_JUCE_APPLICATION(PluginDNAApplication)
