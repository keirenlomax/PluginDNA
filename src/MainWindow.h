#pragma once

#include "JuceHeader.h"
#include "MainComponent.h"

class MainWindow : public juce::DocumentWindow {
public:
    MainWindow(juce::String name);
    ~MainWindow() override;

    void closeButtonPressed() override;
    void resized() override;

private:
    juce::Viewport* viewport = nullptr;
    MainComponent* mainComponent = nullptr;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainWindow)
};
