#include "MainWindow.h"

MainWindow::MainWindow(juce::String name)
    : DocumentWindow(
          name,
          juce::Desktop::getInstance().getDefaultLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId),
          DocumentWindow::allButtons) {
    setUsingNativeTitleBar(true);

    auto* newViewport = new juce::Viewport();
    auto* newMainComponent = new MainComponent();

    viewport = newViewport;
    mainComponent = newMainComponent;

    viewport->setScrollBarsShown(true, false);
    viewport->setScrollBarThickness(14);
    viewport->setViewedComponent(mainComponent, true);

    setContentOwned(viewport, true);

#if JUCE_IOS || JUCE_ANDROID
    setFullScreen(true);
#else
    setResizable(true, true);
    setResizeLimits(900, 600, 2200, 1800);
    centreWithSize(1200, 800);
#endif

    resized();
    setVisible(true);
    toFront(true);
}

MainWindow::~MainWindow() {}

void MainWindow::resized() {
    DocumentWindow::resized();

    if (viewport == nullptr || mainComponent == nullptr)
        return;

    const int availableWidth = juce::jmax(900, viewport->getMaximumVisibleWidth());
    const int contentHeight = 2420; // V22.1: include DNA Map below summary

    mainComponent->setSize(availableWidth, contentHeight);
}

void MainWindow::closeButtonPressed() {
    juce::JUCEApplication::getInstance()->systemRequestedQuit();
}
