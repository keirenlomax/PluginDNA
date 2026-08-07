#pragma once

#include "JuceHeader.h"
#include <cstdint>

struct SineGenerator {
    double sampleRate = 48000.0;
    double frequency = 1000.0;
    double phase = 0.0;
    float amplitude = 0.5f;

    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct NoiseGenerator {
    float amplitude = 0.5f;
    juce::Random rng;

    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct SweepGenerator {
    double sampleRate = 48000.0;
    double startHz = 20.0;
    double endHz = 20000.0;
    double duration = 5.0;
    float amplitude = 0.5f;

    double currentPhase = 0.0;
    double currentFreq = 20.0;
    int64_t currentSample = 0;

    void reset();
    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct InteractionGenerator {
    double sampleRate = 96000.0;
    double duration = 5.0;
    float amplitude = 0.5f;
    int64_t currentSample = 0;
    std::vector<double> phasesBroadband = std::vector<double>(5, 0.0);
    std::vector<double> phasesHighPair = std::vector<double>(2, 0.0);

    void reset();
    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct DeltaGenerator {
    int64_t currentSample = 0;
    int64_t impulseSample = 8192;
    float amplitude = 0.5f;

    void reset() { currentSample = 0; }
    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct BoundaryGenerator {
    double sampleRate = 96000.0;
    double duration = 10.0;
    float amplitude = 0.5f;
    int64_t currentSample = 0;
    std::vector<double> infraPhases = std::vector<double>(6, 0.0);
    std::vector<double> ultraPhases = std::vector<double>(4, 0.0);
    void reset();
    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct StereoTestGenerator {
    double sampleRate = 96000.0;
    double duration = 4.0;
    float amplitude = 0.5f;
    int64_t currentSample = 0;
    double phase = 0.0;
    void reset() { currentSample = 0; phase = 0.0; }
    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct SummingTestGenerator {
    double sampleRate = 96000.0;
    double duration = 3.0;
    float amplitude = 0.5f;
    int64_t currentSample = 0;
    double phaseA = 0.0, phaseB = 0.0;
    void reset() { currentSample = 0; phaseA = 0.0; phaseB = 0.0; }
    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};

struct AliasStressGenerator {
    double sampleRate = 96000.0;
    double duration = 5.0;
    float amplitude = 0.5f;
    int64_t currentSample = 0;
    double phase = 0.0;
    void reset() { currentSample = 0; phase = 0.0; }
    void fillBlock(juce::AudioBuffer<float>& buffer, int numSamples);
};
