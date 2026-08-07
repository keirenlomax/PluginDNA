#include "SignalGenerator.h"
#include <cmath>
#include <algorithm>

void SineGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples) {
    const double phaseIncrement = 2.0 * juce::MathConstants<double>::pi * frequency / sampleRate;

    for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
        auto* channelData = buffer.getWritePointer(ch);
        double currentPhase = phase;

        for (int i = 0; i < numSamples; ++i) {
            channelData[i] = amplitude * (float)std::sin(currentPhase);
            currentPhase += phaseIncrement;

            if (currentPhase > 2.0 * juce::MathConstants<double>::pi)
                currentPhase -= 2.0 * juce::MathConstants<double>::pi;
        }
    }

    phase += phaseIncrement * numSamples;
    if (phase > 2.0 * juce::MathConstants<double>::pi)
        phase -= 2.0 * juce::MathConstants<double>::pi;
}

void NoiseGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples) {
    for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
        auto* channelData = buffer.getWritePointer(ch);
        for (int i = 0; i < numSamples; ++i) {
            // Generate white noise in range [-amplitude, amplitude]
            channelData[i] = amplitude * (2.0f * rng.nextFloat() - 1.0f);
        }
    }
}

void SweepGenerator::reset() {
    currentPhase = 0.0;
    currentFreq = startHz;
    currentSample = 0;
}

void SweepGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples) {
    const int64_t totalSamples = (int64_t)(duration * sampleRate);
    const double logStart = std::log(startHz);
    const double logEnd = std::log(endHz);

    // Generate one sweep sample and copy it to every channel. The old
    // channel-first loop advanced currentSample separately for L and R,
    // causing the left channel to skip half the sweep positions.
    for (int i = 0; i < numSamples; ++i) {
        float value = 0.0f;

        if (currentSample < totalSamples) {
            const double t = (double)currentSample / (double)totalSamples;
            const double logFreq = logStart + t * (logEnd - logStart);
            currentFreq = std::exp(logFreq);

            const double phaseIncrement = 2.0 * juce::MathConstants<double>::pi * currentFreq / sampleRate;
            value = amplitude * (float)std::sin(currentPhase);

            currentPhase += phaseIncrement;
            if (currentPhase > 2.0 * juce::MathConstants<double>::pi)
                currentPhase -= 2.0 * juce::MathConstants<double>::pi;

            currentSample++;
        }

        for (int ch = 0; ch < buffer.getNumChannels(); ++ch)
            buffer.getWritePointer(ch)[i] = value;
    }
}


void InteractionGenerator::reset()
{
    currentSample = 0;
    std::fill(phasesBroadband.begin(), phasesBroadband.end(), 0.0);
    std::fill(phasesHighPair.begin(), phasesHighPair.end(), 0.0);
}

void InteractionGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples)
{
    static const std::vector<double> broadband { 60.0, 250.0, 1000.0, 4000.0, 10000.0 };
    static const std::vector<double> highPair { 19000.0, 20000.0 };
    const int64_t halfSample = static_cast<int64_t>(duration * sampleRate * 0.5);

    for (int i = 0; i < numSamples; ++i)
    {
        const bool firstHalf = currentSample < halfSample;
        const auto& frequencies = firstHalf ? broadband : highPair;
        auto& phases = firstHalf ? phasesBroadband : phasesHighPair;
        const double perTone = amplitude / static_cast<double>(frequencies.size());
        double value = 0.0;
        for (std::size_t tone = 0; tone < frequencies.size(); ++tone)
        {
            value += perTone * std::sin(phases[tone]);
            phases[tone] += 2.0 * juce::MathConstants<double>::pi * frequencies[tone] / sampleRate;
            if (phases[tone] > 2.0 * juce::MathConstants<double>::pi)
                phases[tone] -= 2.0 * juce::MathConstants<double>::pi;
        }
        for (int ch = 0; ch < buffer.getNumChannels(); ++ch)
            buffer.getWritePointer(ch)[i] = static_cast<float>(value);
        ++currentSample;
    }
}

void DeltaGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples)
{
    for (int i = 0; i < numSamples; ++i)
    {
        const float value = currentSample == impulseSample ? amplitude : 0.0f;
        for (int ch = 0; ch < buffer.getNumChannels(); ++ch)
            buffer.getWritePointer(ch)[i] = value;
        ++currentSample;
    }
}

void BoundaryGenerator::reset()
{
    currentSample = 0;
    std::fill(infraPhases.begin(), infraPhases.end(), 0.0);
    std::fill(ultraPhases.begin(), ultraPhases.end(), 0.0);
}

void BoundaryGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples)
{
    static const double infra[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0};
    static const double ultra[] = {22000.0, 26000.0, 32000.0, 40000.0};
    const int64_t split = static_cast<int64_t>(duration * sampleRate * 0.5);
    for (int i = 0; i < numSamples; ++i, ++currentSample)
    {
        double sample = 0.0;
        if (currentSample < split)
        {
            for (int t = 0; t < 6; ++t)
            {
                sample += std::sin(infraPhases[(std::size_t)t]);
                infraPhases[(std::size_t)t] += 2.0 * juce::MathConstants<double>::pi * infra[t] / sampleRate;
            }
            sample *= amplitude / 6.0;
        }
        else
        {
            for (int t = 0; t < 4; ++t)
            {
                sample += std::sin(ultraPhases[(std::size_t)t]);
                ultraPhases[(std::size_t)t] += 2.0 * juce::MathConstants<double>::pi * ultra[t] / sampleRate;
            }
            sample *= amplitude / 4.0;
        }
        for (int ch = 0; ch < buffer.getNumChannels(); ++ch)
            buffer.setSample(ch, i, static_cast<float>(sample));
    }
}

void StereoTestGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples)
{
    const int64_t section = std::max<int64_t>(1, static_cast<int64_t>(duration * sampleRate / 4.0));
    for (int i=0;i<numSamples;++i,++currentSample)
    {
        const int mode = juce::jlimit(0,3,(int)(currentSample/section));
        const int64_t local = currentSample % section;
        const float v = amplitude * static_cast<float>(std::sin(2.0 * juce::MathConstants<double>::pi * 997.0 * local / sampleRate));
        float l=0.0f,r=0.0f;
        if(mode==0) l=v; else if(mode==1) r=v; else if(mode==2) l=r=v; else {l=v;r=-v;}
        if(buffer.getNumChannels()>0) buffer.setSample(0,i,l);
        if(buffer.getNumChannels()>1) buffer.setSample(1,i,r);
    }
}

void SummingTestGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples)
{
    const int64_t section = std::max<int64_t>(1, static_cast<int64_t>(duration * sampleRate / 3.0));
    for(int i=0;i<numSamples;++i,++currentSample)
    {
        const int mode=juce::jlimit(0,2,(int)(currentSample/section));
        const int64_t local=currentSample%section;
        const double pa=2.0*juce::MathConstants<double>::pi*997.0*local/sampleRate;
        const double pb=2.0*juce::MathConstants<double>::pi*1531.0*local/sampleRate;
        const float a=0.5f*amplitude*(float)(std::sin(pa)+0.35*std::sin(2.0*pa));
        const float b=0.5f*amplitude*(float)(std::sin(pb)+0.25*std::sin(3.0*pb));
        const float v=mode==0?a:(mode==1?b:a+b);
        for(int ch=0;ch<buffer.getNumChannels();++ch)buffer.setSample(ch,i,v);
    }
}

void AliasStressGenerator::fillBlock(juce::AudioBuffer<float>& buffer, int numSamples)
{
    static const std::vector<double> frequencies {
        500.0, 1000.0, 2000.0, 5000.0, 8000.0,
        10000.0, 12000.0, 15000.0, 18000.0, 20000.0
    };
    const int64_t segment = std::max<int64_t>(1, static_cast<int64_t>(duration * sampleRate / frequencies.size()));
    for (int i = 0; i < numSamples; ++i, ++currentSample)
    {
        const int index = juce::jlimit(0, static_cast<int>(frequencies.size()) - 1,
                                      static_cast<int>(currentSample / segment));
        const double frequency = std::min(frequencies[static_cast<std::size_t>(index)], sampleRate * 0.49);
        const float value = amplitude * static_cast<float>(std::sin(phase));
        phase += 2.0 * juce::MathConstants<double>::pi * frequency / sampleRate;
        if (phase >= 2.0 * juce::MathConstants<double>::pi)
            phase = std::fmod(phase, 2.0 * juce::MathConstants<double>::pi);
        for (int ch = 0; ch < buffer.getNumChannels(); ++ch)
            buffer.setSample(ch, i, value);
    }
}
