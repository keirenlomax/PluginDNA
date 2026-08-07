#pragma once

#include "Analyzer.h"
#include "MeasurementResult.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <vector>

class InteractionAnalyzer final : public Analyzer
{
public:
    InteractionAnalyzer(double sampleRateIn, double durationSecondsIn,
                        std::vector<juce::String> parameterNamesIn)
        : sampleRate(sampleRateIn), durationSeconds(durationSecondsIn),
          parameterNames(std::move(parameterNamesIn))
    {
        result.analyserName = "InteractionDNA";
    }

    void processBlock(const BlockContext& ctx) override
    {
        auto& run = runs[ctx.runId];
        if (!run.initialised)
        {
            run.initialised = true;
            run.inputGainDb = ctx.inputGainDb;
            run.paramValues = ctx.paramNamedValues;
            run.sectionSample = std::max<int64_t>(1, static_cast<int64_t>(durationSeconds * sampleRate / 3.0));
            initialiseBanks(run);
        }

        for (int i = 0; i < ctx.numSamples; ++i)
        {
            const int64_t absoluteSample = ctx.firstSample + i;
            const int test = juce::jlimit(0, 2, static_cast<int>(absoluteSample / run.sectionSample));
            const float input = ctx.inL != nullptr ? ctx.inL[i] : 0.0f;
            const float output = ctx.outL != nullptr ? ctx.outL[i] : 0.0f;
            auto& bank = test == 0 ? run.broadband : (test == 1 ? run.highPair : run.lowImd);
            bank.inputEnergy += static_cast<double>(input) * input;
            bank.sampleCount++;
            for (auto& detector : bank.detectors)
                detector.process(output);
        }
    }

    void finish(const juce::File&) override
    {
        buildResult();
    }

    const MeasurementDataset& getResult() const noexcept { return result; }
    MeasurementDataset takeResult() noexcept { return std::move(result); }

private:
    struct Detector
    {
        double frequency = 0.0;
        int classId = 0;
        double stepCos = 1.0;
        double stepSin = 0.0;
        double phaseCos = 1.0;
        double phaseSin = 0.0;
        int64_t n = 0;

        void configure(double frequencyIn, int classIdIn, double sr)
        {
            frequency = frequencyIn;
            classId = classIdIn;
            const double omega = 2.0 * juce::MathConstants<double>::pi * frequency / sr;
            stepCos = std::cos(omega);
            stepSin = std::sin(omega);
        }

        void process(double sample)
        {
            real += sample * phaseCos;
            imag -= sample * phaseSin;
            const double nextCos = phaseCos * stepCos - phaseSin * stepSin;
            phaseSin = phaseSin * stepCos + phaseCos * stepSin;
            phaseCos = nextCos;
            if ((n & 4095) == 4095)
            {
                const double norm = std::sqrt(phaseCos * phaseCos + phaseSin * phaseSin);
                if (norm > 0.0) { phaseCos /= norm; phaseSin /= norm; }
            }
            ++n;
        }

        double magnitude() const
        {
            return n > 0 ? (2.0 / static_cast<double>(n)) * std::sqrt(real * real + imag * imag) : 0.0;
        }

        double real = 0.0;
        double imag = 0.0;
    };

    struct Bank
    {
        std::vector<Detector> detectors;
        double inputEnergy = 0.0;
        int64_t sampleCount = 0;
    };

    struct RunState
    {
        bool initialised = false;
        int64_t sectionSample = 0;
        float inputGainDb = 0.0f;
        std::map<juce::String, float> paramValues;
        Bank broadband;
        Bank highPair;
        Bank lowImd;
    };

    static std::vector<std::pair<double, int>> buildBroadbandProducts()
    {
        const std::vector<double> tones { 60.0, 250.0, 1000.0, 4000.0, 10000.0 };
        std::vector<std::pair<double, int>> products;
        std::set<int> used;
        const auto add = [&](double frequency, int classId, auto& out, auto& seen)
        {
            if (frequency < 20.0 || frequency > 20000.0) return;
            for (double tone : tones)
                if (std::abs(frequency - tone) < 2.0) return;
            const int key = static_cast<int>(std::lround(frequency));
            if (seen.insert(key).second) out.emplace_back(frequency, classId);
        };
        for (std::size_t i = 0; i < tones.size(); ++i)
            for (std::size_t j = i + 1; j < tones.size(); ++j)
            {
                const double a = tones[i], b = tones[j];
                add(std::abs(b - a), 1, products, used); // difference
                add(a + b, 2, products, used);          // sum
                add(std::abs(2.0 * a - b), 3, products, used);
                add(std::abs(2.0 * b - a), 3, products, used);
            }
        return products;
    }

    static std::vector<std::pair<double, int>> buildHighPairProducts()
    {
        return { {1000.0, 1}, {18000.0, 3}, {21000.0, 3},
                 {38000.0, 4}, {39000.0, 2}, {40000.0, 4} };
    }

    static std::vector<std::pair<double, int>> buildLowImdProducts()
    {
        std::vector<std::pair<double,int>> p;
        for (int n=1;n<=5;++n) { p.emplace_back(7000.0-60.0*n, 5); p.emplace_back(7000.0+60.0*n, 5); }
        p.emplace_back(120.0,3); p.emplace_back(14000.0,3);
        return p;
    }

    void initialiseBank(Bank& bank, const std::vector<std::pair<double, int>>& products)
    {
        bank.detectors.reserve(products.size());
        for (const auto& [frequency, classId] : products)
        {
            if (frequency >= sampleRate * 0.49) continue;
            Detector detector;
            detector.configure(frequency, classId, sampleRate);
            bank.detectors.push_back(detector);
        }
    }

    void initialiseBanks(RunState& run)
    {
        initialiseBank(run.broadband, buildBroadbandProducts());
        initialiseBank(run.highPair, buildHighPairProducts());
        initialiseBank(run.lowImd, buildLowImdProducts());
    }

    static double amplitudeToDb(double amplitude)
    {
        return 20.0 * std::log10(std::max(amplitude, 1.0e-12));
    }

    void appendBankRows(int runId, int testId, const Bank& bank, const RunState& run)
    {
        const double inputRms = bank.sampleCount > 0
            ? std::sqrt(bank.inputEnergy / static_cast<double>(bank.sampleCount)) : 0.0;
        for (const auto& detector : bank.detectors)
        {
            std::vector<double> row;
            row.push_back(runId);
            row.push_back(testId);
            row.push_back(detector.frequency);
            row.push_back(detector.classId);
            row.push_back(amplitudeToDb(detector.magnitude() / std::max(inputRms, 1.0e-12)));
            row.push_back(run.inputGainDb);
            for (const auto& name : parameterNames)
            {
                const auto it = run.paramValues.find(name);
                row.push_back(it != run.paramValues.end() ? it->second : 0.0);
            }
            result.rows.push_back(std::move(row));
        }
    }

    void buildResult()
    {
        result.columns = { "runId", "testId", "productHz", "productClass", "magDbRelativeToInput", "inputGainDb" };
        for (const auto& name : parameterNames)
            result.columns.push_back(name.toStdString());
        for (const auto& [runId, run] : runs)
        {
            appendBankRows(runId, 1, run.broadband, run);
            appendBankRows(runId, 2, run.highPair, run);
            appendBankRows(runId, 3, run.lowImd, run);
        }
    }

    double sampleRate = 96000.0;
    double durationSeconds = 5.0;
    std::vector<juce::String> parameterNames;
    std::map<int, RunState> runs;
    MeasurementDataset result;
};

inline std::unique_ptr<Analyzer> createInteractionAnalyzer(double sampleRate, double durationSeconds,
                                                           const std::vector<juce::String>& parameterNames)
{
    return std::make_unique<InteractionAnalyzer>(sampleRate, durationSeconds, parameterNames);
}
