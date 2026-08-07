#include "MainComponent.h"
#include "RmsPeakAnalyzer.h"
#include "ThdAnalyzer.h"
#include "LinearResponseAnalyzer.h"
#include "InteractionAnalyzer.h"
#include "StereoAnalyzer.h"
#include "SummingAnalyzer.h"
#include "AliasAnalyzer.h"
#include "TimingAnalyzer.h"
#include "ResidualAnalyzer.h"
#include "BoundaryAnalyzer.h"
#include "BucketSpec.h"
#include "TransferCurveAnalyzer.h"
#include "RawCsvAnalyzer.h"
#include "MeasurementConfigComponent.h"
#include "ParameterConfigComponent.h"
#include <thread>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

#include <fstream>
#include <iomanip>

namespace
{
    int findColumn(const MeasurementDataset& dataset, const std::string& name)
    {
        for (std::size_t i = 0; i < dataset.columns.size(); ++i)
            if (dataset.columns[i] == name)
                return static_cast<int>(i);

        return -1;
    }

    double amplitudeToDb(double amplitude)
    {
        return 20.0 * std::log10(std::max(amplitude, 1.0e-12));
    }

    struct LevelSummary
    {
        juce::String level;
        juce::String peaks;
        juce::String dynamics;
    };

    LevelSummary makeLevelSummary(const MeasurementDataset& dataset)
    {
        const int rmsIn = findColumn(dataset, "rmsInL");
        const int rmsOut = findColumn(dataset, "rmsOutL");
        const int peakIn = findColumn(dataset, "peakInL");
        const int peakOut = findColumn(dataset, "peakOutL");

        if (rmsIn < 0 || rmsOut < 0 || peakIn < 0 || peakOut < 0 || dataset.rows.empty())
            return { "LEVEL MATCH  unavailable", "TRANSIENT IMPACT  unavailable", "DYNAMIC RANGE  unavailable" };

        double levelChange = 0.0;
        double crestChange = 0.0;
        int validRows = 0;

        for (const auto& row : dataset.rows)
        {
            const int lastIndex = std::max(std::max(rmsIn, rmsOut), std::max(peakIn, peakOut));
            if (static_cast<int>(row.size()) <= lastIndex)
                continue;

            const double inRms = row[static_cast<std::size_t>(rmsIn)];
            const double outRms = row[static_cast<std::size_t>(rmsOut)];
            const double inPeak = row[static_cast<std::size_t>(peakIn)];
            const double outPeak = row[static_cast<std::size_t>(peakOut)];

            if (inRms <= 0.0 || outRms <= 0.0 || inPeak <= 0.0 || outPeak <= 0.0)
                continue;

            const double rowLevelChange = amplitudeToDb(outRms) - amplitudeToDb(inRms);
            const double rowPeakChange = amplitudeToDb(outPeak) - amplitudeToDb(inPeak);

            levelChange += rowLevelChange;
            crestChange += rowPeakChange - rowLevelChange;
            ++validRows;
        }

        if (validRows == 0)
            return { "LEVEL MATCH  unavailable", "TRANSIENT IMPACT  unavailable", "DYNAMIC RANGE  unavailable" };

        levelChange /= validRows;
        crestChange /= validRows;

        juce::String levelText;
        if (std::abs(levelChange) < 0.25)
        {
            levelText = "LEVEL MATCH  No adjustment needed";
        }
        else
        {
            const double suggestedTrim = -levelChange;
            const juce::String sign = suggestedTrim > 0.0 ? "+" : "";
            levelText = "LEVEL MATCH  Suggested output trim  " + sign
                      + juce::String(suggestedTrim, 2) + " dB";
        }

        juce::String transientWords;
        if (crestChange > 0.5)
            transientWords = "Punchier";
        else if (crestChange < -0.5)
            transientWords = "Softer transients";
        else
            transientWords = "Little change";

        juce::String rangeWords;
        if (crestChange > 0.5)
            rangeWords = "More open";
        else if (crestChange < -0.5)
            rangeWords = "More controlled";
        else
            rangeWords = "Largely unchanged";

        const juce::String changeSign = crestChange > 0.0 ? "+" : "";
        const juce::String measuredChange = changeSign + juce::String(crestChange, 2) + " dB";

        return {
            levelText,
            "TRANSIENT IMPACT  " + transientWords + "  (" + measuredChange + ")",
            "DYNAMIC RANGE  " + rangeWords + "  (" + measuredChange + ")"
        };
    }

    struct HarmonicSummary
    {
        juce::String identity;
        juce::String balance;
        juce::String growth;
    };

    HarmonicSummary makeHarmonicSummary(const MeasurementDataset& dataset)
    {
        const int thdColumn = findColumn(dataset, "thd");
        const int gainColumn = findColumn(dataset, "inputGainDb");
        std::array<int, 9> harmonicColumns {};
        for (int harmonic = 2; harmonic <= 10; ++harmonic)
            harmonicColumns[static_cast<std::size_t>(harmonic - 2)] =
                findColumn(dataset, "h" + std::to_string(harmonic));

        if (thdColumn < 0 || gainColumn < 0 || harmonicColumns[0] < 0 || dataset.rows.empty())
            return { "HARMONIC DNA  unavailable", "BALANCE  unavailable", "WHEN PUSHED  unavailable" };

        struct GainProfile
        {
            std::array<double, 9> sums {};
            int count = 0;
        };

        std::map<int, GainProfile> profiles;
        const int lastColumn = std::max(gainColumn, harmonicColumns.back());

        for (const auto& row : dataset.rows)
        {
            if (static_cast<int>(row.size()) <= lastColumn)
                continue;

            const int gain = static_cast<int>(std::lround(row[static_cast<std::size_t>(gainColumn)]));
            std::array<double, 9> rowHarmonics {};
            bool valid = true;
            for (std::size_t i = 0; i < harmonicColumns.size(); ++i)
            {
                const double value = row[static_cast<std::size_t>(harmonicColumns[i])];
                if (!std::isfinite(value) || value < 0.0)
                {
                    valid = false;
                    break;
                }
                rowHarmonics[i] = value;
            }

            if (valid)
            {
                auto& profile = profiles[gain];
                for (std::size_t i = 0; i < rowHarmonics.size(); ++i)
                    profile.sums[i] += rowHarmonics[i];
                ++profile.count;
            }
        }

        auto averageProfile = [&profiles](int gain) -> std::array<double, 9>
        {
            std::array<double, 9> result {};
            const auto it = profiles.find(gain);
            if (it == profiles.end() || it->second.count == 0)
                return result;
            for (std::size_t i = 0; i < result.size(); ++i)
                result[i] = it->second.sums[i] / it->second.count;
            return result;
        };

        // -18 dBFS is the normal gain-staging reference. If unavailable, use
        // the least-driven measured profile as the honest fallback.
        int referenceGain = -18;
        if (profiles.find(referenceGain) == profiles.end() && !profiles.empty())
            referenceGain = profiles.begin()->first;

        int pushedGain = -3;
        if (profiles.find(pushedGain) == profiles.end() && !profiles.empty())
            pushedGain = profiles.rbegin()->first;

        const auto reference = averageProfile(referenceGain);
        const auto pushed = averageProfile(pushedGain);

        auto sumRange = [](const std::array<double, 9>& values, int firstHarmonic, int lastHarmonic)
        {
            double sum = 0.0;
            for (int harmonic = firstHarmonic; harmonic <= lastHarmonic; ++harmonic)
                sum += values[static_cast<std::size_t>(harmonic - 2)];
            return sum;
        };

        const double refTotal = sumRange(reference, 2, 10);
        double refPower = 0.0;
        for (double harmonic : reference)
            refPower += harmonic * harmonic;
        const double referenceThd = std::sqrt(refPower);

        if (referenceThd < 0.0005 || refTotal <= 1.0e-10)
            return { "HARMONIC DNA  Very clean", "BALANCE  No clear harmonic pattern", "WHEN PUSHED  Little harmonic growth" };

        int strongestHarmonic = 2;
        for (int harmonic = 3; harmonic <= 10; ++harmonic)
            if (reference[static_cast<std::size_t>(harmonic - 2)] >
                reference[static_cast<std::size_t>(strongestHarmonic - 2)])
                strongestHarmonic = harmonic;

        const double strongestShare =
            reference[static_cast<std::size_t>(strongestHarmonic - 2)] / refTotal;
        const double lowOrder = sumRange(reference, 2, 4);
        const double highOrder = sumRange(reference, 5, 10);

        juce::String identity;
        if (strongestShare >= 0.55)
        {
            const juce::String suffix = strongestHarmonic == 2 ? "nd"
                                      : strongestHarmonic == 3 ? "rd"
                                      : "th";
            identity = "Mostly " + juce::String(strongestHarmonic) + suffix + " harmonic";
        }
        else if (lowOrder > highOrder * 2.5)
            identity = "Smooth low-order harmonics";
        else if (highOrder > lowOrder)
            identity = "Dense high-order harmonics";
        else
            identity = "Mixed harmonic stack";

        double even = 0.0;
        double odd = 0.0;
        for (int harmonic = 2; harmonic <= 10; ++harmonic)
        {
            const double value = reference[static_cast<std::size_t>(harmonic - 2)];
            if ((harmonic % 2) == 0)
                even += value;
            else
                odd += value;
        }

        juce::String balance;
        if (even > odd * 1.6)
            balance = "Even-heavy";
        else if (odd > even * 1.6)
            balance = "Odd-heavy";
        else
            balance = "Even / odd mix";

        const double pushedTotal = sumRange(pushed, 2, 10);
        const double refHighShare = highOrder / std::max(refTotal, 1.0e-12);
        const double pushedHigh = sumRange(pushed, 5, 10);
        const double pushedHighShare = pushedHigh / std::max(pushedTotal, 1.0e-12);

        juce::String growth;
        if (pushedTotal < refTotal * 1.5)
            growth = "Pattern stays restrained";
        else if (pushedHighShare > refHighShare + 0.15)
            growth = "Higher harmonics join in";
        else if (pushed[1] > reference[1] * 2.0 || pushed[3] > reference[3] * 2.0)
            growth = "Odd harmonics rise strongly";
        else
            growth = "Same character, more intense";

        return {
            "HARMONIC DNA  " + identity + "  (at " + juce::String(referenceGain) + " dBFS)",
            "BALANCE  " + balance,
            "WHEN PUSHED  " + growth + "  (to " + juce::String(pushedGain) + " dBFS)"
        };
    }


    struct WaveformSummary
    {
        juce::String shape;
        juce::String detail;
        juce::String style;
    };

    WaveformSummary makeWaveformSummary(const MeasurementDataset& dataset)
    {
        const int runIdColumn = findColumn(dataset, "runId");
        const int inputColumn = findColumn(dataset, "inL");
        const int outputColumn = findColumn(dataset, "outL");

        if (runIdColumn < 0 || inputColumn < 0 || outputColumn < 0 || dataset.rows.empty())
            return { "WAVE SHAPE  unavailable", "WHAT IT DOES  unavailable", "SIMILAR TO  unavailable" };

        struct RunSamples
        {
            std::vector<double> input;
            std::vector<double> output;
        };

        std::map<int, RunSamples> runs;
        const int lastColumn = std::max(runIdColumn, std::max(inputColumn, outputColumn));

        // Keep contiguous samples from each run. The old global stride could lock
        // onto a small set of phases for periodic sine tests (for example every
        // 16th sample of a 48-sample cycle), leaving too few input amplitudes to
        // build a transfer curve. Contiguous samples cover the full waveform.
        constexpr std::size_t samplesPerRun = 24000;
        constexpr std::size_t settleSamples = 2048;
        std::map<int, std::size_t> seenPerRun;

        for (const auto& row : dataset.rows)
        {
            if (static_cast<int>(row.size()) <= lastColumn)
                continue;

            const int runId = static_cast<int>(std::llround(row[static_cast<std::size_t>(runIdColumn)]));
            auto& seen = seenPerRun[runId];
            ++seen;

            // Ignore the first short section so plugin state and filters can settle.
            if (seen <= settleSamples)
                continue;

            auto& run = runs[runId];
            if (run.input.size() >= samplesPerRun)
                continue;

            const double x = row[static_cast<std::size_t>(inputColumn)];
            const double y = row[static_cast<std::size_t>(outputColumn)];
            if (!std::isfinite(x) || !std::isfinite(y))
                continue;

            run.input.push_back(x);
            run.output.push_back(y);
        }

        struct ShapeMetrics
        {
            double score = -1.0;
            double asymmetry = 0.0;
            double edgeRatio = 1.0;
            double slopeVariation = 0.0;
            double scatter = 0.0;
            bool foldback = false;
            bool abrupt = false;
            bool timeDependent = false;
        } best;

        for (const auto& [runId, samples] : runs)
        {
            juce::ignoreUnused(runId);
            const int sampleCount = static_cast<int>(std::min(samples.input.size(), samples.output.size()));
            if (sampleCount < 512)
                continue;

            // Align output to input before making x -> y pairs. A small plugin delay or
            // phase shift otherwise looks like a loop and can be mistaken for foldback.
            int bestLag = 0;
            double bestCorrelation = -1.0;
            for (int lag = -64; lag <= 64; ++lag)
            {
                double sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
                int count = 0;
                const int begin = std::max(0, -lag);
                const int finish = std::min(sampleCount, sampleCount - lag);
                for (int i = begin; i < finish; i += 4)
                {
                    const double x = samples.input[static_cast<std::size_t>(i)];
                    const double y = samples.output[static_cast<std::size_t>(i + lag)];
                    sumXX += x * x;
                    sumYY += y * y;
                    sumXY += x * y;
                    ++count;
                }

                if (count > 0 && sumXX > 1.0e-12 && sumYY > 1.0e-12)
                {
                    const double correlation = std::abs(sumXY / std::sqrt(sumXX * sumYY));
                    if (correlation > bestCorrelation)
                    {
                        bestCorrelation = correlation;
                        bestLag = lag;
                    }
                }
            }

            double maxInput = 0.0;
            const int alignedBegin = std::max(0, -bestLag);
            const int alignedEnd = std::min(sampleCount, sampleCount - bestLag);
            for (int i = alignedBegin; i < alignedEnd; ++i)
                maxInput = std::max(maxInput, std::abs(samples.input[static_cast<std::size_t>(i)]));
            if (maxInput <= 1.0e-9)
                continue;

            constexpr int binCount = 129;
            std::vector<double> sums(binCount, 0.0);
            std::vector<double> sumSquares(binCount, 0.0);
            std::vector<int> counts(binCount, 0);

            for (int i = alignedBegin; i < alignedEnd; ++i)
            {
                const double x = samples.input[static_cast<std::size_t>(i)] / maxInput;
                const double y = samples.output[static_cast<std::size_t>(i + bestLag)];
                const int bin = juce::jlimit(0, binCount - 1,
                    static_cast<int>(std::llround((x + 1.0) * 0.5 * (binCount - 1))));
                sums[static_cast<std::size_t>(bin)] += y;
                sumSquares[static_cast<std::size_t>(bin)] += y * y;
                ++counts[static_cast<std::size_t>(bin)];
            }

            std::vector<double> curve(binCount, std::numeric_limits<double>::quiet_NaN());
            double scatterSum = 0.0;
            int scatterCount = 0;
            for (int bin = 0; bin < binCount; ++bin)
            {
                const int count = counts[static_cast<std::size_t>(bin)];
                if (count < 6)
                    continue;

                const double mean = sums[static_cast<std::size_t>(bin)] / count;
                curve[static_cast<std::size_t>(bin)] = mean;
                const double variance = std::max(0.0,
                    sumSquares[static_cast<std::size_t>(bin)] / count - mean * mean);
                scatterSum += std::sqrt(variance);
                ++scatterCount;
            }

            // A sine wave only visits a finite set of amplitudes each cycle.
            // At 96 kHz with a 1 kHz test tone that is roughly 49 distinct
            // magnitudes across each half-cycle, so many of the 129 bins are
            // legitimately empty. Interpolate between measured bins before
            // calculating slopes; otherwise valid sine runs appear unavailable.
            int firstValidBin = -1;
            int lastValidBin = -1;
            for (int bin = 0; bin < binCount; ++bin)
            {
                if (std::isfinite(curve[static_cast<std::size_t>(bin)]))
                {
                    if (firstValidBin < 0)
                        firstValidBin = bin;
                    lastValidBin = bin;
                }
            }

            if (firstValidBin < 0 || lastValidBin <= firstValidBin)
                continue;

            // Extend the measured end values to the outer bins. A properly
            // generated sine normally reaches both ends, but this also handles
            // phase/sample-rate combinations that miss the exact peak sample.
            for (int bin = 0; bin < firstValidBin; ++bin)
                curve[static_cast<std::size_t>(bin)] = curve[static_cast<std::size_t>(firstValidBin)];
            for (int bin = lastValidBin + 1; bin < binCount; ++bin)
                curve[static_cast<std::size_t>(bin)] = curve[static_cast<std::size_t>(lastValidBin)];

            int left = firstValidBin;
            while (left < lastValidBin)
            {
                int right = left + 1;
                while (right <= lastValidBin &&
                       !std::isfinite(curve[static_cast<std::size_t>(right)]))
                    ++right;

                if (right > lastValidBin)
                    break;

                const double y0 = curve[static_cast<std::size_t>(left)];
                const double y1 = curve[static_cast<std::size_t>(right)];
                const int gap = right - left;
                for (int bin = left + 1; bin < right; ++bin)
                {
                    const double t = static_cast<double>(bin - left) / gap;
                    curve[static_cast<std::size_t>(bin)] = y0 + t * (y1 - y0);
                }
                left = right;
            }

            // Smooth the completed transfer curve. Tiny local reversals from
            // dither, modulation or numerical noise must not count as foldback.
            std::vector<double> smooth = curve;
            for (int bin = 2; bin < binCount - 2; ++bin)
            {
                double sum = 0.0;
                int count = 0;
                for (int offset = -2; offset <= 2; ++offset)
                {
                    const double value = curve[static_cast<std::size_t>(bin + offset)];
                    if (std::isfinite(value))
                    {
                        sum += value;
                        ++count;
                    }
                }
                if (count >= 3)
                    smooth[static_cast<std::size_t>(bin)] = sum / count;
            }

            double minY = std::numeric_limits<double>::infinity();
            double maxY = -std::numeric_limits<double>::infinity();
            for (double value : smooth)
                if (std::isfinite(value))
                {
                    minY = std::min(minY, value);
                    maxY = std::max(maxY, value);
                }

            const double outputRange = maxY - minY;
            if (!std::isfinite(outputRange) || outputRange <= 1.0e-9)
                continue;

            std::vector<double> slopes(binCount - 1, std::numeric_limits<double>::quiet_NaN());
            for (int bin = 0; bin < binCount - 1; ++bin)
                if (std::isfinite(smooth[static_cast<std::size_t>(bin)]) &&
                    std::isfinite(smooth[static_cast<std::size_t>(bin + 1)]))
                    slopes[static_cast<std::size_t>(bin)] =
                        (smooth[static_cast<std::size_t>(bin + 1)] - smooth[static_cast<std::size_t>(bin)]) * 64.0;

            auto meanSlope = [&](int first, int last)
            {
                double sum = 0.0;
                int count = 0;
                for (int i = first; i <= last; ++i)
                {
                    const double slope = slopes[static_cast<std::size_t>(i)];
                    if (std::isfinite(slope))
                    {
                        sum += slope;
                        ++count;
                    }
                }
                return count > 0 ? sum / count : 0.0;
            };

            const double centreSlope = meanSlope(48, 79);
            const double edgeSlope = 0.5 * (meanSlope(5, 23) + meanSlope(104, 122));
            if (std::abs(centreSlope) <= 1.0e-9)
                continue;

            const double centreY = std::isfinite(smooth[64]) ? smooth[64] : 0.0;
            double asymmetrySum = 0.0;
            int asymmetryCount = 0;
            for (int offset = 8; offset <= 56; offset += 4)
            {
                const double negative = smooth[static_cast<std::size_t>(64 - offset)];
                const double positive = smooth[static_cast<std::size_t>(64 + offset)];
                if (std::isfinite(negative) && std::isfinite(positive))
                {
                    asymmetrySum += std::abs((positive + negative) - 2.0 * centreY) / outputRange;
                    ++asymmetryCount;
                }
            }

            double slopeMean = 0.0;
            int slopeCount = 0;
            for (double slope : slopes)
                if (std::isfinite(slope))
                {
                    slopeMean += slope;
                    ++slopeCount;
                }
            slopeMean = slopeCount > 0 ? slopeMean / slopeCount : 0.0;

            double slopeVariance = 0.0;
            for (double slope : slopes)
                if (std::isfinite(slope))
                    slopeVariance += (slope - slopeMean) * (slope - slopeMean);

            int transitionBins = 0;
            for (int bin = 65; bin < 124; ++bin)
            {
                const double slope = slopes[static_cast<std::size_t>(bin)];
                if (!std::isfinite(slope))
                    continue;
                const double ratio = std::abs(slope / centreSlope);
                if (ratio < 0.85 && ratio > 0.20)
                    ++transitionBins;
            }

            // Require a sustained reversal with a meaningful drop. One or two noisy
            // negative slopes are not foldback.
            int longestNegativeRun = 0;
            int currentNegativeRun = 0;
            double largestReversal = 0.0;
            double reversalStart = 0.0;
            bool inReversal = false;
            for (int bin = 8; bin < binCount - 9; ++bin)
            {
                const double slope = slopes[static_cast<std::size_t>(bin)];
                const bool clearlyNegative = std::isfinite(slope) &&
                    slope < -0.18 * std::abs(centreSlope);

                if (clearlyNegative)
                {
                    if (!inReversal)
                    {
                        reversalStart = smooth[static_cast<std::size_t>(bin)];
                        inReversal = true;
                    }
                    ++currentNegativeRun;
                    longestNegativeRun = std::max(longestNegativeRun, currentNegativeRun);
                    if (std::isfinite(reversalStart) && std::isfinite(smooth[static_cast<std::size_t>(bin + 1)]))
                        largestReversal = std::max(largestReversal,
                            reversalStart - smooth[static_cast<std::size_t>(bin + 1)]);
                }
                else
                {
                    currentNegativeRun = 0;
                    inReversal = false;
                }
            }

            ShapeMetrics metrics;
            metrics.asymmetry = asymmetryCount > 0 ? asymmetrySum / asymmetryCount : 0.0;
            metrics.edgeRatio = std::abs(edgeSlope / centreSlope);
            metrics.slopeVariation = slopeCount > 0
                ? std::sqrt(slopeVariance / slopeCount) / std::abs(centreSlope)
                : 0.0;
            metrics.scatter = scatterCount > 0
                ? (scatterSum / scatterCount) / outputRange
                : 0.0;
            metrics.foldback = longestNegativeRun >= 6 && largestReversal > 0.06 * outputRange;
            metrics.timeDependent = !metrics.foldback && metrics.scatter > 0.035;
            metrics.abrupt = transitionBins <= 4 && metrics.edgeRatio < 0.30;
            metrics.score = std::abs(1.0 - metrics.edgeRatio)
                + metrics.slopeVariation
                + metrics.asymmetry * 2.0
                + (metrics.foldback ? 2.5 : 0.0)
                + (metrics.timeDependent ? 0.8 : 0.0);

            if (metrics.score > best.score)
                best = metrics;
        }

        if (best.score < 0.0)
            return { "WAVE SHAPE  unavailable", "WHAT IT DOES  unavailable", "SIMILAR TO  unavailable" };

        const bool asymmetric = best.asymmetry > 0.035;
        juce::String shape;
        juce::String action;
        juce::String style;

        if (best.foldback)
        {
            shape = asymmetric ? "Asymmetric foldback shaping" : "Balanced foldback shaping";
            action = "Wave clearly folds back after reaching its peak";
            style = "Creative foldback distortion / wavefolder behaviour";
        }
        else if (best.timeDependent)
        {
            shape = asymmetric ? "Moving asymmetric saturation" : "Time-dependent saturation";
            action = "Shape changes through the waveform instead of following one fixed curve";
            style = asymmetric
                ? "Vintage preamp / dynamic analogue-style behaviour"
                : "Tape / dynamic analogue-style behaviour";
        }
        else if (best.edgeRatio < 0.12 && best.abrupt)
        {
            shape = asymmetric ? "Asymmetric hard clipping" : "Hard clipping";
            action = asymmetric ? "Cuts one side of the peaks harder" : "Cuts peaks sharply";
            style = "Modern digital clipper / brickwall limiter behaviour";
        }
        else if (best.edgeRatio < 0.12)
        {
            shape = asymmetric ? "Asymmetric heavy saturation" : "Heavy smooth saturation";
            action = asymmetric ? "Rounds one half harder, then flattens peaks" : "Rounds and flattens peaks strongly";
            style = asymmetric
                ? "Driven vintage preamp / transformer-style behaviour"
                : "Tape-style soft limiting / driven analogue stage";
        }
        else if (best.edgeRatio < 0.65)
        {
            shape = asymmetric ? "Asymmetric smooth saturation" : "Smooth balanced saturation";
            action = asymmetric ? "Rounds the two halves differently" : "Rounds peaks gradually";
            style = asymmetric
                ? "Tube / vintage preamp-style behaviour"
                : "Tube / tape-style saturation behaviour";
        }
        else if (best.edgeRatio > 1.25)
        {
            shape = asymmetric ? "Asymmetric expansion" : "Expanding waveshaping";
            action = "Pushes peaks outward before the edges";
            style = "Creative waveshaper / exciter behaviour";
        }
        else
        {
            shape = asymmetric ? "Slightly uneven, mostly clean" : "Clean and mostly unchanged";
            action = asymmetric ? "Treats each half slightly differently" : "Leaves the waveform shape largely intact";
            style = asymmetric ? "Very subtle analogue-style imbalance" : "Clean gain / utility behaviour";
        }

        return {
            "WAVE SHAPE  " + shape,
            "WHAT IT DOES  " + action,
            "SIMILAR TO  " + style
        };
    }


    struct TemporalSummary
    {
        juce::String response;
        juce::String attack;
        juce::String after;
    };

    TemporalSummary makeTemporalSummary(const MeasurementDataset& dataset)
    {
        const int runColumn = findColumn(dataset, "runId");
        const int timeColumn = findColumn(dataset, "time_sec");
        const int inputColumn = findColumn(dataset, "inL");
        const int outputColumn = findColumn(dataset, "outL");
        if (runColumn < 0 || timeColumn < 0 || inputColumn < 0 || outputColumn < 0 || dataset.rows.empty())
            return { "TEMPORAL RESPONSE  unavailable", "ATTACK HANDLING  unavailable", "AFTER-EFFECTS  unavailable" };

        struct Samples { std::vector<double> t, x, y; };
        std::map<int, Samples> runs;
        constexpr std::size_t maxSamplesPerRun = 96000;
        const int required = std::max({ runColumn, timeColumn, inputColumn, outputColumn });
        for (const auto& row : dataset.rows)
        {
            if (static_cast<int>(row.size()) <= required)
                continue;
            const int runId = static_cast<int>(std::llround(row[static_cast<std::size_t>(runColumn)]));
            auto& r = runs[runId];
            if (r.x.size() >= maxSamplesPerRun)
                continue;
            const double t = row[static_cast<std::size_t>(timeColumn)];
            const double x = row[static_cast<std::size_t>(inputColumn)];
            const double y = row[static_cast<std::size_t>(outputColumn)];
            if (!std::isfinite(t) || !std::isfinite(x) || !std::isfinite(y))
                continue;
            r.t.push_back(t); r.x.push_back(x); r.y.push_back(y);
        }

        double strongestMemory = 0.0;
        double strongestAttackChange = 0.0;
        double strongestOvershoot = 0.0;
        double strongestSettling = 0.0;
        bool found = false;

        for (const auto& [runId, r] : runs)
        {
            juce::ignoreUnused(runId);
            const int n = static_cast<int>(std::min({ r.t.size(), r.x.size(), r.y.size() }));
            if (n < 4096)
                continue;

            std::vector<int> crossings;
            for (int i = 1; i < n; ++i)
                if (r.x[static_cast<std::size_t>(i - 1)] <= 0.0 && r.x[static_cast<std::size_t>(i)] > 0.0)
                    crossings.push_back(i);
            if (crossings.size() < 12)
                continue;

            const int period = crossings[2] - crossings[1];
            if (period < 8)
                continue;

            // Find the small delay that best aligns the settled output with the input.
            const int steadyBegin = crossings[std::min<std::size_t>(8, crossings.size() - 2)];
            const int steadyEnd = std::min(n, steadyBegin + period * 24);
            int bestLag = 0;
            double bestCorr = -1.0;
            for (int lag = -128; lag <= 128; ++lag)
            {
                double xx = 0.0, yy = 0.0, xy = 0.0;
                int count = 0;
                const int begin = std::max(steadyBegin, -lag);
                const int end = std::min(steadyEnd, n - lag);
                for (int i = begin; i < end; i += 2)
                {
                    const double x = r.x[static_cast<std::size_t>(i)];
                    const double y = r.y[static_cast<std::size_t>(i + lag)];
                    xx += x * x; yy += y * y; xy += x * y; ++count;
                }
                if (count > 0 && xx > 1.0e-12 && yy > 1.0e-12)
                {
                    const double corr = std::abs(xy / std::sqrt(xx * yy));
                    if (corr > bestCorr) { bestCorr = corr; bestLag = lag; }
                }
            }

            auto cycleGain = [&](int a, int b) -> double
            {
                a = std::max(a, -bestLag);
                b = std::min({ b, n, n - bestLag });
                double inPeak = 0.0, outPeak = 0.0;
                for (int i = a; i < b; ++i)
                {
                    inPeak = std::max(inPeak, std::abs(r.x[static_cast<std::size_t>(i)]));
                    outPeak = std::max(outPeak, std::abs(r.y[static_cast<std::size_t>(i + bestLag)]));
                }
                return inPeak > 1.0e-9 ? outPeak / inPeak : 0.0;
            };

            std::vector<double> steadyGains;
            for (std::size_t c = 8; c + 1 < crossings.size() && c < 24; ++c)
                steadyGains.push_back(cycleGain(crossings[c], crossings[c + 1]));
            if (steadyGains.empty())
                continue;
            std::sort(steadyGains.begin(), steadyGains.end());
            const std::size_t mid = steadyGains.size() / 2;
            const double steadyGain = (steadyGains.size() % 2 == 0)
                ? 0.5 * (steadyGains[mid - 1] + steadyGains[mid])
                : steadyGains[mid];
            if (!std::isfinite(steadyGain) || steadyGain <= 1.0e-9)
                continue;

            const double firstGain = cycleGain(crossings[0], crossings[1]);
            const double attackChange = firstGain / steadyGain - 1.0;
            strongestAttackChange = std::abs(attackChange) > std::abs(strongestAttackChange)
                ? attackChange : strongestAttackChange;
            strongestOvershoot = std::max(strongestOvershoot, attackChange);

            // Watch how quickly the cycle gain settles. Alternating or persistent
            // deviations suggest ringing/overshoot rather than a single static curve.
            int unsettledCycles = 0;
            double settlingEnergy = 0.0;
            for (std::size_t c = 0; c + 1 < crossings.size() && c < 8; ++c)
            {
                const double g = cycleGain(crossings[c], crossings[c + 1]);
                const double deviation = std::abs(g / steadyGain - 1.0);
                settlingEnergy += deviation;
                if (deviation > 0.035)
                    ++unsettledCycles;
            }
            strongestSettling = std::max(strongestSettling,
                settlingEnergy / 8.0 + 0.02 * unsettledCycles);

            // Quantify history dependence after delay alignment. For a memoryless
            // processor, the same input amplitude maps tightly to one output value.
            constexpr int bins = 65;
            std::vector<double> sum(bins, 0.0), sumSq(bins, 0.0);
            std::vector<int> count(bins, 0);
            double maxX = 0.0, maxY = 0.0;
            for (int i = steadyBegin; i < steadyEnd && i + bestLag < n; ++i)
            {
                maxX = std::max(maxX, std::abs(r.x[static_cast<std::size_t>(i)]));
                maxY = std::max(maxY, std::abs(r.y[static_cast<std::size_t>(i + bestLag)]));
            }
            if (maxX <= 1.0e-9 || maxY <= 1.0e-9)
                continue;
            for (int i = steadyBegin; i < steadyEnd && i + bestLag < n; ++i)
            {
                const double x = r.x[static_cast<std::size_t>(i)] / maxX;
                const double y = r.y[static_cast<std::size_t>(i + bestLag)];
                const int bin = juce::jlimit(0, bins - 1,
                    static_cast<int>(std::llround((x + 1.0) * 0.5 * (bins - 1))));
                sum[static_cast<std::size_t>(bin)] += y;
                sumSq[static_cast<std::size_t>(bin)] += y * y;
                ++count[static_cast<std::size_t>(bin)];
            }
            double scatter = 0.0; int scatterBins = 0;
            for (int b = 0; b < bins; ++b)
            {
                if (count[static_cast<std::size_t>(b)] < 8)
                    continue;
                const double mean = sum[static_cast<std::size_t>(b)] / count[static_cast<std::size_t>(b)];
                const double var = std::max(0.0,
                    sumSq[static_cast<std::size_t>(b)] / count[static_cast<std::size_t>(b)] - mean * mean);
                scatter += std::sqrt(var) / maxY;
                ++scatterBins;
            }
            if (scatterBins > 0)
                strongestMemory = std::max(strongestMemory, scatter / scatterBins);
            found = true;
        }

        if (!found)
            return { "TEMPORAL RESPONSE  unavailable", "ATTACK HANDLING  unavailable", "AFTER-EFFECTS  unavailable" };

        juce::String response;
        if (strongestMemory > 0.055)
            response = "Clearly stateful - depends on recent samples";
        else if (strongestMemory > 0.022)
            response = "Some signal memory / time-dependent behaviour";
        else
            response = "Mostly immediate - little sample memory detected";

        juce::String attack;
        if (strongestAttackChange > 0.10)
            attack = "Preserves or lifts the first attack";
        else if (strongestAttackChange < -0.10)
            attack = "Softens the first attack";
        else
            attack = "Leaves the initial attack broadly intact";

        juce::String after;
        if (strongestOvershoot > 0.18 && strongestSettling > 0.10)
            after = "Overshoot with short settling / ringing";
        else if (strongestOvershoot > 0.12)
            after = "Brief overshoot detected";
        else if (strongestSettling > 0.10)
            after = "Short settling / ringing detected";
        else
            after = "No obvious overshoot or short ringing in the sine onset";

        return {
            "TEMPORAL RESPONSE  " + response,
            "ATTACK HANDLING  " + attack,
            "AFTER-EFFECTS  " + after
        };
    }

    struct ToneSummary
    {
        juce::String bass;
        juce::String midrange;
        juce::String treble;
        juce::String largest;
    };

    juce::String formatFrequency(double frequencyHz)
    {
        if (frequencyHz >= 1000.0)
            return juce::String(frequencyHz / 1000.0, frequencyHz >= 10000.0 ? 1 : 2) + " kHz";

        return juce::String(static_cast<int>(std::round(frequencyHz))) + " Hz";
    }

    juce::String formatSignedDb(double valueDb)
    {
        return juce::String(valueDb > 0.0 ? "+" : "")
             + juce::String(valueDb, 2) + " dB";
    }

    ToneSummary makeToneSummary(const MeasurementDataset& dataset)
    {
        const int runIdColumn = findColumn(dataset, "runId");
        const int frequencyColumn = findColumn(dataset, "freqHz");
        const int magnitudeColumn = findColumn(dataset, "magDb");

        if (runIdColumn < 0 || frequencyColumn < 0 || magnitudeColumn < 0 || dataset.rows.empty())
            return { "TONE SHAPE  unavailable", "LOW END  unavailable", "HIGH END  unavailable", "STRONGEST FEATURE  unavailable" };

        struct FrequencyPoint
        {
            double frequency = 0.0;
            double magnitude = 0.0;
        };

        std::map<int, std::vector<FrequencyPoint>> runs;
        const int lastColumn = std::max(runIdColumn, std::max(frequencyColumn, magnitudeColumn));

        for (const auto& row : dataset.rows)
        {
            if (static_cast<int>(row.size()) <= lastColumn)
                continue;

            const double frequency = row[static_cast<std::size_t>(frequencyColumn)];
            const double magnitude = row[static_cast<std::size_t>(magnitudeColumn)];

            if (!std::isfinite(frequency) || !std::isfinite(magnitude) || frequency < 20.0 || frequency > 20000.0)
                continue;

            const int runId = static_cast<int>(std::llround(row[static_cast<std::size_t>(runIdColumn)]));
            runs[runId].push_back({ frequency, magnitude });
        }

        // Thirty-one logarithmic bands give approximately third-octave spacing from
        // 20 Hz to 20 kHz. Each perceptual region therefore contributes equally,
        // instead of dense high-frequency FFT bins dominating the average.
        constexpr int bandCount = 31;
        const double minFrequency = 20.0;
        const double maxFrequency = 20000.0;
        const double logRange = std::log(maxFrequency / minFrequency);

        std::vector<double> bandSums(bandCount, 0.0);
        std::vector<int> bandContributions(bandCount, 0);

        for (const auto& [runId, points] : runs)
        {
            juce::ignoreUnused(runId);

            std::vector<double> runSums(bandCount, 0.0);
            std::vector<int> runCounts(bandCount, 0);

            for (const auto& point : points)
            {
                const double normalised = std::log(point.frequency / minFrequency) / logRange;
                const int band = juce::jlimit(0, bandCount - 1,
                    static_cast<int>(std::floor(normalised * static_cast<double>(bandCount))));

                runSums[static_cast<std::size_t>(band)] += point.magnitude;
                ++runCounts[static_cast<std::size_t>(band)];
            }

            std::vector<double> runBands(bandCount, std::numeric_limits<double>::quiet_NaN());
            for (int band = 0; band < bandCount; ++band)
            {
                if (runCounts[static_cast<std::size_t>(band)] > 0)
                    runBands[static_cast<std::size_t>(band)] =
                        runSums[static_cast<std::size_t>(band)] / runCounts[static_cast<std::size_t>(band)];
            }

            // Reference each run to its central midrange (roughly 500 Hz-2 kHz).
            // This removes simple broadband gain while preserving shelves, tilts,
            // filters and resonant features.
            double referenceSum = 0.0;
            int referenceCount = 0;
            for (int band = 0; band < bandCount; ++band)
            {
                const double centre = minFrequency * std::exp((band + 0.5) * logRange / bandCount);
                const double value = runBands[static_cast<std::size_t>(band)];
                if (centre >= 500.0 && centre <= 2000.0 && std::isfinite(value))
                {
                    referenceSum += value;
                    ++referenceCount;
                }
            }

            if (referenceCount == 0)
                continue;

            const double reference = referenceSum / static_cast<double>(referenceCount);
            for (int band = 0; band < bandCount; ++band)
            {
                const double value = runBands[static_cast<std::size_t>(band)];
                if (!std::isfinite(value))
                    continue;

                bandSums[static_cast<std::size_t>(band)] += value - reference;
                ++bandContributions[static_cast<std::size_t>(band)];
            }
        }

        std::vector<double> response(bandCount, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> centres(bandCount, 0.0);
        for (int band = 0; band < bandCount; ++band)
        {
            centres[static_cast<std::size_t>(band)] =
                minFrequency * std::exp((band + 0.5) * logRange / bandCount);

            if (bandContributions[static_cast<std::size_t>(band)] > 0)
                response[static_cast<std::size_t>(band)] =
                    bandSums[static_cast<std::size_t>(band)] /
                    static_cast<double>(bandContributions[static_cast<std::size_t>(band)]);
        }

        // Light three-band smoothing rejects isolated FFT-bin noise while retaining
        // broad shelves, tilts, filters and bell-shaped changes.
        std::vector<double> smoothed = response;
        for (int band = 1; band < bandCount - 1; ++band)
        {
            const double left = response[static_cast<std::size_t>(band - 1)];
            const double centre = response[static_cast<std::size_t>(band)];
            const double right = response[static_cast<std::size_t>(band + 1)];
            if (std::isfinite(left) && std::isfinite(centre) && std::isfinite(right))
                smoothed[static_cast<std::size_t>(band)] = (left + 2.0 * centre + right) / 4.0;
        }

        auto averageRange = [&](double lowHz, double highHz) -> double
        {
            double sum = 0.0;
            int count = 0;
            for (int band = 0; band < bandCount; ++band)
            {
                const double frequency = centres[static_cast<std::size_t>(band)];
                const double value = smoothed[static_cast<std::size_t>(band)];
                if (frequency >= lowHz && frequency <= highHz && std::isfinite(value))
                {
                    sum += value;
                    ++count;
                }
            }
            return count > 0 ? sum / static_cast<double>(count)
                             : std::numeric_limits<double>::quiet_NaN();
        };

        const double subBass = averageRange(20.0, 60.0);
        const double bass = averageRange(60.0, 180.0);
        const double mid = averageRange(500.0, 2000.0);
        const double presence = averageRange(2000.0, 6000.0);
        const double treble = averageRange(6000.0, 12000.0);
        const double air = averageRange(12000.0, 20000.0);

        if (!std::isfinite(bass) || !std::isfinite(mid) || !std::isfinite(treble))
            return { "TONE SHAPE  unavailable", "LOW END  unavailable", "HIGH END  unavailable", "STRONGEST FEATURE  unavailable" };

        double strongestValue = 0.0;
        double strongestFrequency = 0.0;
        int strongestBand = -1;
        for (int band = 0; band < bandCount; ++band)
        {
            const double value = smoothed[static_cast<std::size_t>(band)];
            if (std::isfinite(value) && (strongestBand < 0 || std::abs(value) > std::abs(strongestValue)))
            {
                strongestValue = value;
                strongestFrequency = centres[static_cast<std::size_t>(band)];
                strongestBand = band;
            }
        }

        const double lowAverage = std::isfinite(subBass) ? (subBass + bass) * 0.5 : bass;
        const double highAverage = std::isfinite(air) ? (treble + air) * 0.5 : treble;
        const double lowVsMid = lowAverage - mid;
        const double highVsMid = highAverage - mid;
        const double tiltDifference = lowAverage - highAverage;

        juce::String shape = "Neutral / flat";

        // Strong edge attenuation takes priority because it describes filter-like
        // behaviour more usefully than a generic shelf or tilt label.
        if (std::isfinite(air) && air <= -3.0 && treble < -1.0)
            shape = "Low-pass behaviour";
        else if (std::isfinite(subBass) && subBass <= -3.0 && bass < -1.0)
            shape = "High-pass behaviour";
        else
        {
            // Detect a broad interior peak or dip before shelves/tilts.
            double interiorPeak = 0.0;
            double interiorPeakFrequency = 0.0;
            for (int band = 3; band < bandCount - 3; ++band)
            {
                const double value = smoothed[static_cast<std::size_t>(band)];
                if (!std::isfinite(value))
                    continue;

                const double left = smoothed[static_cast<std::size_t>(band - 3)];
                const double right = smoothed[static_cast<std::size_t>(band + 3)];
                if (!std::isfinite(left) || !std::isfinite(right))
                    continue;

                const double prominence = value - (left + right) * 0.5;
                if (std::abs(prominence) > std::abs(interiorPeak))
                {
                    interiorPeak = prominence;
                    interiorPeakFrequency = centres[static_cast<std::size_t>(band)];
                }
            }

            if (interiorPeak >= 1.25)
                shape = "Bell boost around " + formatFrequency(interiorPeakFrequency);
            else if (interiorPeak <= -1.25)
                shape = "Bell cut around " + formatFrequency(interiorPeakFrequency);
            else if (tiltDifference >= 1.5 && lowVsMid > 0.4 && highVsMid < -0.4)
                shape = "Downward tilt";
            else if (tiltDifference <= -1.5 && lowVsMid < -0.4 && highVsMid > 0.4)
                shape = "Upward tilt";
            else if (lowVsMid >= 0.75 && std::abs(highVsMid) < 0.75)
                shape = "Low shelf lift";
            else if (lowVsMid <= -0.75 && std::abs(highVsMid) < 0.75)
                shape = "Low shelf cut";
            else if (highVsMid >= 0.75 && std::abs(lowVsMid) < 0.75)
                shape = "High shelf lift";
            else if (highVsMid <= -0.75 && std::abs(lowVsMid) < 0.75)
                shape = "High shelf cut";
            else if (std::abs(strongestValue) >= 0.5)
                shape = "Subtle shaped response";
        }

        auto describeRegion = [](double value, const juce::String& positive, const juce::String& negative) -> juce::String
        {
            if (value >= 1.0) return positive;
            if (value >= 0.35) return "Slight " + positive.toLowerCase();
            if (value <= -1.0) return negative;
            if (value <= -0.35) return "Slight " + negative.toLowerCase();
            return "Neutral";
        };

        const juce::String lowWords = describeRegion(lowVsMid, "Lift", "Reduction");
        const juce::String highWords = describeRegion(highVsMid, "Lift", "Roll-off");

        return {
            "TONE SHAPE  " + shape,
            "LOW END  " + lowWords + "  (" + formatSignedDb(lowVsMid) + " vs mids)",
            "HIGH END  " + highWords + "  (" + formatSignedDb(highVsMid) + " vs mids)",
            "STRONGEST FEATURE  " + formatSignedDb(strongestValue) + " @ " + formatFrequency(strongestFrequency)
        };
    }


    struct NonlinearitySummary
    {
        juce::String type;
        juce::String detail;
    };

    double medianValue(std::vector<double> values)
    {
        if (values.empty())
            return std::numeric_limits<double>::quiet_NaN();

        std::sort(values.begin(), values.end());
        const std::size_t middle = values.size() / 2;
        if ((values.size() % 2) == 0)
            return 0.5 * (values[middle - 1] + values[middle]);

        return values[middle];
    }

    NonlinearitySummary makeNonlinearitySummary(const MeasurementDataset& dataset)
    {
        const int runColumn = findColumn(dataset, "runId");
        const int xColumn = findColumn(dataset, "x");
        const int yColumn = findColumn(dataset, "meanY");
        const int gainColumn = findColumn(dataset, "inputGainDb");

        if (runColumn < 0 || xColumn < 0 || yColumn < 0 || dataset.rows.empty())
            return { "NONLINEARITY  unavailable", "CURVE BEHAVIOUR  unavailable" };

        struct Point { double x = 0.0; double y = 0.0; };
        struct RunCurve
        {
            std::vector<Point> points;
            double inputGainDb = 0.0;
            std::vector<std::pair<std::string, double>> parameters;
        };

        std::map<int, RunCurve> curves;
        const int requiredColumn = std::max(runColumn, std::max(xColumn, yColumn));

        for (const auto& row : dataset.rows)
        {
            if (static_cast<int>(row.size()) <= requiredColumn)
                continue;

            const double x = row[static_cast<std::size_t>(xColumn)];
            const double y = row[static_cast<std::size_t>(yColumn)];
            if (!std::isfinite(x) || !std::isfinite(y))
                continue;

            const int runId = static_cast<int>(std::llround(row[static_cast<std::size_t>(runColumn)]));
            auto& curve = curves[runId];
            curve.points.push_back({ x, y });

            if (gainColumn >= 0 && static_cast<int>(row.size()) > gainColumn)
                curve.inputGainDb = row[static_cast<std::size_t>(gainColumn)];

            if (curve.parameters.empty())
            {
                for (std::size_t column = 0; column < dataset.columns.size() && column < row.size(); ++column)
                {
                    const auto& name = dataset.columns[column];
                    if (name == "runId" || name == "binIndex" || name == "x" || name == "meanY"
                        || name == "count" || name == "inputGainDb")
                        continue;

                    curve.parameters.push_back({ name, row[column] });
                }
            }
        }

        struct Metrics
        {
            double score = -1.0;
            double deviation = 0.0;
            double edgeRatio = 1.0;
            double asymmetry = 0.0;
            double onsetFraction = 1.0;
            bool nonMonotonic = false;
            const RunCurve* curve = nullptr;
        };

        Metrics strongest;

        for (auto& [runId, curve] : curves)
        {
            juce::ignoreUnused(runId);
            if (curve.points.size() < 12)
                continue;

            std::sort(curve.points.begin(), curve.points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });

            double maxAbsX = 0.0;
            double maxAbsY = 0.0;
            for (const auto& point : curve.points)
            {
                maxAbsX = std::max(maxAbsX, std::abs(point.x));
                maxAbsY = std::max(maxAbsY, std::abs(point.y));
            }

            if (maxAbsX <= 1.0e-9 || maxAbsY <= 1.0e-9)
                continue;

            // Fit the small-signal centre region. This becomes the straight-line
            // reference against which bending, clipping and asymmetry are judged.
            double xy = 0.0;
            double xx = 0.0;
            for (const auto& point : curve.points)
            {
                const double fraction = std::abs(point.x) / maxAbsX;
                if (fraction >= 0.08 && fraction <= 0.42)
                {
                    xy += point.x * point.y;
                    xx += point.x * point.x;
                }
            }

            if (xx <= 1.0e-12)
                continue;

            const double centreGain = xy / xx;
            if (std::abs(centreGain) <= 1.0e-9)
                continue;

            double errorSq = 0.0;
            double referenceSq = 0.0;
            for (const auto& point : curve.points)
            {
                const double expected = centreGain * point.x;
                const double error = point.y - expected;
                errorSq += error * error;
                referenceSq += expected * expected;
            }
            const double deviation = std::sqrt(errorSq / std::max(referenceSq, 1.0e-12));

            std::vector<double> centreSlopes;
            std::vector<double> edgeSlopes;
            std::vector<std::pair<double, double>> positiveSlopes;
            bool nonMonotonic = false;

            for (std::size_t index = 1; index < curve.points.size(); ++index)
            {
                const auto& previous = curve.points[index - 1];
                const auto& current = curve.points[index];
                if (previous.x < 0.0 || current.x <= 0.0)
                    continue;

                const double dx = current.x - previous.x;
                if (dx <= 1.0e-12)
                    continue;

                const double slope = (current.y - previous.y) / dx;
                const double fraction = 0.5 * (current.x + previous.x) / maxAbsX;
                const double normalisedSlope = slope / centreGain;
                positiveSlopes.push_back({ fraction, normalisedSlope });

                if (fraction >= 0.12 && fraction <= 0.45)
                    centreSlopes.push_back(normalisedSlope);
                if (fraction >= 0.75)
                    edgeSlopes.push_back(normalisedSlope);

                if (normalisedSlope < -0.08)
                    nonMonotonic = true;
            }

            const double centreSlope = medianValue(centreSlopes);
            const double edgeSlope = medianValue(edgeSlopes);
            if (!std::isfinite(centreSlope) || !std::isfinite(edgeSlope) || std::abs(centreSlope) < 1.0e-6)
                continue;

            const double edgeRatio = edgeSlope / centreSlope;
            double onsetFraction = 1.0;
            for (const auto& [fraction, slope] : positiveSlopes)
            {
                if (fraction >= 0.25 && slope < 0.82 * centreSlope)
                {
                    onsetFraction = fraction;
                    break;
                }
            }

            // Compare positive and negative halves at matching amplitudes. An odd,
            // symmetrical transfer curve has y(+x) + y(-x) approximately zero.
            double asymmetrySum = 0.0;
            int asymmetryCount = 0;
            for (const auto& positive : curve.points)
            {
                if (positive.x <= 0.1 * maxAbsX)
                    continue;

                const Point* closest = nullptr;
                double closestDistance = std::numeric_limits<double>::max();
                for (const auto& negative : curve.points)
                {
                    if (negative.x >= 0.0)
                        break;
                    const double distance = std::abs(negative.x + positive.x);
                    if (distance < closestDistance)
                    {
                        closestDistance = distance;
                        closest = &negative;
                    }
                }

                if (closest != nullptr && closestDistance < 0.03 * maxAbsX)
                {
                    asymmetrySum += std::abs(positive.y + closest->y) / (2.0 * maxAbsY);
                    ++asymmetryCount;
                }
            }

            const double asymmetry = asymmetryCount > 0
                ? asymmetrySum / static_cast<double>(asymmetryCount)
                : 0.0;

            const double compressionAmount = std::max(0.0, 1.0 - edgeRatio);
            const double expansionAmount = std::max(0.0, edgeRatio - 1.0);
            const double score = deviation + 0.55 * compressionAmount
                               + 0.35 * expansionAmount + 0.65 * asymmetry
                               + (nonMonotonic ? 0.75 : 0.0);

            if (score > strongest.score)
            {
                strongest = { score, deviation, edgeRatio, asymmetry,
                              onsetFraction, nonMonotonic, &curve };
            }
        }

        if (strongest.curve == nullptr)
            return { "NONLINEARITY  unavailable", "CURVE BEHAVIOUR  insufficient transfer data" };

        juce::String type;
        if (strongest.deviation < 0.015
            && strongest.edgeRatio >= 0.88 && strongest.edgeRatio <= 1.12
            && strongest.asymmetry < 0.025 && !strongest.nonMonotonic)
        {
            type = "Mostly linear";
        }
        else if (strongest.nonMonotonic)
        {
            type = "Complex waveshaping";
        }
        else if (strongest.asymmetry >= 0.055)
        {
            type = "Asymmetric waveshaping";
        }
        else if (strongest.edgeRatio <= 0.16 && strongest.onsetFraction >= 0.62)
        {
            type = "Hard clipping";
        }
        else if (strongest.edgeRatio < 0.78)
        {
            type = "Progressive soft saturation";
        }
        else if (strongest.edgeRatio > 1.22)
        {
            type = "Expanding waveshaper";
        }
        else
        {
            type = "Gentle waveshaping";
        }

        const juce::String symmetry = strongest.asymmetry >= 0.055 ? "asymmetric" : "mostly symmetrical";
        juce::String onset;
        if (strongest.deviation < 0.015)
            onset = "no meaningful bend";
        else if (strongest.onsetFraction >= 0.68)
            onset = "abrupt near the top";
        else if (strongest.onsetFraction >= 0.42)
            onset = "progressive onset";
        else
            onset = "early onset";

        juce::String testedAt = "strongest at " + juce::String(strongest.curve->inputGainDb, 0) + " dBFS";
        if (strongest.curve->parameters.size() == 1)
        {
            testedAt += ", " + juce::String(strongest.curve->parameters[0].first)
                     + " " + juce::String(strongest.curve->parameters[0].second, 2);
        }

        return {
            "NONLINEARITY  " + type,
            "CURVE BEHAVIOUR  " + symmetry + " · " + onset + " · " + testedAt
        };
    }


    struct BehaviourSummary
    {
        juce::String headline;
        juce::String changes;
    };

    std::string makeParameterKey(const MeasurementDataset& dataset,
                                 const std::vector<double>& row,
                                 const std::vector<std::string>& excludedColumns)
    {
        std::string key;
        for (std::size_t column = 0; column < dataset.columns.size() && column < row.size(); ++column)
        {
            const auto& name = dataset.columns[column];
            if (std::find(excludedColumns.begin(), excludedColumns.end(), name) != excludedColumns.end())
                continue;

            key += name + "=" + std::to_string(row[column]) + ";";
        }
        return key;
    }

    double matchedMeanDifference(const std::map<std::string, double>& reference,
                                 const std::map<std::string, double>& comparison)
    {
        double sum = 0.0;
        int count = 0;
        for (const auto& [key, referenceValue] : reference)
        {
            const auto it = comparison.find(key);
            if (it == comparison.end())
                continue;

            sum += std::abs(it->second - referenceValue);
            ++count;
        }
        return count > 0 ? sum / static_cast<double>(count) : 0.0;
    }

    BehaviourSummary makeBehaviourSummary(const MeasurementDataset* rmsDataset,
                                           const MeasurementDataset* thdDataset,
                                           const MeasurementDataset* linearDataset)
    {
        double levelChange = 0.0;
        double dynamicChange = 0.0;
        double harmonicChangePercent = 0.0;
        double toneChange = 0.0;
        bool hasLevel = false;
        bool hasDynamics = false;
        bool hasHarmonics = false;
        bool hasTone = false;

        if (rmsDataset != nullptr)
        {
            const int gainColumn = findColumn(*rmsDataset, "inputGainDb");
            const int rmsIn = findColumn(*rmsDataset, "rmsInL");
            const int rmsOut = findColumn(*rmsDataset, "rmsOutL");
            const int peakIn = findColumn(*rmsDataset, "peakInL");
            const int peakOut = findColumn(*rmsDataset, "peakOutL");
            if (gainColumn >= 0 && rmsIn >= 0 && rmsOut >= 0 && peakIn >= 0 && peakOut >= 0)
            {
                std::map<std::string, double> gainAtReference, gainAtPush;
                std::map<std::string, double> crestAtReference, crestAtPush;
                const std::vector<std::string> excluded = {
                    "runId", "inputGainDb", "rmsInL", "rmsInR", "rmsOutL", "rmsOutR",
                    "peakInL", "peakInR", "peakOutL", "peakOutR"
                };

                for (const auto& row : rmsDataset->rows)
                {
                    const int last = std::max({ gainColumn, rmsIn, rmsOut, peakIn, peakOut });
                    if (static_cast<int>(row.size()) <= last)
                        continue;
                    const double inputGain = row[static_cast<std::size_t>(gainColumn)];
                    const double inRms = row[static_cast<std::size_t>(rmsIn)];
                    const double outRms = row[static_cast<std::size_t>(rmsOut)];
                    const double inPeak = row[static_cast<std::size_t>(peakIn)];
                    const double outPeak = row[static_cast<std::size_t>(peakOut)];
                    if (inRms <= 0.0 || outRms <= 0.0 || inPeak <= 0.0 || outPeak <= 0.0)
                        continue;

                    const auto key = makeParameterKey(*rmsDataset, row, excluded);
                    const double gain = amplitudeToDb(outRms) - amplitudeToDb(inRms);
                    const double crest = (amplitudeToDb(outPeak) - amplitudeToDb(inPeak)) - gain;
                    if (std::abs(inputGain + 18.0) < 0.25)
                    {
                        gainAtReference[key] = gain;
                        crestAtReference[key] = crest;
                    }
                    else if (std::abs(inputGain + 12.0) < 0.25)
                    {
                        gainAtPush[key] = gain;
                        crestAtPush[key] = crest;
                    }
                }

                if (!gainAtReference.empty() && !gainAtPush.empty())
                {
                    levelChange = matchedMeanDifference(gainAtReference, gainAtPush);
                    dynamicChange = matchedMeanDifference(crestAtReference, crestAtPush);
                    hasLevel = true;
                    hasDynamics = true;
                }
            }
        }

        if (thdDataset != nullptr)
        {
            const int gainColumn = findColumn(*thdDataset, "inputGainDb");
            const int thdColumn = findColumn(*thdDataset, "thd");
            if (gainColumn >= 0 && thdColumn >= 0)
            {
                struct Average { double sum = 0.0; int count = 0; };
                std::map<std::string, Average> referenceAverages, pushAverages;
                const std::vector<std::string> excluded = { "runId", "centreSample", "thd", "inputGainDb" };

                for (const auto& row : thdDataset->rows)
                {
                    const int last = std::max(gainColumn, thdColumn);
                    if (static_cast<int>(row.size()) <= last)
                        continue;
                    const double inputGain = row[static_cast<std::size_t>(gainColumn)];
                    const double thd = row[static_cast<std::size_t>(thdColumn)];
                    if (!std::isfinite(thd) || thd < 0.0)
                        continue;
                    const auto key = makeParameterKey(*thdDataset, row, excluded);
                    auto* target = std::abs(inputGain + 18.0) < 0.25 ? &referenceAverages
                                 : std::abs(inputGain + 12.0) < 0.25 ? &pushAverages : nullptr;
                    if (target != nullptr)
                    {
                        (*target)[key].sum += thd * 100.0;
                        ++(*target)[key].count;
                    }
                }

                std::map<std::string, double> reference, push;
                for (const auto& [key, value] : referenceAverages)
                    if (value.count > 0) reference[key] = value.sum / value.count;
                for (const auto& [key, value] : pushAverages)
                    if (value.count > 0) push[key] = value.sum / value.count;

                if (!reference.empty() && !push.empty())
                {
                    harmonicChangePercent = matchedMeanDifference(reference, push);
                    hasHarmonics = true;
                }
            }
        }

        if (linearDataset != nullptr)
        {
            const int runColumn = findColumn(*linearDataset, "runId");
            const int gainColumn = findColumn(*linearDataset, "inputGainDb");
            const int frequencyColumn = findColumn(*linearDataset, "freqHz");
            const int magnitudeColumn = findColumn(*linearDataset, "magDb");
            if (runColumn >= 0 && gainColumn >= 0 && frequencyColumn >= 0 && magnitudeColumn >= 0)
            {
                struct Profile
                {
                    double inputGain = 0.0;
                    std::string key;
                    std::vector<double> sum = std::vector<double>(16, 0.0);
                    std::vector<int> count = std::vector<int>(16, 0);
                };
                std::map<int, Profile> profiles;
                const std::vector<std::string> excluded = { "runId", "freqHz", "magDb", "inputGainDb" };
                const double minF = 20.0, maxF = 20000.0;
                const double logRange = std::log(maxF / minF);

                for (const auto& row : linearDataset->rows)
                {
                    const int last = std::max({ runColumn, gainColumn, frequencyColumn, magnitudeColumn });
                    if (static_cast<int>(row.size()) <= last)
                        continue;
                    const double frequency = row[static_cast<std::size_t>(frequencyColumn)];
                    const double magnitude = row[static_cast<std::size_t>(magnitudeColumn)];
                    if (frequency < minF || frequency > maxF || !std::isfinite(magnitude))
                        continue;
                    const int runId = static_cast<int>(std::llround(row[static_cast<std::size_t>(runColumn)]));
                    auto& profile = profiles[runId];
                    profile.inputGain = row[static_cast<std::size_t>(gainColumn)];
                    profile.key = makeParameterKey(*linearDataset, row, excluded);
                    const int band = juce::jlimit(0, 15, static_cast<int>(std::floor(
                        std::log(frequency / minF) / logRange * 16.0)));
                    profile.sum[static_cast<std::size_t>(band)] += magnitude;
                    ++profile.count[static_cast<std::size_t>(band)];
                }

                std::map<std::string, std::vector<double>> reference, push;
                for (auto& [runId, profile] : profiles)
                {
                    juce::ignoreUnused(runId);
                    std::vector<double> bands(16, 0.0);
                    double midSum = 0.0; int midCount = 0;
                    for (int band = 0; band < 16; ++band)
                    {
                        if (profile.count[static_cast<std::size_t>(band)] > 0)
                            bands[static_cast<std::size_t>(band)] = profile.sum[static_cast<std::size_t>(band)] /
                                profile.count[static_cast<std::size_t>(band)];
                        const double centre = minF * std::exp((band + 0.5) * logRange / 16.0);
                        if (centre >= 500.0 && centre <= 2000.0 && profile.count[static_cast<std::size_t>(band)] > 0)
                        { midSum += bands[static_cast<std::size_t>(band)]; ++midCount; }
                    }
                    const double mid = midCount > 0 ? midSum / midCount : 0.0;
                    for (auto& value : bands) value -= mid;
                    if (std::abs(profile.inputGain + 18.0) < 0.25) reference[profile.key] = bands;
                    else if (std::abs(profile.inputGain + 12.0) < 0.25) push[profile.key] = bands;
                }

                double sum = 0.0; int count = 0;
                for (const auto& [key, refBands] : reference)
                {
                    const auto it = push.find(key);
                    if (it == push.end()) continue;
                    for (std::size_t band = 0; band < refBands.size(); ++band)
                    { sum += std::abs(it->second[band] - refBands[band]); ++count; }
                }
                if (count > 0)
                {
                    toneChange = sum / count;
                    hasTone = true;
                }
            }
        }

        double score = 0.0;
        if (hasLevel) score = std::max(score, levelChange / 1.5);
        if (hasDynamics) score = std::max(score, dynamicChange / 1.0);
        if (hasHarmonics) score = std::max(score, harmonicChangePercent / 1.0);
        if (hasTone) score = std::max(score, toneChange / 1.5);

        juce::String response;
        if (score < 0.20) response = "Mostly static";
        else if (score < 0.55) response = "Mildly responsive to drive";
        else if (score < 1.20) response = "Responsive to drive";
        else response = "Strongly responsive to drive";

        juce::StringArray changed;
        if (hasHarmonics && harmonicChangePercent >= 0.10) changed.add("harmonics");
        if (hasDynamics && dynamicChange >= 0.20) changed.add("dynamics");
        if (hasTone && toneChange >= 0.20) changed.add("tone");
        if (hasLevel && levelChange >= 0.35) changed.add("level");

        const juce::String detail = changed.isEmpty()
            ? "CHANGES WITH DRIVE  Little meaningful change from -18 to -12 dBFS"
            : "CHANGES WITH DRIVE  " + changed.joinIntoString(" · ");

        return { "BEHAVIOUR  " + response + "  (-18 to -12 dBFS)", detail };
    }

    struct OperatingRangeSummary
    {
        juce::String range;
        juce::String onset;
        juce::String stress;
    };

    struct DriveStepMetrics
    {
        double level = 0.0;
        double crest = 0.0;
        double thdPercent = 0.0;
        double tone = 0.0;
        bool hasLevel = false;
        bool hasCrest = false;
        bool hasThd = false;
        bool hasTone = false;
    };

    OperatingRangeSummary makeOperatingRangeSummary(const MeasurementDataset* rmsDataset,
                                                     const MeasurementDataset* thdDataset,
                                                     const MeasurementDataset* linearDataset)
    {
        const std::vector<double> levels = { -18.0, -12.0, -6.0, -3.0 };
        std::map<double, DriveStepMetrics> metrics;
        for (double level : levels)
            metrics[level] = {};

        if (rmsDataset != nullptr)
        {
            const int gainColumn = findColumn(*rmsDataset, "inputGainDb");
            const int rmsIn = findColumn(*rmsDataset, "rmsInL");
            const int rmsOut = findColumn(*rmsDataset, "rmsOutL");
            const int peakIn = findColumn(*rmsDataset, "peakInL");
            const int peakOut = findColumn(*rmsDataset, "peakOutL");
            if (gainColumn >= 0 && rmsIn >= 0 && rmsOut >= 0 && peakIn >= 0 && peakOut >= 0)
            {
                struct Average { double gain = 0.0; double crest = 0.0; int count = 0; };
                std::map<double, Average> averages;
                for (const auto& row : rmsDataset->rows)
                {
                    const int last = std::max({ gainColumn, rmsIn, rmsOut, peakIn, peakOut });
                    if (static_cast<int>(row.size()) <= last)
                        continue;
                    const double inputGain = row[static_cast<std::size_t>(gainColumn)];
                    const double inRms = row[static_cast<std::size_t>(rmsIn)];
                    const double outRms = row[static_cast<std::size_t>(rmsOut)];
                    const double inPeak = row[static_cast<std::size_t>(peakIn)];
                    const double outPeak = row[static_cast<std::size_t>(peakOut)];
                    if (inRms <= 0.0 || outRms <= 0.0 || inPeak <= 0.0 || outPeak <= 0.0)
                        continue;
                    for (double level : levels)
                    {
                        if (std::abs(inputGain - level) < 0.25)
                        {
                            const double gain = amplitudeToDb(outRms) - amplitudeToDb(inRms);
                            const double crest = (amplitudeToDb(outPeak) - amplitudeToDb(inPeak)) - gain;
                            averages[level].gain += gain;
                            averages[level].crest += crest;
                            ++averages[level].count;
                            break;
                        }
                    }
                }
                for (const auto& [level, average] : averages)
                {
                    if (average.count > 0)
                    {
                        metrics[level].level = average.gain / average.count;
                        metrics[level].crest = average.crest / average.count;
                        metrics[level].hasLevel = true;
                        metrics[level].hasCrest = true;
                    }
                }
            }
        }

        if (thdDataset != nullptr)
        {
            const int gainColumn = findColumn(*thdDataset, "inputGainDb");
            const int thdColumn = findColumn(*thdDataset, "thd");
            if (gainColumn >= 0 && thdColumn >= 0)
            {
                struct Average { double sum = 0.0; int count = 0; };
                std::map<double, Average> averages;
                for (const auto& row : thdDataset->rows)
                {
                    if (static_cast<int>(row.size()) <= std::max(gainColumn, thdColumn))
                        continue;
                    const double inputGain = row[static_cast<std::size_t>(gainColumn)];
                    const double thd = row[static_cast<std::size_t>(thdColumn)];
                    if (!std::isfinite(thd) || thd < 0.0)
                        continue;
                    for (double level : levels)
                    {
                        if (std::abs(inputGain - level) < 0.25)
                        {
                            averages[level].sum += thd * 100.0;
                            ++averages[level].count;
                            break;
                        }
                    }
                }
                for (const auto& [level, average] : averages)
                {
                    if (average.count > 0)
                    {
                        metrics[level].thdPercent = average.sum / average.count;
                        metrics[level].hasThd = true;
                    }
                }
            }
        }

        if (linearDataset != nullptr)
        {
            const int runColumn = findColumn(*linearDataset, "runId");
            const int gainColumn = findColumn(*linearDataset, "inputGainDb");
            const int frequencyColumn = findColumn(*linearDataset, "freqHz");
            const int magnitudeColumn = findColumn(*linearDataset, "magDb");
            if (runColumn >= 0 && gainColumn >= 0 && frequencyColumn >= 0 && magnitudeColumn >= 0)
            {
                struct Profile
                {
                    double inputGain = 0.0;
                    std::vector<double> sum = std::vector<double>(16, 0.0);
                    std::vector<int> count = std::vector<int>(16, 0);
                };
                std::map<int, Profile> profiles;
                const double minF = 20.0, maxF = 20000.0;
                const double logRange = std::log(maxF / minF);
                for (const auto& row : linearDataset->rows)
                {
                    const int last = std::max({ runColumn, gainColumn, frequencyColumn, magnitudeColumn });
                    if (static_cast<int>(row.size()) <= last)
                        continue;
                    const double frequency = row[static_cast<std::size_t>(frequencyColumn)];
                    const double magnitude = row[static_cast<std::size_t>(magnitudeColumn)];
                    if (frequency < minF || frequency > maxF || !std::isfinite(magnitude))
                        continue;
                    const int runId = static_cast<int>(std::llround(row[static_cast<std::size_t>(runColumn)]));
                    auto& profile = profiles[runId];
                    profile.inputGain = row[static_cast<std::size_t>(gainColumn)];
                    const int band = juce::jlimit(0, 15, static_cast<int>(std::floor(
                        std::log(frequency / minF) / logRange * 16.0)));
                    profile.sum[static_cast<std::size_t>(band)] += magnitude;
                    ++profile.count[static_cast<std::size_t>(band)];
                }

                struct ToneAverage { double sum = 0.0; int count = 0; };
                std::map<double, ToneAverage> toneAverages;
                for (const auto& [runId, profile] : profiles)
                {
                    juce::ignoreUnused(runId);
                    std::vector<double> bands(16, 0.0);
                    double midSum = 0.0; int midCount = 0;
                    for (int band = 0; band < 16; ++band)
                    {
                        if (profile.count[static_cast<std::size_t>(band)] > 0)
                            bands[static_cast<std::size_t>(band)] = profile.sum[static_cast<std::size_t>(band)] /
                                profile.count[static_cast<std::size_t>(band)];
                        const double centre = minF * std::exp((band + 0.5) * logRange / 16.0);
                        if (centre >= 500.0 && centre <= 2000.0 && profile.count[static_cast<std::size_t>(band)] > 0)
                        { midSum += bands[static_cast<std::size_t>(band)]; ++midCount; }
                    }
                    if (midCount == 0)
                        continue;
                    const double mid = midSum / midCount;
                    double shape = 0.0; int shapeCount = 0;
                    for (int band = 0; band < 16; ++band)
                    {
                        if (profile.count[static_cast<std::size_t>(band)] > 0)
                        { shape += std::abs(bands[static_cast<std::size_t>(band)] - mid); ++shapeCount; }
                    }
                    if (shapeCount == 0)
                        continue;
                    for (double level : levels)
                    {
                        if (std::abs(profile.inputGain - level) < 0.25)
                        {
                            toneAverages[level].sum += shape / shapeCount;
                            ++toneAverages[level].count;
                            break;
                        }
                    }
                }
                for (const auto& [level, average] : toneAverages)
                {
                    if (average.count > 0)
                    {
                        metrics[level].tone = average.sum / average.count;
                        metrics[level].hasTone = true;
                    }
                }
            }
        }

        struct StepResult
        {
            double score = 0.0;
            juce::StringArray changes;
        };
        std::vector<StepResult> steps;
        for (std::size_t i = 1; i < levels.size(); ++i)
        {
            const auto& previous = metrics[levels[i - 1]];
            const auto& current = metrics[levels[i]];
            StepResult step;
            if (previous.hasThd && current.hasThd)
            {
                const double delta = current.thdPercent - previous.thdPercent;
                step.score = std::max(step.score, std::abs(delta) / 1.0);
                if (delta >= 0.10) step.changes.add("harmonics rise");
            }
            if (previous.hasCrest && current.hasCrest)
            {
                const double delta = current.crest - previous.crest;
                step.score = std::max(step.score, std::abs(delta) / 1.0);
                if (delta <= -0.20) step.changes.add("transients soften");
                else if (delta >= 0.20) step.changes.add("transients open up");
            }
            if (previous.hasTone && current.hasTone)
            {
                const double delta = current.tone - previous.tone;
                step.score = std::max(step.score, std::abs(delta) / 1.5);
                if (std::abs(delta) >= 0.20) step.changes.add("tone shifts");
            }
            if (previous.hasLevel && current.hasLevel)
            {
                const double delta = current.level - previous.level;
                step.score = std::max(step.score, std::abs(delta) / 1.5);
                if (std::abs(delta) >= 0.35) step.changes.add("gain changes");
            }
            steps.push_back(step);
        }

        int firstMeaningful = -1;
        for (std::size_t i = 0; i < steps.size(); ++i)
        {
            if (steps[i].score >= 0.20)
            { firstMeaningful = static_cast<int>(i); break; }
        }

        juce::String range;
        juce::String onset;
        if (firstMeaningful < 0)
        {
            range = "OPERATING RANGE  Stable through -3 dBFS";
            onset = "CHARACTER STARTS  No major change detected";
        }
        else
        {
            const double stableTo = levels[static_cast<std::size_t>(firstMeaningful)];
            const double startsAt = levels[static_cast<std::size_t>(firstMeaningful + 1)];
            range = "OPERATING RANGE  Stable through " + juce::String(stableTo, 0) + " dBFS";
            onset = "CHARACTER STARTS  Around " + juce::String(startsAt, 0) + " dBFS";
        }

        juce::StringArray stressChanges;
        if (!steps.empty())
        {
            for (const auto& change : steps.back().changes)
                stressChanges.addIfNotAlreadyThere(change);
        }
        const juce::String stress = stressChanges.isEmpty()
            ? "UNDER STRESS  No major change from -6 to -3 dBFS"
            : "UNDER STRESS  " + stressChanges.joinIntoString(" · ");

        return { range, onset, stress };
    }

}


    struct IdentitySummary
    {
        juce::String role;
        juce::String why;
        juce::String watch;
    };

    bool containsIgnoreCase(const juce::String& text, const juce::String& token)
    {
        return text.toLowerCase().contains(token.toLowerCase());
    }


    struct HarmonicFingerprintEvidence
    {
        bool available = false;
        double lowColourDb = -120.0;
        double midColourDb = -120.0;
        double highColourDb = -120.0;
        double overallColourDb = -120.0;
        double highOrderVsLowOrderDb = 0.0;
        double secondVsThirdDb = 0.0;
    };

    HarmonicFingerprintEvidence makeHarmonicFingerprintEvidence(const MeasurementDataset& data)
    {
        HarmonicFingerprintEvidence evidence;
        if (data.isEmpty())
            return evidence;

        const int runCol = findColumn(data, "runId");
        const int freqCol = findColumn(data, "fundamentalHz");
        const int orderCol = findColumn(data, "harmonicOrder");
        const int magCol = findColumn(data, "magDbRelativeToInput");
        const int gainCol = findColumn(data, "inputGainDb");
        if (runCol < 0 || freqCol < 0 || orderCol < 0 || magCol < 0)
            return evidence;

        struct Accumulator
        {
            double lowPower = 0.0, midPower = 0.0, highPower = 0.0, allPower = 0.0;
            double lowOrderPower = 0.0, highOrderPower = 0.0;
            double secondPower = 0.0, thirdPower = 0.0;
            int lowCount = 0, midCount = 0, highCount = 0, allCount = 0;
            int lowOrderCount = 0, highOrderCount = 0, secondCount = 0, thirdCount = 0;
        };

        std::map<int, Accumulator> perRun;
        double closestGainDistance = std::numeric_limits<double>::max();
        double referenceGain = -18.0;
        if (gainCol >= 0)
        {
            for (const auto& row : data.rows)
            {
                if ((int) row.size() <= gainCol) continue;
                const double distance = std::abs(row[(std::size_t) gainCol] + 18.0);
                if (distance < closestGainDistance)
                {
                    closestGainDistance = distance;
                    referenceGain = row[(std::size_t) gainCol];
                }
            }
        }

        for (const auto& row : data.rows)
        {
            const int required = std::max({ runCol, freqCol, orderCol, magCol, gainCol });
            if ((int) row.size() <= required) continue;
            if (gainCol >= 0 && std::abs(row[(std::size_t) gainCol] - referenceGain) > 0.25)
                continue;

            const int order = (int) std::lround(row[(std::size_t) orderCol]);
            if (order < 2 || order > 7) continue;
            const double frequency = row[(std::size_t) freqCol];
            const double db = juce::jlimit(-120.0, 24.0, row[(std::size_t) magCol]);
            const double power = std::pow(10.0, db / 10.0);
            auto& a = perRun[(int) std::lround(row[(std::size_t) runCol])];

            a.allPower += power; ++a.allCount;
            if (frequency < 250.0) { a.lowPower += power; ++a.lowCount; }
            else if (frequency < 4000.0) { a.midPower += power; ++a.midCount; }
            else { a.highPower += power; ++a.highCount; }

            if (order <= 3) { a.lowOrderPower += power; ++a.lowOrderCount; }
            if (order >= 5) { a.highOrderPower += power; ++a.highOrderCount; }
            if (order == 2) { a.secondPower += power; ++a.secondCount; }
            if (order == 3) { a.thirdPower += power; ++a.thirdCount; }
        }

        const auto meanDb = [](double power, int count)
        {
            return count > 0 ? 10.0 * std::log10(std::max(power / (double) count, 1.0e-12)) : -120.0;
        };

        // Use the run with the strongest meaningful nonlinear contribution at the
        // reference gain. This describes the processor's available colour without
        // letting a bypass/neutral setting dilute the fingerprint.
        double strongest = -120.0;
        for (const auto& [runId, a] : perRun)
        {
            juce::ignoreUnused(runId);
            const double overall = meanDb(a.allPower, a.allCount);
            if (overall <= strongest) continue;
            strongest = overall;
            evidence.available = true;
            evidence.overallColourDb = overall;
            evidence.lowColourDb = meanDb(a.lowPower, a.lowCount);
            evidence.midColourDb = meanDb(a.midPower, a.midCount);
            evidence.highColourDb = meanDb(a.highPower, a.highCount);
            evidence.highOrderVsLowOrderDb = meanDb(a.highOrderPower, a.highOrderCount)
                                              - meanDb(a.lowOrderPower, a.lowOrderCount);
            evidence.secondVsThirdDb = meanDb(a.secondPower, a.secondCount)
                                        - meanDb(a.thirdPower, a.thirdCount);
        }
        return evidence;
    }

    struct ProducerVocabularySummary
    {
        juce::String character;
        juce::String behaviour;
    };

    ProducerVocabularySummary makeProducerVocabularySummary(const IdentitySummary& identity,
                                                              const LevelSummary& level,
                                                              const HarmonicSummary& harmonic,
                                                              const ToneSummary& tone,
                                                              const BehaviourSummary& drive,
                                                              const TemporalSummary& temporal,
                                                              const HarmonicFingerprintEvidence& fingerprint)
    {
        juce::StringArray characterWords;
        juce::StringArray behaviourWords;

        const auto addUnique = [](juce::StringArray& words, const juce::String& word)
        {
            if (!words.contains(word, true))
                words.add(word);
        };

        const auto all = identity.role + " " + identity.why + " " + harmonic.identity + " "
                       + harmonic.balance + " " + tone.bass + " " + tone.midrange + " "
                       + tone.treble + " " + drive.headline + " " + temporal.response + " "
                       + temporal.attack + " " + temporal.after + " " + level.level + " "
                       + level.peaks + " " + level.dynamics;

        struct WordScore { juce::String word; double score = 0.0; };
        std::vector<WordScore> scores {
            { "Sweet", 0.0 }, { "Warm", 0.0 }, { "Dense", 0.0 },
            { "Full", 0.0 }, { "Focused", 0.0 }, { "Airy", 0.0 },
            { "Grindy", 0.0 }, { "Smooth", 0.0 }, { "Present", 0.0 },
            { "Articulate", 0.0 }, { "Punchy", 0.0 }, { "Forward", 0.0 },
            { "Subtle", 0.0 }
        };
        const auto score = [&scores](const juce::String& word, double amount)
        {
            for (auto& item : scores)
                if (item.word == word) { item.score += amount; return; }
        };

        // Existing sine, tone and temporal evidence.
        if (all.containsIgnoreCase("2nd harmonic") || all.containsIgnoreCase("even-heavy"))
            score("Sweet", 70.0);
        if (all.containsIgnoreCase("dense") || all.containsIgnoreCase("high-order"))
            score("Dense", 55.0);
        if (all.containsIgnoreCase("preserves") || all.containsIgnoreCase("lifts the first attack"))
        {
            score("Present", 65.0);
            score("Articulate", 70.0);
        }
        if (all.containsIgnoreCase("punchier") || all.containsIgnoreCase("more open"))
            score("Punchy", 55.0);
        // Tone-summary text is only a fallback. When the H2-H7 fingerprint
        // exists, frequency-specific harmonic evidence makes these decisions.
        if (!fingerprint.available)
        {
            if (all.containsIgnoreCase("roll-off") || all.containsIgnoreCase("softens"))
                score("Smooth", 40.0);
            if (all.containsIgnoreCase("bass lift") || all.containsIgnoreCase("low end"))
                score("Full", 45.0);
            if (all.containsIgnoreCase("midrange forward") || all.containsIgnoreCase("mids forward"))
                score("Forward", 55.0);
            if (all.containsIgnoreCase("treble lift") || all.containsIgnoreCase("air")
                || all.containsIgnoreCase("high-frequency lift"))
                score("Airy", 55.0);
        }

        // Frequency-specific H2-H7 evidence from the sweep fingerprint. These rules
        // compare regions, so they generalise to any processor rather than any name.
        if (fingerprint.available)
        {
            const double lowVsMid = fingerprint.lowColourDb - fingerprint.midColourDb;
            const double highVsMid = fingerprint.highColourDb - fingerprint.midColourDb;
            const double highVsLow = fingerprint.highColourDb - fingerprint.lowColourDb;

            if (lowVsMid > 2.5) score("Warm", juce::jlimit(35.0, 90.0, 45.0 + lowVsMid * 2.0));
            // Full requires substantial colour, not merely a low-frequency bias.
            if (lowVsMid > 1.5 && fingerprint.overallColourDb > -30.0)
                score("Full", juce::jlimit(20.0, 58.0, 25.0 + lowVsMid * 0.8));
            if (lowVsMid > 3.0 && highVsMid > 1.5) score("Focused", 75.0);
            else if (lowVsMid > 5.0) score("Focused", 65.0);

            if (highVsLow > 2.5 && fingerprint.highOrderVsLowOrderDb < -3.0)
                score("Airy", 65.0);
            if (highVsMid > 2.5 && fingerprint.highOrderVsLowOrderDb > -5.0)
                score("Grindy", 75.0);
            if (fingerprint.highOrderVsLowOrderDb < -8.0)
                score("Smooth", 50.0);
            if (fingerprint.secondVsThirdDb > 3.0)
                score("Sweet", 35.0);
            if (fingerprint.overallColourDb > -32.0
                && std::abs(lowVsMid) < 8.0 && std::abs(highVsMid) < 8.0)
                score("Dense", 45.0);
            if (fingerprint.overallColourDb < -30.0)
                score("Subtle", 35.0);
        }

        if (all.containsIgnoreCase("largely unchanged")
            || all.containsIgnoreCase("no adjustment needed")
            || all.containsIgnoreCase("mostly static"))
            score("Subtle", 35.0);

        std::stable_sort(scores.begin(), scores.end(), [](const WordScore& a, const WordScore& b)
        {
            return a.score > b.score;
        });
        for (const auto& item : scores)
        {
            if (item.score < 50.0) continue;
            addUnique(characterWords, item.word);
            if (characterWords.size() == 3) break;
        }

        if (all.containsIgnoreCase("stateful") || all.containsIgnoreCase("signal memory"))
            addUnique(behaviourWords, "Stateful");
        if (all.containsIgnoreCase("responsive"))
            addUnique(behaviourWords, "Responsive");
        if (all.containsIgnoreCase("progressive"))
            addUnique(behaviourWords, "Progressive");
        if (all.containsIgnoreCase("harmonic") || all.containsIgnoreCase("saturat") || all.containsIgnoreCase("clip"))
            addUnique(behaviourWords, "Nonlinear");
        if (all.containsIgnoreCase("mostly static") || all.containsIgnoreCase("little meaningful change"))
            addUnique(behaviourWords, "Static");
        if (identity.role.containsIgnoreCase("clean gain"))
            addUnique(behaviourWords, "Transparent");

        while (characterWords.size() > 3)
            characterWords.remove(characterWords.size() - 1);
        while (behaviourWords.size() > 3)
            behaviourWords.remove(behaviourWords.size() - 1);

        const auto character = characterWords.isEmpty()
            ? juce::String("CHARACTER  Subtle")
            : juce::String("CHARACTER  ") + characterWords.joinIntoString(" · ");
        const auto behaviour = behaviourWords.isEmpty()
            ? juce::String("BEHAVIOUR  Static")
            : juce::String("BEHAVIOUR  ") + behaviourWords.joinIntoString(" · ");

        return { character, behaviour };
    }

    IdentitySummary makeIdentitySummary(const LevelSummary& level,
                                        const HarmonicSummary& harmonic,
                                        const ToneSummary& tone,
                                        const BehaviourSummary& behaviour,
                                        const OperatingRangeSummary& operating,
                                        const NonlinearitySummary& nonlinearity,
                                        const WaveformSummary& waveform,
                                        const TemporalSummary& temporal)
    {
        const bool harmonicAvailable = !containsIgnoreCase(harmonic.identity, "unavailable")
                                    && !containsIgnoreCase(harmonic.identity, "run thd");
        const bool toneAvailable = !containsIgnoreCase(tone.bass, "run linear")
                                && !containsIgnoreCase(tone.bass, "unavailable");
        const bool mostlySecond = containsIgnoreCase(harmonic.identity, "mostly 2nd");
        const bool mostlyThird = containsIgnoreCase(harmonic.identity, "mostly 3rd");
        const bool highOrder = containsIgnoreCase(harmonic.identity, "high-order")
                            || containsIgnoreCase(harmonic.growth, "higher harmonics");
        const bool evenHeavy = containsIgnoreCase(harmonic.balance, "even-heavy");
        const bool oddHeavy = containsIgnoreCase(harmonic.balance, "odd-heavy");
        const bool veryClean = containsIgnoreCase(harmonic.identity, "very clean");
        const bool toneNeutral = containsIgnoreCase(tone.bass, "neutral")
                              && containsIgnoreCase(tone.midrange, "neutral")
                              && containsIgnoreCase(tone.treble, "neutral");
        const bool toneChanged = toneAvailable && !toneNeutral;
        const bool dynamicsControlled = containsIgnoreCase(level.dynamics, "controlled")
                                     || containsIgnoreCase(level.peaks, "softer");
        const bool dynamicsOpen = containsIgnoreCase(level.dynamics, "more open")
                               || containsIgnoreCase(level.peaks, "punchier");
        const bool dynamicsStatic = containsIgnoreCase(level.dynamics, "unchanged")
                                 || containsIgnoreCase(level.peaks, "little change");
        const bool hardClip = containsIgnoreCase(waveform.shape, "hard clipping")
                           || containsIgnoreCase(nonlinearity.type, "hard clipping");
        const bool foldback = containsIgnoreCase(waveform.shape, "foldback");
        const bool softSaturation = containsIgnoreCase(waveform.shape, "smooth")
                                 || containsIgnoreCase(nonlinearity.type, "soft saturation")
                                 || containsIgnoreCase(nonlinearity.type, "complex waveshaping");
        const bool responsive = containsIgnoreCase(behaviour.headline, "responsive");
        const bool stronglyResponsive = containsIgnoreCase(behaviour.headline, "strongly responsive");
        const bool staticBehaviour = containsIgnoreCase(behaviour.headline, "static");
        const bool stateful = containsIgnoreCase(temporal.response, "stateful")
                           || containsIgnoreCase(temporal.response, "signal memory");
        const bool attackPreserved = containsIgnoreCase(temporal.attack, "preserves")
                                  || containsIgnoreCase(temporal.attack, "lifts");
        const bool attackSoftened = containsIgnoreCase(temporal.attack, "softens");
        const bool overshoot = containsIgnoreCase(temporal.after, "overshoot");
        const bool ringing = containsIgnoreCase(temporal.after, "ringing");

        struct Scores
        {
            int clean = 0;
            int harmonicEnhancer = 0;
            int saturator = 0;
            int transientPreservingSaturator = 0;
            int transientSmoother = 0;
            int transientEnhancer = 0;
            int toneShaper = 0;
            int clipper = 0;
            int wavefolder = 0;
        } score;

        // Harmonic evidence carries the most weight when the plugin is primarily
        // creating colour. Small crest changes cannot overrule a clear harmonic identity.
        if (mostlySecond) { score.harmonicEnhancer += 60; score.saturator += 15; }
        if (mostlyThird)  { score.saturator += 35; }
        if (highOrder)    { score.saturator += 25; score.clipper += 12; score.wavefolder += 10; }
        if (evenHeavy)    { score.harmonicEnhancer += 18; score.saturator += 8; }
        if (oddHeavy)     { score.saturator += 15; score.clipper += 8; }
        if (veryClean)    { score.clean += 55; }

        // Dynamics evidence describes the audible consequence, but is deliberately
        // weaker than a dominant harmonic signature.
        if (dynamicsControlled) { score.saturator += 24; score.transientSmoother += 32; }
        if (dynamicsOpen)       { score.transientEnhancer += 22; score.transientPreservingSaturator += 16; }
        if (dynamicsStatic)     { score.clean += 8; score.harmonicEnhancer += 12; }

        // Temporal evidence is the key distinction for processors such as Spiral2.
        if (stateful)          { score.transientPreservingSaturator += 55; score.saturator += 10; }
        if (attackPreserved)   { score.transientPreservingSaturator += 45; score.transientEnhancer += 18; }
        if (attackSoftened)    { score.transientSmoother += 35; score.saturator += 10; }
        if (overshoot)         { score.transientPreservingSaturator += 12; score.transientEnhancer += 8; }
        if (ringing)           { score.transientPreservingSaturator += 6; }

        if (toneChanged)       { score.toneShaper += 45; }
        if (toneNeutral)       { score.clean += 6; score.harmonicEnhancer += 8; score.saturator += 5; }
        if (responsive)        { score.saturator += 16; score.transientPreservingSaturator += 8; }
        if (stronglyResponsive){ score.saturator += 10; }
        if (staticBehaviour)   { score.clean += 10; score.harmonicEnhancer += 8; }
        if (softSaturation)    { score.saturator += 20; }
        if (hardClip)          { score.clipper += 70; }

        // Foldback is supporting evidence only. It becomes the headline only when
        // sustained reversal, hard high-order content and no stronger temporal or
        // harmonic role all agree.
        if (foldback)          { score.wavefolder += 18; }
        if (foldback && highOrder) score.wavefolder += 18;
        if (foldback && hardClip)  score.wavefolder += 14;
        if (stateful || mostlySecond || dynamicsControlled)
            score.wavefolder = std::max(0, score.wavefolder - 22);

        if (toneNeutral && dynamicsStatic && staticBehaviour && veryClean)
            score.clean += 35;
        if (mostlySecond && dynamicsStatic)
            score.harmonicEnhancer += 25;
        if (dynamicsControlled && harmonicAvailable)
            score.saturator += 18;

        struct Candidate { int score; juce::String role; juce::String why; };
        std::vector<Candidate> candidates {
            { score.clean, "PRIMARY ROLE  Clean gain / utility",
              "WHY  Stays neutral, adds little harmonic colour and barely changes dynamics" },
            { score.harmonicEnhancer, "PRIMARY ROLE  Harmonic enhancer",
              mostlySecond
                ? "WHY  Adds mainly 2nd harmonic colour with little change to tone or transients"
                : "WHY  Adds harmonic colour while leaving the original signal largely intact" },
            { score.transientPreservingSaturator, "PRIMARY ROLE  Transient-preserving saturator",
              "WHY  Adds dense harmonic colour while preserving or lifting onset attacks" },
            { score.saturator, "PRIMARY ROLE  Saturation / bus control",
              "WHY  Harmonic colour grows while peaks become smoother, denser or more controlled" },
            { score.transientSmoother, "PRIMARY ROLE  Transient smoother",
              "WHY  Softens peak contrast more than it changes tone or harmonic character" },
            { score.transientEnhancer, "PRIMARY ROLE  Transient enhancer",
              "WHY  Preserves or increases attack contrast with little peak smoothing" },
            { score.toneShaper, "PRIMARY ROLE  Tone shaper",
              "WHY  Its clearest contribution is a measurable change in frequency balance" },
            { score.clipper, "PRIMARY ROLE  Clipper / peak cutter",
              "WHY  Peaks are cut abruptly and higher odd harmonics build strongly" },
            { score.wavefolder, "PRIMARY ROLE  Creative wavefolder distortion",
              "WHY  Sustained transfer reversal and dense high-order harmonics dominate the result" }
        };

        const auto best = *std::max_element(candidates.begin(), candidates.end(),
            [](const Candidate& a, const Candidate& b) { return a.score < b.score; });

        juce::String watch = "WATCH OUT  No major issue detected in the tested range";
        if (containsIgnoreCase(level.level, "suggested output trim"))
            watch = "WATCH OUT  Output level changes significantly; level-match before judging";
        if (containsIgnoreCase(operating.stress, "highs soften") || containsIgnoreCase(tone.treble, "roll-off"))
            watch = "WATCH OUT  Top end softens as the processor is driven";
        if (containsIgnoreCase(operating.stress, "harmonics rise") || highOrder)
            watch = "WATCH OUT  Higher harmonics increase when pushed hard";
        if (stateful && attackPreserved && (overshoot || ringing))
            watch = "WATCH OUT  Strong settings can exaggerate attack, air or short overshoot";
        if (hardClip)
            watch = "WATCH OUT  Peak flattening becomes abrupt at stronger settings";
        if (best.role.containsIgnoreCase("wavefolder"))
            watch = "WATCH OUT  Extreme settings create sustained foldback and dense high-order distortion";

        return { best.role, best.why, watch };
    }


    juce::var makeNumericStats(double minimum, double maximum, double mean, std::int64_t count)
    {
        auto* object = new juce::DynamicObject();
        object->setProperty("min", minimum);
        object->setProperty("max", maximum);
        object->setProperty("mean", mean);
        object->setProperty("count", static_cast<juce::int64>(count));
        return juce::var(object);
    }


    bool isEvidenceSettingColumn(const juce::String& name, const Config& config)
    {
        // Evidence settings are deliberately strict: only the drive bucket and
        // parameters that PluginDNA intentionally changed for this run belong
        // here. Every analyser-produced column is a measurement, even when a
        // future analyser introduces a metric name we have never seen before.
        if (name == "inputGainDb")
            return true;
        for (const auto& bucket : config.parameterBuckets)
            if (name == bucket.paramName)
                return true;
        return false;
    }

    juce::var deriveRunEvidenceMetrics(const juce::String& filename,
                                       const MeasurementDataset& dataset,
                                       const std::vector<const std::vector<double>*>& rows)
    {
        auto* derived = new juce::DynamicObject();
        const auto index = [&dataset](const char* name) { return findColumn(dataset, name); };
        const auto value = [](const std::vector<double>& row, int column)
        {
            return column >= 0 && static_cast<std::size_t>(column) < row.size()
                ? row[static_cast<std::size_t>(column)]
                : std::numeric_limits<double>::quiet_NaN();
        };
        const auto db = [](double amplitude)
        {
            return 20.0 * std::log10(std::max(std::abs(amplitude), 1.0e-15));
        };
        const auto collect = [&](const char* name, const std::function<bool(const std::vector<double>&)>& include = {})
        {
            std::vector<double> values;
            const int c = index(name);
            if (c < 0) return values;
            for (const auto* row : rows)
            {
                if (include && !include(*row)) continue;
                const double v = value(*row, c);
                if (std::isfinite(v)) values.push_back(v);
            }
            return values;
        };
        const auto mean = [](const std::vector<double>& values)
        {
            if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
            return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
        };
        const auto stdev = [&](const std::vector<double>& values)
        {
            if (values.size() < 2) return 0.0;
            const double m = mean(values);
            double sum = 0.0;
            for (double v : values) sum += (v - m) * (v - m);
            return std::sqrt(sum / static_cast<double>(values.size()));
        };
        const auto add = [derived](const juce::Identifier& name, double v)
        {
            if (std::isfinite(v)) derived->setProperty(name, v);
        };
        const auto addStats = [&](const juce::Identifier& name, const std::vector<double>& values)
        {
            if (values.empty()) return;
            const auto mm = std::minmax_element(values.begin(), values.end());
            auto* stats = new juce::DynamicObject();
            stats->setProperty("min", *mm.first);
            stats->setProperty("max", *mm.second);
            stats->setProperty("mean", mean(values));
            stats->setProperty("standardDeviation", stdev(values));
            stats->setProperty("count", static_cast<juce::int64>(values.size()));
            derived->setProperty(name, juce::var(stats));
        };

        if (filename.containsIgnoreCase("rms_peak"))
        {
            const int ri=index("rmsInL"), ro=index("rmsOutL"), pi=index("peakInL"), po=index("peakOutL");
            std::vector<double> level, peak, crestIn, crestOut, crestChange;
            for (const auto* row : rows)
            {
                const double a=value(*row,ri), b=value(*row,ro), c=value(*row,pi), d=value(*row,po);
                if (!(a>0&&b>0&&c>0&&d>0)) continue;
                const double li=db(b)-db(a), pk=db(d)-db(c), ci=db(c)-db(a), co=db(d)-db(b);
                level.push_back(li); peak.push_back(pk); crestIn.push_back(ci); crestOut.push_back(co); crestChange.push_back(co-ci);
            }
            addStats("levelChangeDb", level); addStats("peakChangeDb", peak);
            addStats("crestFactorInDb", crestIn); addStats("crestFactorOutDb", crestOut);
            addStats("crestFactorChangeDb", crestChange);
            if (!level.empty() && !peak.empty())
            {
                add("peakMinusLevelChangeDb", mean(peak)-mean(level));
                add("averageGainDb", mean(level));
                add("gainVarianceDb", stdev(level));
                add("peakGainVarianceDb", stdev(peak));
                add("dynamicHeadroomChangeDb", -mean(peak));
                const double delta = mean(peak) - mean(level);
                add("compressionExpansionIndexDb", delta);
                add("compressionAmountDb", std::max(0.0, -delta));
                add("expansionAmountDb", std::max(0.0, delta));
                add("levelTrackingErrorDb", std::abs(mean(level)));
                add("peakTrackingErrorDb", std::abs(mean(peak)));
            }
        }
        else if (filename.containsIgnoreCase("grid_thd"))
        {
            const auto thd=collect("thd"); addStats("thd", thd);
            double evenPower=0.0, oddPower=0.0;
            std::vector<double> ladder;
            for (int h=2; h<=10; ++h)
            {
                const auto vals=collect(("h"+std::to_string(h)).c_str());
                if (vals.empty()) continue;
                const double m=mean(vals); add("h"+juce::String(h), m); ladder.push_back(std::max(m,1.0e-15));
                const double pwr=m*m; if ((h%2)==0) evenPower+=pwr; else oddPower+=pwr;
            }
            if (evenPower>0 || oddPower>0)
            {
                const double oddVsEvenDb = 10.0*std::log10(std::max(oddPower,1.0e-30)/std::max(evenPower,1.0e-30));
                add("oddEvenPowerRatioDb", oddVsEvenDb);
                add("evenOddPowerRatioDb", -oddVsEvenDb);
            }
            if (ladder.size()>=2)
            {
                std::vector<double> slopes;
                double weightedOrder = 0.0, totalAmplitude = 0.0, entropy = 0.0;
                int significantCount = 0, highestSignificant = 0;
                for (std::size_t i=1;i<ladder.size();++i) slopes.push_back(db(ladder[i])-db(ladder[i-1]));
                for (std::size_t i=0;i<ladder.size();++i)
                {
                    const double a = ladder[i];
                    totalAmplitude += a;
                    weightedOrder += static_cast<double>(i + 2) * a;
                    if (db(a) > -100.0) { ++significantCount; highestSignificant = static_cast<int>(i + 2); }
                }
                if (totalAmplitude > 0.0)
                {
                    for (double a : ladder) { const double p = a / totalAmplitude; if (p > 0.0) entropy -= p * std::log2(p); }
                    add("harmonicCentroidOrder", weightedOrder / totalAmplitude);
                    add("harmonicEntropyBits", entropy);
                }
                add("harmonicLadderMeanStepDb", mean(slopes));
                add("harmonicDecaySlopeDbPerOrder", mean(slopes));
                add("harmonicLadderStepVariationDb", stdev(slopes));
                add("harmonicLadderSmoothness", 1.0 / (1.0 + stdev(slopes)));
                add("harmonicDensity", static_cast<double>(significantCount) / static_cast<double>(ladder.size()));
                add("highestSignificantHarmonic", highestSignificant);
                add("ladderContinuity", static_cast<double>(significantCount) / std::max(1, highestSignificant - 1));
            }
        }
        else if (filename.containsIgnoreCase("transfer_curves"))
        {
            const int xc=index("x"), yc=index("meanY");
            std::vector<std::pair<double,double>> points;
            for(const auto* row:rows){const double x=value(*row,xc),y=value(*row,yc);if(std::isfinite(x)&&std::isfinite(y))points.push_back({x,y});}
            std::sort(points.begin(),points.end());
            if(points.size()>=3)
            {
                double maxDeviation=0.0, posErr=0.0, negErr=0.0; int posN=0,negN=0, reversals=0, inflections=0;
                std::vector<double> slopes;
                for(const auto& pt:points){maxDeviation=std::max(maxDeviation,std::abs(pt.second-pt.first));if(pt.first>=0){posErr+=pt.second-pt.first;++posN;}else{negErr+=pt.second-pt.first;++negN;}}
                for(std::size_t i=1;i<points.size();++i){const double dx=points[i].first-points[i-1].first;if(std::abs(dx)>1e-12)slopes.push_back((points[i].second-points[i-1].second)/dx);}
                for(std::size_t i=1;i<slopes.size();++i){if(slopes[i]<0&&slopes[i-1]>=0)++reversals;if((slopes[i]-slopes[i-1])*(i>1?slopes[i-1]-slopes[i-2]:0.0)<0)++inflections;}
                const double meanLocalSlope = mean(slopes);
                add("maxStaticDeviation",maxDeviation); add("meanSlope",meanLocalSlope); add("slopeVariation",stdev(slopes));
                add("localSlopeMean", meanLocalSlope);
                add("localSlopeMinimum",slopes.empty()?NAN:*std::min_element(slopes.begin(),slopes.end()));
                add("localSlopeMaximum",slopes.empty()?NAN:*std::max_element(slopes.begin(),slopes.end()));
                add("minimumLocalSlope",slopes.empty()?NAN:*std::min_element(slopes.begin(),slopes.end()));
                add("maximumLocalSlope",slopes.empty()?NAN:*std::max_element(slopes.begin(),slopes.end()));
                double localSlopeAtZero = std::numeric_limits<double>::quiet_NaN();
                double closestZero = std::numeric_limits<double>::max();
                for (std::size_t i=0; i<slopes.size(); ++i)
                {
                    const double midpoint = 0.5 * (points[i].first + points[i+1].first);
                    if (std::abs(midpoint) < closestZero) { closestZero = std::abs(midpoint); localSlopeAtZero = slopes[i]; }
                }
                add("localSlopeAtZero", localSlopeAtZero);
                const double asymmetry=(posN?posErr/posN:0.0)+(negN?negErr/negN:0.0);
                const double symmetryIndex = 1.0/(1.0+std::abs(asymmetry));
                add("asymmetryError",asymmetry); add("symmetryScore",symmetryIndex); add("symmetryIndex",symmetryIndex);
                double deviationPower = 0.0, gainErrorPower = 0.0;
                for (const auto& pt : points)
                {
                    const double nonlinearDeviation = pt.second - meanLocalSlope * pt.first;
                    deviationPower += nonlinearDeviation * nonlinearDeviation;
                    const double gainError = (meanLocalSlope - 1.0) * pt.first;
                    gainErrorPower += gainError * gainError;
                }
                const double nonlinearRms = std::sqrt(deviationPower / points.size());
                const double gainErrorRms = std::sqrt(gainErrorPower / points.size());
                add("saturationEfficiency", nonlinearRms / std::max(1.0e-15, nonlinearRms + gainErrorRms));
                add("foldbackReversalCount",reversals); add("inflectionCount",inflections);
                int nearLinear=0, compressing=0, expanding=0, dead=0; double kneeInput=std::numeric_limits<double>::quiet_NaN();
                for(std::size_t i=0;i<slopes.size();++i)
                {
                    const double sl=slopes[i];
                    if(std::abs(sl-1.0)<0.05) ++nearLinear;
                    if(sl<0.95) ++compressing;
                    if(sl>1.05) ++expanding;
                    if(std::abs(sl)<0.05) ++dead;
                    if(!std::isfinite(kneeInput) && std::abs(sl-1.0)>0.10) kneeInput=std::abs(points[i].first);
                }
                if(!slopes.empty())
                {
                    add("linearRegionFraction",static_cast<double>(nearLinear)/slopes.size());
                    add("compressionRegionFraction",static_cast<double>(compressing)/slopes.size());
                    add("expansionRegionFraction",static_cast<double>(expanding)/slopes.size());
                    add("deadZoneFraction",static_cast<double>(dead)/slopes.size());
                    add("meanCompressionRatio", mean(slopes)>0.0 ? 1.0/mean(slopes) : 0.0);
                }
                add("estimatedKneeInput",kneeInput);
            }
        }
        else if (filename.containsIgnoreCase("linear_response"))
        {
            const int fc=index("freqHz"), mc=index("magDb");
            struct Band{double lo,hi;const char* name;};
            const Band bands[]={{20,60,"subDb"},{60,200,"bassDb"},{200,800,"lowMidDb"},{800,3000,"midDb"},{3000,7000,"presenceDb"},{7000,14000,"highDb"},{14000,24000,"airDb"}};
            std::vector<std::pair<double,double>> response;
            for(const auto* row:rows){const double f=value(*row,fc),m=value(*row,mc);if(f>0&&std::isfinite(m))response.push_back({f,m});}
            for(const auto& band:bands){std::vector<double> vals;for(const auto& p:response)if(p.first>=band.lo&&p.first<band.hi)vals.push_back(p.second);add(band.name,mean(vals));}
            if(!response.empty())
            {
                auto maxIt=std::max_element(response.begin(),response.end(),[](auto&a,auto&b){return a.second<b.second;});
                auto minIt=std::min_element(response.begin(),response.end(),[](auto&a,auto&b){return a.second<b.second;});
                add("largestBoostDb",maxIt->second);add("largestBoostFrequencyHz",maxIt->first);add("largestCutDb",minIt->second);add("largestCutFrequencyHz",minIt->first);
                std::vector<double> mags;for(auto&p:response)mags.push_back(p.second);add("responseRippleDb",*std::max_element(mags.begin(),mags.end())-*std::min_element(mags.begin(),mags.end()));
                std::sort(response.begin(),response.end());
                const double logSpan=std::log2(response.back().first/response.front().first);
                if(logSpan>0.0) add("overallTiltDbPerOctave",(response.back().second-response.front().second)/logSpan);
                double weightSum=0.0, centroid=0.0, spread=0.0;
                for(const auto& p:response){const double w=std::pow(10.0,p.second/20.0);weightSum+=w;centroid+=p.first*w;}
                if(weightSum>0.0){centroid/=weightSum;for(const auto& p:response){const double w=std::pow(10.0,p.second/20.0);spread+=(p.first-centroid)*(p.first-centroid)*w;}add("spectralCentroidHz",centroid);add("spectralSpreadHz",std::sqrt(spread/weightSum));}
                int resonances=0; double strongestProminence=0.0, strongestQ=0.0;
                for(std::size_t i=1;i+1<response.size();++i)
                {
                    const double prominence=response[i].second-0.5*(response[i-1].second+response[i+1].second);
                    if(prominence>0.25){++resonances;strongestProminence=std::max(strongestProminence,prominence);const double bw=std::max(1.0,response[i+1].first-response[i-1].first);strongestQ=std::max(strongestQ,response[i].first/bw);}
                }
                add("resonanceCount",resonances);add("strongestResonanceProminenceDb",strongestProminence);add("strongestResonanceQEstimate",strongestQ);
            }
        }
        else if (filename.containsIgnoreCase("harmonic_fingerprint"))
        {
            const int fc=index("fundamentalHz"), hc=index("harmonicOrder"), mc=index("magDbRelativeToInput");
            double lowP=0,highP=0; std::vector<double> low,mid,high;
            for(int h=2;h<=7;++h)
            {
                std::vector<double> vals;
                for(const auto* row:rows)if(std::llround(value(*row,hc))==h){const double m=value(*row,mc),f=value(*row,fc);if(std::isfinite(m)){vals.push_back(m);if(f<250)low.push_back(m);else if(f<3000)mid.push_back(m);else high.push_back(m);}}
                addStats("h"+juce::String(h)+"AcrossFrequencyDb",vals);
                for(double m:vals){double pwr=std::pow(10.0,m/10.0);if(h<=3)lowP+=pwr;else highP+=pwr;}
            }
            add("lowBandHarmonicDb",mean(low));add("midBandHarmonicDb",mean(mid));add("highBandHarmonicDb",mean(high));
            if(lowP>0)add("highOrderVsLowOrderDb",10.0*std::log10(std::max(highP,1e-30)/lowP));
        }
        else if (filename.containsIgnoreCase("residual_dna"))
        {
            addStats("overallResidualDb",collect("residualDbRelativeToOutput")); addStats("lowResidualDb",collect("lowResidualDb"));
            addStats("midResidualDb",collect("midResidualDb")); addStats("highResidualDb",collect("highResidualDb")); addStats("ultrasonicResidualDb",collect("ultrasonicResidualDb"));
            addStats("alignmentLagSamples",collect("lagSamples")); addStats("alignmentGainDb",collect("gainDb"));
            const auto overall=collect("residualDbRelativeToOutput"), low=collect("lowResidualDb"), mid=collect("midResidualDb"), high=collect("highResidualDb"), ultra=collect("ultrasonicResidualDb");
            add("lowVsMidResidualDb",mean(low)-mean(mid));
            add("highVsMidResidualDb",mean(high)-mean(mid));
            add("ultrasonicVsHighResidualDb",mean(ultra)-mean(high));
            if(!overall.empty()) add("residualDynamicVariationDb",*std::max_element(overall.begin(),overall.end())-*std::min_element(overall.begin(),overall.end()));
            const double la=std::pow(10.0,mean(low)/20.0), ma=std::pow(10.0,mean(mid)/20.0), ha=std::pow(10.0,mean(high)/20.0), ua=std::pow(10.0,mean(ultra)/20.0);
            const double total=la+ma+ha+ua;
            if(total>0.0){add("residualSpectralCentroidHz",(100.0*la+1000.0*ma+8000.0*ha+24000.0*ua)/total);add("residualLowEnergyFraction",la/total);add("residualMidEnergyFraction",ma/total);add("residualHighEnergyFraction",ha/total);add("residualUltrasonicEnergyFraction",ua/total);}
            addStats("residualCentroidHz", collect("residualSpectralCentroidHz"));
            addStats("residualBandwidthHz", collect("residualEffectiveBandwidthHz"));
            addStats("residualEntropyBits", collect("residualSpectralEntropyBits"));
            addStats("residualEntropyNormalized", collect("residualSpectralEntropyNormalized"));
            add("residualSpectralTiltDb",mean(high)-mean(low));
        }
        else if (filename.containsIgnoreCase("interaction_dna"))
        {
            const int cc=index("productClass"); auto all=collect("magDbRelativeToInput"); addStats("allProductsDb",all);
            for(int cls=1;cls<=4;++cls)addStats("productClass"+juce::String(cls)+"Db",collect("magDbRelativeToInput",[&](const std::vector<double>& r){return std::llround(value(r,cc))==cls;}));
            if(!all.empty())
            {
                add("worstProductDb",*std::max_element(all.begin(),all.end()));
                add("productDensityAboveMinus80Db",static_cast<double>(std::count_if(all.begin(),all.end(),[](double v){return v>-80.0;}))/all.size());
                add("productDensityAboveMinus100Db",static_cast<double>(std::count_if(all.begin(),all.end(),[](double v){return v>-100.0;}))/all.size());
                add("interactionComplexity",stdev(all));
                double sumPower=0.0;for(double v:all)sumPower+=std::pow(10.0,v/10.0);add("totalInteractionPowerDb",10.0*std::log10(std::max(sumPower,1.0e-30)));
            }
        }
        else if (filename.containsIgnoreCase("timing_dna"))
        {
            const int rc=index("recordType"),fc=index("frequencyHz"),pc=index("phaseDegrees");
            const auto summary=[&](const std::vector<double>& r){return std::llround(value(r,rc))==2;};
            addStats("latencySamples",collect("latencySamples",summary));addStats("latencyMs",collect("latencyMs",summary));addStats("preRingingDb",collect("preRingingDb",summary));
            addStats("postRingingDb",collect("postRingingDb",summary));addStats("settlingMs",collect("settlingMs",summary));addStats("overshootDb",collect("overshootDb",summary));
            std::vector<std::pair<double,double>> phase;for(const auto* row:rows)if(std::llround(value(*row,rc))==1){double f=value(*row,fc),p=value(*row,pc);if(f>0&&std::isfinite(p))phase.push_back({f,p});}
            std::sort(phase.begin(),phase.end());std::vector<double> gd, phaseValues, lowPhase, midPhase, highPhase, lowGd, midGd, highGd;
            for(const auto& p:phase){phaseValues.push_back(p.second);if(p.first<200)lowPhase.push_back(p.second);else if(p.first<3000)midPhase.push_back(p.second);else highPhase.push_back(p.second);}
            for(std::size_t i=1;i<phase.size();++i){double df=phase[i].first-phase[i-1].first;if(df<=0)continue;double dp=phase[i].second-phase[i-1].second;while(dp>180)dp-=360;while(dp<-180)dp+=360;const double g=-dp/360.0/df*1000.0;gd.push_back(g);const double f=0.5*(phase[i].first+phase[i-1].first);if(f<200)lowGd.push_back(g);else if(f<3000)midGd.push_back(g);else highGd.push_back(g);}
            addStats("phaseRotationDegrees",phaseValues);addStats("lowPhaseDegrees",lowPhase);addStats("midPhaseDegrees",midPhase);addStats("highPhaseDegrees",highPhase);
            addStats("groupDelayMs",gd);addStats("lowGroupDelayMs",lowGd);addStats("midGroupDelayMs",midGd);addStats("highGroupDelayMs",highGd);
            if(!phaseValues.empty()){add("maximumAbsolutePhaseRotationDegrees",std::abs(*std::max_element(phaseValues.begin(),phaseValues.end(),[](double a,double b){return std::abs(a)<std::abs(b);})));add("phaseRotationVariance",stdev(phaseValues));}
        }
        else if (filename.containsIgnoreCase("alias_stress"))
        {
            const auto aliasDb = collect("aliasPowerDb");
            const auto aliasRatio = collect("aliasPowerRatio");
            const auto legalRatio = collect("legalHarmonicPowerRatio");
            const auto folds = collect("detectedFoldCount");
            const auto predicted = collect("predictedFoldCount");
            const auto strongestFreq = collect("strongestAliasFrequencyHz");
            const auto strongestRatio = collect("strongestAliasRatio");
            const auto integrity = collect("ladderIntegrity");
            const auto clean = collect("cleanAtMinus90Db");
            addStats("aliasPowerDb", aliasDb);
            addStats("aliasPowerRatio", aliasRatio);
            addStats("legalHarmonicPowerRatio", legalRatio);
            addStats("detectedFoldCount", folds);
            addStats("predictedFoldCount", predicted);
            addStats("strongestAliasFrequencyHz", strongestFreq);
            addStats("strongestAliasRatio", strongestRatio);
            addStats("ladderIntegrity", integrity);
            if (!aliasDb.empty()) add("worstAliasDb", *std::max_element(aliasDb.begin(), aliasDb.end()));
            if (!folds.empty()) add("maximumDetectedFoldCount", *std::max_element(folds.begin(), folds.end()));
            if (!integrity.empty()) add("minimumLadderIntegrity", *std::min_element(integrity.begin(), integrity.end()));

            const int fc = index("fundamentalHz"), ac = index("aliasPowerDb");
            double highestClean = std::numeric_limits<double>::quiet_NaN();
            double firstDirty = std::numeric_limits<double>::quiet_NaN();
            for (const auto* row : rows)
            {
                const double f = value(*row, fc), a = value(*row, ac);
                if (!(f > 0.0) || !std::isfinite(a)) continue;
                if (a <= -90.0) highestClean = std::isfinite(highestClean) ? std::max(highestClean, f) : f;
                else firstDirty = std::isfinite(firstDirty) ? std::min(firstDirty, f) : f;
            }
            add("highestCleanFundamentalHzAtMinus90Db", highestClean);
            add("firstAliasFundamentalHzAboveMinus90Db", firstDirty);
            if (!clean.empty()) add("cleanTestFractionAtMinus90Db", mean(clean));
        }
        else if (filename.containsIgnoreCase("stereo_dna"))
        {
            addStats("leftOnlyLeakDb", collect("leftOnlyLeakDb"));
            addStats("rightOnlyLeakDb", collect("rightOnlyLeakDb"));
            addStats("midChannelMismatchDb", collect("midChannelMismatchDb"));
            addStats("sideChannelMismatchDb", collect("sideChannelMismatchDb"));
            addStats("midToSideLeakDb", collect("midToSideLeakDb"));
            addStats("sideToMidLeakDb", collect("sideToMidLeakDb"));
            addStats("midCorrelation", collect("midCorrelation"));
            addStats("sideCorrelation", collect("sideCorrelation"));
            const auto lr = collect("leftOnlyLeakDb"), rl = collect("rightOnlyLeakDb");
            if (!lr.empty() && !rl.empty()) add("crosstalkAsymmetryDb", std::abs(mean(lr) - mean(rl)));
            const auto mm = collect("midChannelMismatchDb"), sm = collect("sideChannelMismatchDb");
            if (!mm.empty() && !sm.empty()) add("maximumChannelMismatchDb", std::max(*std::max_element(mm.begin(), mm.end()), *std::max_element(sm.begin(), sm.end())));
        }
        else if (filename.containsIgnoreCase("summing_dna"))
        {
            addStats("separateRms", collect("separateRms"));
            addStats("summedRms", collect("summedRms"));
            addStats("interactionResidualRms", collect("interactionResidualRms"));
            addStats("interactionResidualDbRelative", collect("interactionResidualDbRelative"));
            addStats("summedCrestDb", collect("summedCrestDb"));
            addStats("nonAdditivityPercent", collect("nonAdditivityPercent"));
            const auto residual = collect("interactionResidualDbRelative");
            if (!residual.empty()) add("strongestSummingInteractionDb", *std::max_element(residual.begin(), residual.end()));
        }
        else if (filename.containsIgnoreCase("boundary_dna"))
        {
            const int cc=index("componentClass"), gc=index("inputGainDb"), fc=index("frequencyHz");
            const auto products=collect("outputDbRelativeToInputRms");
            addStats("dcOffset",collect("dcOffset"));addStats("transferDb",collect("transferDb"));addStats("phaseDegrees",collect("phaseDegrees"));addStats("boundaryProductsDb",products);
            for(int cls=0;cls<=8;++cls)addStats("componentClass"+juce::String(cls)+"Db",collect("outputDbRelativeToInputRms",[&](const std::vector<double>& r){return std::llround(value(r,cc))==cls;}));
            if(!products.empty()){add("worstBoundaryProductDb",*std::max_element(products.begin(),products.end()));add("boundaryProductDensityAboveMinus80Db",static_cast<double>(std::count_if(products.begin(),products.end(),[](double v){return v>-80.0;}))/products.size());}
            std::vector<std::pair<double,double>> byGain;
            for(const auto* row:rows){const double g=value(*row,gc),v=value(*row,index("outputDbRelativeToInputRms"));if(std::isfinite(g)&&std::isfinite(v))byGain.push_back({g,v});}
            std::sort(byGain.begin(),byGain.end());double onset=std::numeric_limits<double>::quiet_NaN();for(const auto& p:byGain)if(p.second>-80.0){onset=p.first;break;}add("firstInstabilityInputDb",onset);
            const auto freqs=collect("frequencyHz");if(!freqs.empty()){add("lowestBoundaryFrequencyHz",*std::min_element(freqs.begin(),freqs.end()));add("highestBoundaryFrequencyHz",*std::max_element(freqs.begin(),freqs.end()));}
        }
        return juce::var(derived);
    }

    juce::var summariseDatasetForEvidence(const juce::String& filename, const MeasurementDataset& dataset, const Config& config)
    {
        auto* result = new juce::DynamicObject();
        result->setProperty("source", filename);
        result->setProperty("analyser", juce::String(dataset.analyserName));
        result->setProperty("rowCount", static_cast<juce::int64>(dataset.rows.size()));

        juce::Array<juce::var> columnNames;
        for (const auto& column : dataset.columns)
            columnNames.add(juce::String(column));
        result->setProperty("columns", columnNames);

        // Compact Evidence v2: dataset-wide bookkeeping lives in the raw package.
        // The compact pack keeps settings plus derived measurements only.

        const int runIdColumn = findColumn(dataset, "runId");
        if (runIdColumn >= 0)
        {
            std::map<int, std::vector<const std::vector<double>*>> groupedRows;
            for (const auto& row : dataset.rows)
                if (static_cast<std::size_t>(runIdColumn) < row.size())
                    groupedRows[static_cast<int>(std::llround(row[static_cast<std::size_t>(runIdColumn)]))].push_back(&row);

            juce::Array<juce::var> runs;
            for (const auto& group : groupedRows)
            {
                const auto runId = group.first;
                const auto& rows = group.second;
                auto* runObject = new juce::DynamicObject();
                runObject->setProperty("runId", runId);
                runObject->setProperty("rowCount", static_cast<juce::int64>(rows.size()));
                auto* constants = new juce::DynamicObject();
                auto* metrics = new juce::DynamicObject();

                for (std::size_t column = 0; column < dataset.columns.size(); ++column)
                {
                    if (static_cast<int>(column) == runIdColumn)
                        continue;
                    double minimum = std::numeric_limits<double>::infinity();
                    double maximum = -std::numeric_limits<double>::infinity();
                    double sum = 0.0;
                    std::int64_t count = 0;
                    for (const auto* row : rows)
                    {
                        if (column >= row->size() || !std::isfinite((*row)[column]))
                            continue;
                        minimum = std::min(minimum, (*row)[column]);
                        maximum = std::max(maximum, (*row)[column]);
                        sum += (*row)[column];
                        ++count;
                    }
                    if (count == 0)
                        continue;
                    const auto name = juce::String(dataset.columns[column]);
                    if (isEvidenceSettingColumn(name, config))
                    {
                        if (std::abs(maximum - minimum) < 1.0e-12)
                            constants->setProperty(name, minimum);
                        else
                            constants->setProperty(name, makeNumericStats(minimum, maximum, sum / static_cast<double>(count), count));
                    }
                    else
                    {
                        // Raw analyser columns are intentionally omitted from the compact pack.
                        // deriveRunEvidenceMetrics() below retains the musically/engineering-useful
                        // fingerprint while Export Raw Data remains lossless.
                    }
                }

                metrics->setProperty("derived", deriveRunEvidenceMetrics(filename, dataset, rows));
                runObject->setProperty("settings", juce::var(constants));
                runObject->setProperty("metrics", juce::var(metrics));
                runs.add(juce::var(runObject));
            }
            result->setProperty("runs", runs);
        }
        return juce::var(result);
    }

    juce::String compactFamilyForSource(const juce::String& source)
    {
        if (source.containsIgnoreCase("rms_peak_sine")) return "levelDynamicsSine";
        if (source.containsIgnoreCase("rms_peak_noise")) return "levelDynamicsNoise";
        if (source.containsIgnoreCase("rms_peak_sweep")) return "levelDynamicsSweep";
        if (source.containsIgnoreCase("grid_thd")) return "harmonics";
        if (source.containsIgnoreCase("harmonic_fingerprint")) return "harmonicFrequencyProfile";
        if (source.containsIgnoreCase("transfer_curves_sine")) return "transferSine";
        if (source.containsIgnoreCase("transfer_curves_noise")) return "transferNoise";
        if (source.containsIgnoreCase("transfer_curves_sweep")) return "transferSweep";
        if (source.containsIgnoreCase("linear_response_sine")) return "frequencyResponseSine";
        if (source.containsIgnoreCase("linear_response_noise")) return "frequencyResponseNoise";
        if (source.containsIgnoreCase("linear_response_sweep")) return "frequencyResponseSweep";
        if (source.containsIgnoreCase("residual_dna")) return "residual";
        if (source.containsIgnoreCase("interaction_dna")) return "interaction";
        if (source.containsIgnoreCase("timing_dna")) return "timingPhase";
        if (source.containsIgnoreCase("alias_stress")) return "alias";
        if (source.containsIgnoreCase("stereo_dna")) return "stereo";
        if (source.containsIgnoreCase("summing_dna")) return "summing";
        if (source.containsIgnoreCase("boundary_dna")) return "boundaries";
        return source.upToLastOccurrenceOf(".csv", false, true);
    }

    juce::var buildOperatingPointDNA(const juce::Array<juce::var>& evidence)
    {
        std::map<int, juce::DynamicObject*> points;
        for (const auto& evidenceVar : evidence)
        {
            auto* evidenceObject = evidenceVar.getDynamicObject();
            if (evidenceObject == nullptr) continue;
            const auto source = evidenceObject->getProperty("source").toString();
            const auto family = compactFamilyForSource(source);
            auto* runs = evidenceObject->getProperty("runs").getArray();
            if (runs == nullptr) continue;

            for (const auto& runVar : *runs)
            {
                auto* run = runVar.getDynamicObject();
                if (run == nullptr) continue;
                const int runId = static_cast<int>(run->getProperty("runId"));
                auto it = points.find(runId);
                if (it == points.end())
                {
                    auto* point = new juce::DynamicObject();
                    point->setProperty("runId", runId);
                    point->setProperty("settings", run->getProperty("settings"));
                    point->setProperty("measurements", juce::var(new juce::DynamicObject()));
                    points.emplace(runId, point);
                    it = points.find(runId);
                }
                auto* measurements = it->second->getProperty("measurements").getDynamicObject();
                auto* metrics = run->getProperty("metrics").getDynamicObject();
                if (measurements == nullptr || metrics == nullptr) continue;
                const auto derived = metrics->getProperty("derived");
                if (derived.getDynamicObject() != nullptr && derived.getDynamicObject()->getProperties().size() > 0)
                    measurements->setProperty(family, derived);
            }
        }

        juce::Array<juce::var> array;
        for (const auto& item : points) array.add(juce::var(item.second));
        return juce::var(array);
    }

    juce::var buildCompactThresholds(const juce::var& operatingPointVar)
    {
        auto* result = new juce::DynamicObject();
        auto* extrema = new juce::DynamicObject();
        double worstAlias = -300.0, highestClean = 0.0, firstAlias = std::numeric_limits<double>::infinity();
        double maxCrestReduction = 0.0, maxResponseDeviation = 0.0, maxPhase = 0.0;
        double maxResidual = -300.0, maxInteraction = -300.0;
        bool haveAlias=false, haveResidual=false, haveInteraction=false;

        auto statMean = [](const juce::var& v) -> double
        {
            if (v.isDouble() || v.isInt() || v.isInt64()) return static_cast<double>(v);
            if (auto* o=v.getDynamicObject())
            {
                const auto m=o->getProperty("mean");
                if (m.isDouble() || m.isInt() || m.isInt64()) return static_cast<double>(m);
            }
            return std::numeric_limits<double>::quiet_NaN();
        };
        auto number = [](juce::DynamicObject* o, const char* key) -> double
        {
            if (o==nullptr) return std::numeric_limits<double>::quiet_NaN();
            const auto v=o->getProperty(key);
            return (v.isDouble()||v.isInt()||v.isInt64()) ? static_cast<double>(v) : std::numeric_limits<double>::quiet_NaN();
        };

        if (auto* points=operatingPointVar.getArray())
        {
            for (const auto& pointVar:*points)
            {
                auto* point=pointVar.getDynamicObject(); if(point==nullptr)continue;
                auto* measurements=point->getProperty("measurements").getDynamicObject(); if(measurements==nullptr)continue;
                if(auto* a=measurements->getProperty("alias").getDynamicObject())
                {
                    const double w=number(a,"worstAliasDb"); if(std::isfinite(w)){worstAlias=std::max(worstAlias,w);haveAlias=true;}
                    const double h=number(a,"highestCleanFundamentalHzAtMinus90Db"); if(std::isfinite(h))highestClean=std::max(highestClean,h);
                    const double f=number(a,"firstAliasFundamentalHzAboveMinus90Db"); if(std::isfinite(f))firstAlias=std::min(firstAlias,f);
                }
                if(auto* d=measurements->getProperty("levelDynamicsSine").getDynamicObject())
                {
                    const double c=statMean(d->getProperty("crestFactorChangeDb")); if(std::isfinite(c))maxCrestReduction=std::max(maxCrestReduction,-c);
                }
                for(const char* fam:{"frequencyResponseSine","frequencyResponseNoise","frequencyResponseSweep"})
                    if(auto* f=measurements->getProperty(fam).getDynamicObject())
                    {
                        const double r=number(f,"responseRippleDb"); if(std::isfinite(r))maxResponseDeviation=std::max(maxResponseDeviation,std::abs(r));
                    }
                if(auto* t=measurements->getProperty("timingPhase").getDynamicObject())
                { const double p=number(t,"maximumAbsolutePhaseRotationDegrees"); if(std::isfinite(p))maxPhase=std::max(maxPhase,p); }
                if(auto* r=measurements->getProperty("residual").getDynamicObject())
                { const double v=statMean(r->getProperty("overallResidualDb")); if(std::isfinite(v)){maxResidual=std::max(maxResidual,v);haveResidual=true;} }
                if(auto* i=measurements->getProperty("interaction").getDynamicObject())
                { const double v=number(i,"worstProductDb"); if(std::isfinite(v)){maxInteraction=std::max(maxInteraction,v);haveInteraction=true;} }
            }
        }
        if(haveAlias)extrema->setProperty("worstAliasDb",worstAlias);
        if(highestClean>0)extrema->setProperty("highestCleanFundamentalHzAtMinus90Db",highestClean);
        if(std::isfinite(firstAlias))extrema->setProperty("firstAliasFundamentalHzAboveMinus90Db",firstAlias);
        extrema->setProperty("maximumCrestReductionDb",maxCrestReduction);
        extrema->setProperty("maximumFrequencyResponseRippleDb",maxResponseDeviation);
        extrema->setProperty("maximumAbsolutePhaseRotationDegrees",maxPhase);
        if(haveResidual)extrema->setProperty("maximumMeanResidualDbRelativeToOutput",maxResidual);
        if(haveInteraction)extrema->setProperty("worstInteractionProductDbRelativeToInput",maxInteraction);
        result->setProperty("extrema",juce::var(extrema));
        result->setProperty("note","Thresholds are measurement-derived extrema. Operating-region labels are intentionally deferred until the policy is validated across a larger plugin corpus.");
        return juce::var(result);
    }

    juce::Array<juce::var> generatedBucketValues(const ParameterBucketConfig& bucket)
    {
        juce::Array<juce::var> values;
        BucketSpec spec;
        spec.paramName = bucket.paramName;
        spec.strategy = BucketSpec::strategyFromString(bucket.strategy);
        spec.min = bucket.min;
        spec.max = bucket.max;
        spec.numBuckets = bucket.numBuckets;
        spec.values = bucket.values;
        auto generated = spec.generateValues();
        if (bucket.includePluginDefault && !bucket.strategy.equalsIgnoreCase("Enumerated"))
        {
            const auto defaultValue = juce::jlimit(0.0f, 1.0f, bucket.pluginDefaultValue);
            const auto alreadyIncluded = std::any_of(generated.begin(), generated.end(), [defaultValue](float value)
            {
                return std::abs(value - defaultValue) < 1.0e-5f;
            });
            if (!alreadyIncluded)
                generated.push_back(defaultValue);
            std::sort(generated.begin(), generated.end());
            generated.erase(std::unique(generated.begin(), generated.end(), [](float a, float b)
            {
                return std::abs(a - b) < 1.0e-5f;
            }), generated.end());
        }
        for (const auto value : generated)
            values.add(value);
        return values;
    }

MainComponent::MainComponent()
    : pluginPathLabel("Plugin Path:", "Plugin Path:"), pluginInfoLabel("", "No plugin loaded"),
      serialPluginLabel("Plugin 2 (serial):", "Plugin 2 (serial):"), serialPluginInfoLabel("", "No serial plugin loaded"),
      parametersLabel("", ""), outputPathLabel("Output Path:", "Output Path:"),
      progressLabel("", "Ready"), progressBar(progress) {
    buildVersionLabel.setText({}, juce::dontSendNotification);
    buildVersionLabel.setFont(juce::Font(juce::FontOptions(18.0f, juce::Font::bold)));
    buildVersionLabel.setVisible(false);
    buildVersionLabel.setJustificationType(juce::Justification::centred);
    addAndMakeVisible(buildVersionLabel);

    // Plugin path section
    addAndMakeVisible(pluginPathLabel);
    // Default plugin path - set to dev plugin in debug builds
#ifdef DEBUG
    pluginPathEditor.setText("/Volumes/External SSD/Plug-Ins/VST3/Acustica/GAINSTATION2.vst3",
                             juce::dontSendNotification);
#else
    // Default plugin path - can be set to common VST3 location or left empty
    pluginPathEditor.setText("", juce::dontSendNotification);
#endif
    pluginPathEditor.addListener(this);
    addAndMakeVisible(pluginPathEditor);

    browseButton.setButtonText("Browse...");
    browseButton.addListener(this);
    addAndMakeVisible(browseButton);

    loadPluginButton.setButtonText("Load Plugin");
    loadPluginButton.addListener(this);
    addAndMakeVisible(loadPluginButton);

    openPluginButton.setButtonText("Open");
    openPluginButton.addListener(this);
    openPluginButton.setEnabled(false);
    addAndMakeVisible(openPluginButton);

    addAndMakeVisible(pluginInfoLabel);

    addAndMakeVisible(serialPluginLabel);
    addAndMakeVisible(serialPluginPathEditor);
    serialBrowseButton.setButtonText("Browse..."); serialBrowseButton.addListener(this); addAndMakeVisible(serialBrowseButton);
    serialLoadButton.setButtonText("Load"); serialLoadButton.addListener(this); addAndMakeVisible(serialLoadButton);
    serialOpenButton.setButtonText("Open"); serialOpenButton.addListener(this); serialOpenButton.setEnabled(false); addAndMakeVisible(serialOpenButton);
    addAndMakeVisible(serialPluginInfoLabel);

    // Expandable serial chain: Plugin 3..8 appear one at a time as the previous stage is loaded.
    for (int stageNumber = 3; stageNumber <= maxSerialStages; ++stageNumber)
    {
        auto stage = std::make_unique<AdditionalSerialStage>();
        stage->stageNumber = stageNumber;
        stage->label = std::make_unique<juce::Label>(juce::String(), "Plugin " + juce::String(stageNumber) + " (serial):");
        stage->pathEditor = std::make_unique<juce::TextEditor>();
        stage->browseButton = std::make_unique<juce::TextButton>("Browse...");
        stage->loadButton = std::make_unique<juce::TextButton>("Load");
        stage->openButton = std::make_unique<juce::TextButton>("Open");
        stage->removeButton = std::make_unique<juce::TextButton>("Remove");
        stage->parameterButton = std::make_unique<juce::TextButton>("Params");
        stage->targetButton = std::make_unique<juce::TextButton>("P" + juce::String(stageNumber));
        stage->infoLabel = std::make_unique<juce::Label>(juce::String(), "No plugin loaded");
        stage->browseButton->addListener(this);
        stage->loadButton->addListener(this);
        stage->openButton->addListener(this);
        stage->removeButton->addListener(this);
        stage->parameterButton->addListener(this);
        stage->targetButton->addListener(this);
        stage->parameterButton->setEnabled(false);
        stage->targetButton->setEnabled(false);
        stage->openButton->setEnabled(false);
        addAndMakeVisible(*stage->label);
        addAndMakeVisible(*stage->pathEditor);
        addAndMakeVisible(*stage->browseButton);
        addAndMakeVisible(*stage->loadButton);
        addAndMakeVisible(*stage->openButton);
        addAndMakeVisible(*stage->removeButton);
        addAndMakeVisible(*stage->parameterButton);
        addAndMakeVisible(*stage->targetButton);
        addAndMakeVisible(*stage->infoLabel);
        additionalSerialStages.push_back(std::move(stage));
    }
    refreshAdditionalSerialStageVisibility();

    // Parameter scan panel
    parameterScanGroup.setText("Parameter Scan");
    addAndMakeVisible(parameterScanGroup);

    parametersLabel.setText("Choose parameters to scan", juce::dontSendNotification);
    parametersLabel.setFont(juce::Font(juce::FontOptions(14.0f, juce::Font::bold)));
    addAndMakeVisible(parametersLabel);

    parameterTargetLabel.setText("Parameter target", juce::dontSendNotification);
    addAndMakeVisible(parameterTargetLabel);
    plugin1ParameterButton.setButtonText("P1");
    plugin1ParameterButton.addListener(this);
    plugin1ParameterButton.setColour(juce::TextButton::buttonColourId, juce::Colours::darkcyan);
    addAndMakeVisible(plugin1ParameterButton);
    plugin2ParameterButton.setButtonText("P2");
    plugin2ParameterButton.addListener(this);
    addAndMakeVisible(plugin2ParameterButton);

    parameterSearchEditor.setTextToShowWhenEmpty("Search parameters...", juce::Colours::grey);
    parameterSearchEditor.addListener(this);
    addAndMakeVisible(parameterSearchEditor);

    parameterListBox.setModel(this);
    parameterListBox.setRowSelectedOnMouseDown(false); // Disable row selection - we handle clicks ourselves
    parameterListBox.setMultipleSelectionEnabled(false);
    addAndMakeVisible(parameterListBox);

    selectAllButton.setButtonText("Select All");
    selectAllButton.addListener(this);
    addAndMakeVisible(selectAllButton);

    deselectAllButton.setButtonText("Deselect All");
    deselectAllButton.addListener(this);
    addAndMakeVisible(deselectAllButton);

    noParameterScanButton.setButtonText("No parameter scan");
    noParameterScanButton.addListener(this);
    noParameterScanButton.setTooltip("Run every enabled analyser through the plugin exactly as it is currently configured, without moving any plugin parameter.");
    addAndMakeVisible(noParameterScanButton);

    // Parameter config viewport
    parameterConfigViewport.setViewedComponent(&parameterConfigContainer, false);
    addAndMakeVisible(parameterConfigViewport);

    // Measurement config
    measurementConfig = std::make_unique<MeasurementConfigComponent>();
    addAndMakeVisible(measurementConfig.get());

    // Output path
    addAndMakeVisible(outputPathLabel);
    outputPathEditor.setText(juce::File::getSpecialLocation(juce::File::userHomeDirectory)
                                 .getChildFile("plugin_measure_output")
                                 .getFullPathName(),
                             juce::dontSendNotification);
    outputPathEditor.addListener(this);
    addAndMakeVisible(outputPathEditor);

    browseOutputButton.setButtonText("Browse...");
    browseOutputButton.addListener(this);
    addAndMakeVisible(browseOutputButton);

    exportDataButton.setButtonText("Export Raw Data");
    exportDataButton.setTooltip("Create a lossless full-evidence folder containing every analyser CSV plus the compact evidence JSON.");
    exportDataButton.addListener(this);
    exportDataButton.setEnabled(false);
    addAndMakeVisible(exportDataButton);

    exportEvidenceButton.setButtonText("Export Evidence");
    exportEvidenceButton.addListener(this);
    exportEvidenceButton.setEnabled(false);
    addAndMakeVisible(exportEvidenceButton);

    // Human-readable result summary
    resultsSummaryGroup.setText("Signal Summary");
    addAndMakeVisible(resultsSummaryGroup);

    copySummaryButton.setButtonText("Copy Summary");
    copySummaryButton.addListener(this);
    copySummaryButton.setEnabled(false);
    addAndMakeVisible(copySummaryButton);

    for (auto* label : { &identityRoleLabel, &identityWhyLabel, &identityWatchLabel, &levelSummaryLabel, &peakSummaryLabel, &dynamicsSummaryLabel, &harmonicSummaryLabel, &harmonicBalanceLabel, &harmonicGrowthLabel,
                         &toneBassLabel, &toneMidLabel, &toneTrebleLabel, &toneLargestLabel,
                         &behaviourSummaryLabel, &behaviourChangesLabel,
                         &operatingRangeLabel, &characterStartsLabel, &underStressLabel,
                         &nonlinearityLabel, &curveBehaviourLabel,
                         &waveformShapeLabel, &waveformKneeLabel, &waveformStyleLabel,
                         &temporalResponseLabel, &temporalAttackLabel, &temporalAfterLabel })
    {
        label->setFont(juce::Font(15.0f));
        label->setJustificationType(juce::Justification::centredLeft);
        addAndMakeVisible(label);
    }

    clearResultsSummary();

    // Run button
    runMeasurementButton.setButtonText("Run Measurement");
    runMeasurementButton.addListener(this);
    runMeasurementButton.setEnabled(false);
    addAndMakeVisible(runMeasurementButton);
    stopMeasurementButton.setButtonText("Stop Measurement");
    stopMeasurementButton.addListener(this);
    stopMeasurementButton.setEnabled(false);
    addAndMakeVisible(stopMeasurementButton);

    // Progress
    addAndMakeVisible(progressLabel);
    addAndMakeVisible(progressBar);
    progressBar.setPercentageDisplay(false);

    setSize(1240, 1900);
}

MainComponent::~MainComponent() {}

void MainComponent::paint(juce::Graphics& g) {
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));
}

void MainComponent::resized() {
    auto bounds = getLocalBounds().reduced(10);

    buildVersionLabel.setBounds(bounds.removeFromTop(32));
    bounds.removeFromTop(4);

    // Plugin path section (top)
    int visibleExtraStages = 0;
    for (const auto& stage : additionalSerialStages) if (stage->label->isVisible()) ++visibleExtraStages;
    auto pluginSection = bounds.removeFromTop(160 + visibleExtraStages * 74);
    pluginPathLabel.setBounds(pluginSection.removeFromTop(25));
    auto pluginRow = pluginSection.removeFromTop(30);
    pluginPathEditor.setBounds(pluginRow.removeFromLeft(600));
    pluginRow.removeFromLeft(10);
    browseButton.setBounds(pluginRow.removeFromLeft(100));
    pluginRow.removeFromLeft(10);
    loadPluginButton.setBounds(pluginRow.removeFromLeft(120));
    pluginRow.removeFromLeft(10);
    openPluginButton.setBounds(pluginRow.removeFromLeft(80));
    pluginInfoLabel.setBounds(pluginSection.removeFromTop(25));
    serialPluginLabel.setBounds(pluginSection.removeFromTop(25));
    auto serialRow = pluginSection.removeFromTop(30);
    serialPluginPathEditor.setBounds(serialRow.removeFromLeft(600)); serialRow.removeFromLeft(10);
    serialBrowseButton.setBounds(serialRow.removeFromLeft(100)); serialRow.removeFromLeft(10);
    serialLoadButton.setBounds(serialRow.removeFromLeft(120)); serialRow.removeFromLeft(10);
    serialOpenButton.setBounds(serialRow.removeFromLeft(80));
    serialPluginInfoLabel.setBounds(pluginSection.removeFromTop(25));

    for (auto& stage : additionalSerialStages)
    {
        if (!stage->label->isVisible())
            continue;
        stage->label->setBounds(pluginSection.removeFromTop(22));
        auto stageRow = pluginSection.removeFromTop(30);
        stage->pathEditor->setBounds(stageRow.removeFromLeft(520)); stageRow.removeFromLeft(8);
        stage->browseButton->setBounds(stageRow.removeFromLeft(90)); stageRow.removeFromLeft(8);
        stage->loadButton->setBounds(stageRow.removeFromLeft(80)); stageRow.removeFromLeft(8);
        stage->openButton->setBounds(stageRow.removeFromLeft(70)); stageRow.removeFromLeft(8);
        stage->parameterButton->setBounds(stageRow.removeFromLeft(76)); stageRow.removeFromLeft(8);
        stage->removeButton->setBounds(stageRow.removeFromLeft(80));
        stage->infoLabel->setBounds(pluginSection.removeFromTop(22));
    }
    bounds.removeFromTop(10);

    // Main laboratory panels (side by side, equal height)
    auto contentArea = bounds.removeFromTop(620);
    const int panelGap = 12;
    const int leftWidth = (contentArea.getWidth() - panelGap) / 2;

    auto leftPanelBounds = contentArea.removeFromLeft(leftWidth);
    contentArea.removeFromLeft(panelGap);
    auto rightPanelBounds = contentArea;

    parameterScanGroup.setBounds(leftPanelBounds);
    auto leftPanel = leftPanelBounds.reduced(12, 28);

    auto targetRow = leftPanel.removeFromTop(28);
    parameterTargetLabel.setBounds(targetRow.removeFromLeft(110));
    const int targetButtonWidth = 42;
    const int targetGap = 4;
    plugin1ParameterButton.setBounds(targetRow.removeFromLeft(targetButtonWidth));
    targetRow.removeFromLeft(targetGap);
    plugin2ParameterButton.setBounds(targetRow.removeFromLeft(targetButtonWidth));
    for (auto& stage : additionalSerialStages)
    {
        targetRow.removeFromLeft(targetGap);
        stage->targetButton->setBounds(targetRow.removeFromLeft(targetButtonWidth));
    }
    leftPanel.removeFromTop(5);
    noParameterScanButton.setBounds(leftPanel.removeFromTop(30));
    leftPanel.removeFromTop(6);
    parametersLabel.setBounds(leftPanel.removeFromTop(24));
    parameterSearchEditor.setBounds(leftPanel.removeFromTop(30));
    leftPanel.removeFromTop(6);

    auto buttonRow = leftPanel.removeFromTop(28);
    selectAllButton.setBounds(buttonRow.removeFromLeft(90));
    buttonRow.removeFromLeft(8);
    deselectAllButton.setBounds(buttonRow.removeFromLeft(90));
    leftPanel.removeFromTop(6);

    parameterListBox.setBounds(leftPanel.removeFromTop(170));
    leftPanel.removeFromTop(8);
    parameterConfigViewport.setBounds(leftPanel);

    measurementConfig->setBounds(rightPanelBounds);

    bounds.removeFromTop(10);

    // Bottom: Output and run
    auto outputSection = bounds.removeFromTop(50);
    outputPathLabel.setBounds(outputSection.removeFromTop(25));
    auto outputRow = outputSection.removeFromTop(30);
    outputPathEditor.setBounds(outputRow.removeFromLeft(600));
    outputRow.removeFromLeft(10);
    browseOutputButton.setBounds(outputRow.removeFromLeft(100));
    outputRow.removeFromLeft(10);
    runMeasurementButton.setBounds(outputRow.removeFromLeft(150));
    outputRow.removeFromLeft(8);
    stopMeasurementButton.setBounds(outputRow.removeFromLeft(150));
    outputRow.removeFromLeft(10);
    exportDataButton.setBounds(outputRow.removeFromLeft(135));
    outputRow.removeFromLeft(10);
    exportEvidenceButton.setBounds(outputRow.removeFromLeft(145));

    bounds.removeFromTop(10);

    // Progress
    progressLabel.setBounds(bounds.removeFromTop(25));
    progressBar.setBounds(bounds.removeFromTop(24));
    bounds.removeFromTop(8);

    constexpr int summaryLabelHeight = 25;
    constexpr int summaryGap = 8;
    constexpr int summaryGroupInsets = 78;
    constexpr int summaryLabelCount = 26;
    constexpr int summaryGapCount = 7;
    constexpr int summaryHeight = summaryGroupInsets
                                + summaryLabelCount * summaryLabelHeight
                                + summaryGapCount * summaryGap;

    auto summaryBounds = bounds.removeFromTop(summaryHeight);
    resultsSummaryGroup.setBounds(summaryBounds);
    auto summaryHeader = summaryBounds.removeFromTop(44).reduced(12, 8);
    copySummaryButton.setBounds(summaryHeader.removeFromRight(130));
    summaryBounds = summaryBounds.reduced(12, 0);
    identityRoleLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    identityWhyLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    identityWatchLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    summaryBounds.removeFromTop(summaryGap);
    levelSummaryLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    peakSummaryLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    dynamicsSummaryLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    harmonicSummaryLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    harmonicBalanceLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    harmonicGrowthLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    summaryBounds.removeFromTop(summaryGap);
    toneBassLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    toneMidLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    toneTrebleLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    toneLargestLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    summaryBounds.removeFromTop(summaryGap);
    behaviourSummaryLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    behaviourChangesLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    summaryBounds.removeFromTop(summaryGap);
    operatingRangeLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    characterStartsLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    underStressLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    summaryBounds.removeFromTop(summaryGap);
    nonlinearityLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    curveBehaviourLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    summaryBounds.removeFromTop(summaryGap);
    waveformShapeLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    waveformKneeLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    waveformStyleLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    summaryBounds.removeFromTop(summaryGap);
    temporalResponseLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    temporalAttackLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
    temporalAfterLabel.setBounds(summaryBounds.removeFromTop(summaryLabelHeight));
}

void MainComponent::buttonClicked(juce::Button* button) {
    for (int i = 0; i < (int)additionalSerialStages.size(); ++i)
    {
        auto& stage = *additionalSerialStages[(size_t)i];
        if (button == stage.browseButton.get())
        {
            auto chooser = std::make_shared<juce::FileChooser>("Select Plugin " + juce::String(stage.stageNumber) + " (VST3 or Audio Unit)", juce::File(), "*.vst3;*.component");
            auto chooserFlags = juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles;
            chooser->launchAsync(chooserFlags, [this, chooser, i](const juce::FileChooser& fc) {
                if (fc.getResults().size() > 0 && i < (int)additionalSerialStages.size())
                    additionalSerialStages[(size_t)i]->pathEditor->setText(fc.getResult().getFullPathName(), juce::dontSendNotification);
            });
            return;
        }
        if (button == stage.loadButton.get())
        {
            loadAdditionalSerialStage(i);
            return;
        }
        if (button == stage.openButton.get())
        {
            openPluginEditor(stage.pluginInstance.get(), stage.editorWindow,
                             stage.pluginInstance != nullptr ? stage.pluginInstance->getName() : "Serial Plugin");
            return;
        }
        if (button == stage.parameterButton.get() || button == stage.targetButton.get())
        {
            if (stage.pluginInstance == nullptr) { showError("Load Plugin " + juce::String(stage.stageNumber) + " first."); return; }
            saveCurrentParameterPanelState();
            editingParameterStage = stage.stageNumber;
            editingSerialParameters = false;
            availableParameters = stage.availableParameters;
            selectedParameters = stage.selectedParameters;
            parameterMap = stage.parameterMap;
            plugin1ParameterButton.setColour(juce::TextButton::buttonColourId, juce::Colours::darkgrey);
            plugin2ParameterButton.setColour(juce::TextButton::buttonColourId, juce::Colours::darkgrey);
            for (auto& targetStage : additionalSerialStages)
                targetStage->targetButton->setColour(juce::TextButton::buttonColourId,
                    targetStage.get() == &stage ? juce::Colours::darkcyan : juce::Colours::darkgrey);
            parametersLabel.setText("Plugin " + juce::String(stage.stageNumber) + " parameters", juce::dontSendNotification);
            filteredParameterIndices.clear(); for (int k=0;k<(int)availableParameters.size();++k) filteredParameterIndices.push_back(k);
            parameterSearchEditor.clear(); parameterListBox.updateContent(); parameterListBox.repaint();
            updateParameterList();
            for (auto& comp : parameterConfigComponents) for (const auto& bucket : stage.savedBuckets) if (comp->getConfig().paramName == bucket.paramName) comp->setConfig(bucket);
            return;
        }
        if (button == stage.removeButton.get())
        {
            stage.editorWindow.reset();
            stage.pluginInstance.reset();
            stage.pathEditor->clear();
            stage.infoLabel->setText("No plugin loaded", juce::dontSendNotification);
            stage.openButton->setEnabled(false);
            stage.parameterButton->setEnabled(false);
            stage.targetButton->setEnabled(false);
            stage.availableParameters.clear(); stage.selectedParameters.clear(); stage.parameterMap.clear(); stage.savedBuckets.clear();
            // Removing a stage also clears every later stage so the serial chain remains contiguous.
            for (int j = i + 1; j < (int)additionalSerialStages.size(); ++j)
            {
                auto& later = *additionalSerialStages[(size_t)j];
                later.editorWindow.reset(); later.pluginInstance.reset(); later.pathEditor->clear();
                later.infoLabel->setText("No plugin loaded", juce::dontSendNotification);
                later.openButton->setEnabled(false); later.parameterButton->setEnabled(false); later.targetButton->setEnabled(false);
                later.availableParameters.clear(); later.selectedParameters.clear(); later.parameterMap.clear(); later.savedBuckets.clear();
            }
            refreshAdditionalSerialStageVisibility();
            resized();
            return;
        }
    }
    if (button == &plugin1ParameterButton) {
        saveCurrentParameterPanelState();
        editingParameterStage = 1;
        loadParameterPanelState(false);
        return;
    }
    if (button == &plugin2ParameterButton) {
        if (serialPluginInstance == nullptr) {
            showError("Load Plugin 2 before selecting its parameters.");
            return;
        }
        saveCurrentParameterPanelState();
        editingParameterStage = 2;
        loadParameterPanelState(true);
        return;
    }
    if (button == &browseButton) {
        auto chooser = std::make_shared<juce::FileChooser>("Select Plugin (VST3 or Audio Unit)", juce::File(), "*.vst3;*.component");
        auto chooserFlags = juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles;
        chooser->launchAsync(chooserFlags, [this, chooser](const juce::FileChooser& fc) {
            if (fc.getResults().size() > 0) {
                pluginPathEditor.setText(fc.getResult().getFullPathName(), juce::dontSendNotification);
            }
        });
    } else if (button == &serialBrowseButton) {
        auto chooser = std::make_shared<juce::FileChooser>("Select second Plugin (VST3 or Audio Unit)", juce::File(), "*.vst3;*.component");
        auto chooserFlags = juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles;
        chooser->launchAsync(chooserFlags, [this, chooser](const juce::FileChooser& fc) {
            if (fc.getResults().size() > 0) serialPluginPathEditor.setText(fc.getResult().getFullPathName(), juce::dontSendNotification);
        });
    } else if (button == &browseOutputButton) {
        auto chooser =
            std::make_shared<juce::FileChooser>("Select Output Directory", juce::File(outputPathEditor.getText()));
        auto chooserFlags = juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectDirectories;
        chooser->launchAsync(chooserFlags, [this, chooser](const juce::FileChooser& fc) {
            if (fc.getResults().size() > 0) {
                outputPathEditor.setText(fc.getResult().getFullPathName(), juce::dontSendNotification);
            }
        });
    } else if (button == &loadPluginButton) {
        loadPlugin();
    } else if (button == &openPluginButton) {
        openPluginEditor(pluginInstance.get(), pluginEditorWindow, pluginInstance != nullptr ? pluginInstance->getName() : "Plugin");
    } else if (button == &serialLoadButton) {
        loadSerialPlugin();
    } else if (button == &serialOpenButton) {
        openPluginEditor(serialPluginInstance.get(), serialPluginEditorWindow, serialPluginInstance != nullptr ? serialPluginInstance->getName() : "Serial Plugin");
    } else if (button == &selectAllButton) {
        selectedParameters.resize(availableParameters.size(), true);
        std::fill(selectedParameters.begin(), selectedParameters.end(), true);
        parameterListBox.updateContent();
        // Repaint all visible rows
        int firstVisible = parameterListBox.getRowContainingPosition(0, 0);
        int lastVisible = parameterListBox.getRowContainingPosition(0, parameterListBox.getHeight());
        for (int i = firstVisible; i <= lastVisible && i >= 0 && i < (int)availableParameters.size(); ++i) {
            parameterListBox.repaintRow(i);
        }
        parameterListBox.repaint(); // Also repaint the whole component
        updateParameterList();
    } else if (button == &deselectAllButton) {
        selectedParameters.resize(availableParameters.size(), false);
        std::fill(selectedParameters.begin(), selectedParameters.end(), false);
        parameterListBox.updateContent();
        // Repaint all visible rows
        int firstVisible = parameterListBox.getRowContainingPosition(0, 0);
        int lastVisible = parameterListBox.getRowContainingPosition(0, parameterListBox.getHeight());
        for (int i = firstVisible; i <= lastVisible && i >= 0 && i < (int)availableParameters.size(); ++i) {
            parameterListBox.repaintRow(i);
        }
        parameterListBox.repaint(); // Also repaint the whole component
        updateParameterList();
    } else if (button == &noParameterScanButton) {
        const bool baseline = noParameterScanButton.getToggleState();
        if (baseline)
        {
            std::fill(selectedParameters.begin(), selectedParameters.end(), false);
            parameterListBox.updateContent();
            parameterListBox.repaint();
            updateParameterList();
        }
        parameterListBox.setEnabled(!baseline);
        selectAllButton.setEnabled(!baseline);
        deselectAllButton.setEnabled(!baseline);
        parameterConfigViewport.setEnabled(!baseline);
        progressLabel.setText(baseline ? "Baseline mode: current plugin settings will be tested" : "Parameter scan enabled",
                              juce::dontSendNotification);
    } else if (button == &runMeasurementButton) {
        runMeasurement();
    } else if (button == &stopMeasurementButton) {
        measurementCancelRequested.store(true);
        stopMeasurementButton.setEnabled(false);
        progressLabel.setText("Stopping measurement...", juce::dontSendNotification);
    } else if (button == &copySummaryButton) {
        juce::SystemClipboard::copyTextToClipboard(buildSummaryText());
        progressLabel.setText("Summary copied to clipboard", juce::dontSendNotification);
    } else if (button == &exportDataButton) {
        exportMeasurementData();
    } else if (button == &exportEvidenceButton) {
        exportEvidencePack();
    }
}

void MainComponent::textEditorTextChanged(juce::TextEditor& editor) {
    if (&editor == &parameterSearchEditor)
    {
        filteredParameterIndices.clear();
        const auto query = parameterSearchEditor.getText().trim().toLowerCase();
        for (int i = 0; i < static_cast<int>(availableParameters.size()); ++i)
            if (query.isEmpty() || availableParameters[static_cast<std::size_t>(i)].toLowerCase().contains(query))
                filteredParameterIndices.push_back(i);
        parameterListBox.updateContent();
        parameterListBox.repaint();
    }
}

int MainComponent::getNumRows() {
    return (int)filteredParameterIndices.size();
}

void MainComponent::paintListBoxItem(int rowNumber, juce::Graphics& g, int width, int height, bool rowIsSelected) {
    if (rowNumber < 0 || rowNumber >= (int)filteredParameterIndices.size())
        return;
    const int parameterIndex = filteredParameterIndices[static_cast<std::size_t>(rowNumber)];

    // Background - use subtle alternating colors, ignore rowIsSelected since we disabled selection
    g.fillAll(rowNumber % 2 == 0 ? juce::Colours::white : juce::Colours::lightgrey.withAlpha(0.3f));

    // Draw checkbox
    const int checkboxSize = 18;
    const int checkboxX = 5;
    const int checkboxY = (height - checkboxSize) / 2;
    juce::Rectangle<float> checkboxBounds((float)checkboxX, (float)checkboxY, (float)checkboxSize, (float)checkboxSize);

    bool isChecked = parameterIndex < (int)selectedParameters.size() && selectedParameters[static_cast<std::size_t>(parameterIndex)];

    // Draw checkbox border
    g.setColour(juce::Colours::darkgrey);
    g.drawRect(checkboxBounds, 1.5f);

    // Draw checkbox fill if checked
    if (isChecked) {
        g.setColour(juce::Colours::blue);
        g.fillRect(checkboxBounds.reduced(2.0f));

        // Draw checkmark using a simple path
        g.setColour(juce::Colours::white);
        juce::Path checkmark;
        const float x = checkboxBounds.getX();
        const float y = checkboxBounds.getY();
        const float w = checkboxBounds.getWidth();
        const float h = checkboxBounds.getHeight();

        // Draw checkmark as two connected lines
        checkmark.startNewSubPath(x + w * 0.2f, y + h * 0.5f);
        checkmark.lineTo(x + w * 0.45f, y + h * 0.75f);
        checkmark.lineTo(x + w * 0.8f, y + h * 0.25f);

        g.strokePath(checkmark, juce::PathStrokeType(2.5f, juce::PathStrokeType::curved));
    }

    // Parameter name
    g.setColour(juce::Colours::black);
    g.setFont(14.0f);
    g.drawText(availableParameters[static_cast<std::size_t>(parameterIndex)], checkboxX + checkboxSize + 10, 0, width - checkboxX - checkboxSize - 10,
               height, juce::Justification::centredLeft);
}

void MainComponent::listBoxItemClicked(int row, const juce::MouseEvent&) {
    if (selectedParameters.size() != availableParameters.size())
        selectedParameters.resize(availableParameters.size(), false);

    if (row < 0 || row >= static_cast<int>(filteredParameterIndices.size()))
        return;

    const int parameterIndex = filteredParameterIndices[static_cast<std::size_t>(row)];
    if (parameterIndex < 0 || parameterIndex >= static_cast<int>(selectedParameters.size()))
        return;

    selectedParameters[static_cast<std::size_t>(parameterIndex)] = !selectedParameters[static_cast<std::size_t>(parameterIndex)];
    parameterListBox.repaintRow(row);
    updateParameterList();
}

void MainComponent::loadPlugin() {
    juce::String pluginPath = pluginPathEditor.getText();
    if (pluginPath.isEmpty()) {
        showError("Please specify a plugin path");
        return;
    }

    juce::File pluginFile(pluginPath);
    // VST3 plugins on macOS are bundles (directories), not files
    if (!pluginFile.exists()) {
        showError("Plugin file does not exist: " + pluginPath);
        return;
    }

    pluginEditorWindow.reset();
    progressLabel.setText("Loading plugin...", juce::dontSendNotification);

    // Get sample rate and block size from measurement config
    Config tempConfig;
    measurementConfig->fillConfig(tempConfig);

    juce::String errorMessage;
    pluginInstance = loadPluginInstance(pluginFile, tempConfig.sampleRate, tempConfig.blockSize, errorMessage);

    if (pluginInstance == nullptr) {
        showError(errorMessage.isEmpty() ? "Failed to load plugin" : errorMessage);
        return;
    }

    pluginInfoLabel.setText("Loaded: " + pluginInstance->getName(), juce::dontSendNotification);
    openPluginButton.setEnabled(pluginInstance->hasEditor());
    scanPluginParameters();
    runMeasurementButton.setEnabled(true);
    progressLabel.setText("Plugin loaded successfully", juce::dontSendNotification);
}


void MainComponent::openPluginEditor(juce::AudioPluginInstance* plugin, std::unique_ptr<juce::DocumentWindow>& window, const juce::String& title)
{
    if (plugin == nullptr || !plugin->hasEditor()) return;
    if (window != nullptr) { window->setVisible(true); window->toFront(true); return; }
    auto* editor = plugin->createEditorIfNeeded();
    if (editor == nullptr) return;
    class EditorWindow final : public juce::DocumentWindow
    {
    public:
        EditorWindow(const juce::String& n, juce::Component* c) : juce::DocumentWindow(n, juce::Colours::darkgrey, juce::DocumentWindow::closeButton)
        { setUsingNativeTitleBar(true); setContentNonOwned(c, true); centreWithSize(c->getWidth(), c->getHeight()); setResizable(true, true); }
        void closeButtonPressed() override { setVisible(false); }
    };
    window = std::make_unique<EditorWindow>(title, editor);
    window->setVisible(true);
}

void MainComponent::loadSerialPlugin()
{
    const auto path = serialPluginPathEditor.getText().trim();
    if (path.isEmpty()) {
        serialPluginEditorWindow.reset(); serialPluginInstance.reset(); serialOpenButton.setEnabled(false);
        serialPluginInfoLabel.setText("Plugin 2 disabled", juce::dontSendNotification);
        for (auto& stage : additionalSerialStages) { stage->editorWindow.reset(); stage->pluginInstance.reset(); stage->pathEditor->clear(); stage->openButton->setEnabled(false); }
        refreshAdditionalSerialStageVisibility(); resized(); return;
    }
    serialPluginEditorWindow.reset();
    Config tempConfig; measurementConfig->fillConfig(tempConfig);
    juce::String error;
    serialPluginInstance = loadPluginInstance(juce::File(path), tempConfig.sampleRate, tempConfig.blockSize, error);
    if (serialPluginInstance == nullptr) { showError(error); return; }
    serialPluginInfoLabel.setText("Loaded: " + serialPluginInstance->getName(), juce::dontSendNotification);
    serialOpenButton.setEnabled(serialPluginInstance->hasEditor());
    scanSerialPluginParameters();
    refreshAdditionalSerialStageVisibility();
    resized();
}

void MainComponent::refreshAdditionalSerialStageVisibility()
{
    bool revealNext = (serialPluginInstance != nullptr);
    for (auto& stage : additionalSerialStages)
    {
        const bool visible = revealNext;
        stage->label->setVisible(visible); stage->pathEditor->setVisible(visible);
        stage->browseButton->setVisible(visible); stage->loadButton->setVisible(visible);
        stage->openButton->setVisible(visible); stage->parameterButton->setVisible(visible); stage->removeButton->setVisible(visible);
        stage->targetButton->setVisible(stage->pluginInstance != nullptr);
        stage->infoLabel->setVisible(visible);
        revealNext = visible && stage->pluginInstance != nullptr;
    }
}

void MainComponent::loadAdditionalSerialStage(int index)
{
    if (index < 0 || index >= (int)additionalSerialStages.size()) return;
    auto& stage = *additionalSerialStages[(size_t)index];
    const auto path = stage.pathEditor->getText().trim();
    if (path.isEmpty())
    {
        showError("Choose a plugin for Plugin " + juce::String(stage.stageNumber) + " first.");
        return;
    }
    stage.editorWindow.reset();
    Config tempConfig; measurementConfig->fillConfig(tempConfig);
    juce::String error;
    stage.pluginInstance = loadPluginInstance(juce::File(path), tempConfig.sampleRate, tempConfig.blockSize, error);
    if (stage.pluginInstance == nullptr) { showError(error); return; }
    stage.infoLabel->setText("Loaded: " + stage.pluginInstance->getName(), juce::dontSendNotification);
    stage.openButton->setEnabled(stage.pluginInstance->hasEditor());
    stage.parameterButton->setEnabled(true);
    stage.targetButton->setEnabled(true);
    stage.parameterMap = buildParameterMap(*stage.pluginInstance, true);
    stage.availableParameters.clear();
    stage.selectedParameters.clear();
    for (const auto& [name, param] : stage.parameterMap) { stage.availableParameters.push_back(name); stage.selectedParameters.push_back(false); }
    std::sort(stage.availableParameters.begin(), stage.availableParameters.end());
    refreshAdditionalSerialStageVisibility();
    resized();
}

void MainComponent::scanPluginParameters() {
    if (pluginInstance == nullptr)
        return;

    parameterMap = buildParameterMap(*pluginInstance, true); // Only UI-exposed parameters
    availableParameters.clear();
    selectedParameters.clear();

    for (const auto& [name, param] : parameterMap) {
        availableParameters.push_back(name);
        selectedParameters.push_back(false);
    }

    std::sort(availableParameters.begin(), availableParameters.end());

    // Rebuild the visible row map after every plugin scan. Without this,
    // the ListBox model reports zero rows until the search field changes.
    filteredParameterIndices.clear();
    filteredParameterIndices.reserve(availableParameters.size());
    for (int i = 0; i < static_cast<int>(availableParameters.size()); ++i)
        filteredParameterIndices.push_back(i);

    primaryAvailableParameters = availableParameters;
    primarySelectedParameters = selectedParameters;
    primaryParameterMap = parameterMap;
    if (!editingSerialParameters) {
        parameterListBox.updateContent();
        parameterListBox.repaint();
        updateParameterList();
    }
}

void MainComponent::scanSerialPluginParameters() {
    if (serialPluginInstance == nullptr)
        return;
    serialParameterMap = buildParameterMap(*serialPluginInstance, true);
    serialAvailableParameters.clear();
    serialSelectedParameters.clear();
    for (const auto& [name, param] : serialParameterMap) {
        serialAvailableParameters.push_back(name);
        serialSelectedParameters.push_back(false);
    }
    std::sort(serialAvailableParameters.begin(), serialAvailableParameters.end());
    if (editingSerialParameters)
        loadParameterPanelState(true);
}

void MainComponent::saveCurrentParameterPanelState() {
    std::vector<ParameterBucketConfig> buckets;
    for (const auto& comp : parameterConfigComponents)
        buckets.push_back(comp->getConfig());
    if (editingParameterStage >= 3) {
        const int index = editingParameterStage - 3;
        if (index >= 0 && index < (int)additionalSerialStages.size()) {
            auto& stage = *additionalSerialStages[(size_t)index];
            stage.availableParameters = availableParameters; stage.selectedParameters = selectedParameters;
            stage.parameterMap = parameterMap; stage.savedBuckets = std::move(buckets);
        }
    } else if (editingParameterStage == 2) {
        serialAvailableParameters = availableParameters; serialSelectedParameters = selectedParameters;
        serialParameterMap = parameterMap; serialSavedBuckets = std::move(buckets);
    } else {
        primaryAvailableParameters = availableParameters; primarySelectedParameters = selectedParameters;
        primaryParameterMap = parameterMap; primarySavedBuckets = std::move(buckets);
    }
}

void MainComponent::loadParameterPanelState(bool serial) {
    editingSerialParameters = serial;
    editingParameterStage = serial ? 2 : 1;
    availableParameters = serial ? serialAvailableParameters : primaryAvailableParameters;
    selectedParameters = serial ? serialSelectedParameters : primarySelectedParameters;
    parameterMap = serial ? serialParameterMap : primaryParameterMap;
    plugin1ParameterButton.setColour(juce::TextButton::buttonColourId, serial ? juce::Colours::darkgrey : juce::Colours::darkcyan);
    plugin2ParameterButton.setColour(juce::TextButton::buttonColourId, serial ? juce::Colours::darkcyan : juce::Colours::darkgrey);
    for (auto& stage : additionalSerialStages)
        stage->targetButton->setColour(juce::TextButton::buttonColourId, juce::Colours::darkgrey);
    parametersLabel.setText(serial ? "Plugin 2 parameters" : "Plugin 1 parameters", juce::dontSendNotification);
    filteredParameterIndices.clear();
    for (int i = 0; i < static_cast<int>(availableParameters.size()); ++i)
        filteredParameterIndices.push_back(i);
    parameterSearchEditor.clear();
    parameterListBox.updateContent();
    parameterListBox.repaint();

    const auto saved = serial ? serialSavedBuckets : primarySavedBuckets;
    updateParameterList();
    for (auto& comp : parameterConfigComponents)
        for (const auto& bucket : saved)
            if (comp->getConfig().paramName == bucket.paramName)
                comp->setConfig(bucket);
}

void MainComponent::updateParameterList() {
    // Remove old config components
    parameterConfigComponents.clear();
    parameterConfigContainer.removeAllChildren();

    // Create config components for selected parameters
    int y = 10;
    for (size_t i = 0; i < availableParameters.size(); ++i) {
        if (i < selectedParameters.size() && selectedParameters[i]) {
            auto parameterIt = parameterMap.find(availableParameters[i]);
            auto* comp = new ParameterConfigComponent(availableParameters[i],
                parameterIt != parameterMap.end() ? parameterIt->second : nullptr);
            const int componentHeight = comp->getConfig().strategy == "Enumerated" ? 300 : 180;
            comp->setBounds(10, y, 480, componentHeight);
            parameterConfigContainer.addAndMakeVisible(comp);
            parameterConfigComponents.push_back(std::unique_ptr<ParameterConfigComponent>(comp));
            y += componentHeight + 10;
        }
    }

    parameterConfigContainer.setSize(500, y);
    parameterConfigViewport.setViewPosition(0, 0);
}

void MainComponent::runMeasurement() {
    if (pluginInstance == nullptr) {
        showError("No plugin loaded");
        return;
    }

    saveCurrentParameterPanelState();
    // Count selected parameters across both serial stages
    int selectedCount = 0;
    for (bool selected : primarySelectedParameters) if (selected) selectedCount++;
    for (bool selected : serialSelectedParameters) if (selected) selectedCount++;
    for (const auto& stage : additionalSerialStages) for (bool selected : stage->selectedParameters) if (selected) selectedCount++;

    if (selectedCount == 0 && !noParameterScanButton.getToggleState()) {
        showError("Select a parameter, or enable 'No parameter scan' to test the plugin at its current settings.");
        return;
    }

    // The output path is only used when Export Data is pressed. Measurement
    // itself stays entirely in memory and does not create files or folders.
    juce::File outDir(outputPathEditor.getText());

    // Build config from UI
    Config config = buildConfigFromUI();
    config.pluginPath = pluginPathEditor.getText();

    // Capture the exact loaded editor state (important for shell plugins such as
    // Airwindows Consolidated) and restore it into every fresh measurement instance.
    auto primaryState = std::make_shared<juce::MemoryBlock>();
    pluginInstance->getStateInformation(*primaryState);
    auto serialState = std::make_shared<juce::MemoryBlock>();
    const auto serialPath = serialPluginPathEditor.getText().trim();
    if (serialPluginInstance != nullptr)
        serialPluginInstance->getStateInformation(*serialState);

    std::vector<juce::String> additionalSerialPaths;
    std::vector<std::shared_ptr<juce::MemoryBlock>> additionalSerialStates;
    for (const auto& stage : additionalSerialStages)
    {
        if (stage->pluginInstance == nullptr) break;
        additionalSerialPaths.push_back(stage->pathEditor->getText().trim());
        auto state = std::make_shared<juce::MemoryBlock>();
        stage->pluginInstance->getStateInformation(*state);
        additionalSerialStates.push_back(std::move(state));
    }

    // Run measurement in background thread
    measurementCancelRequested.store(false);
    runMeasurementButton.setEnabled(false);
    stopMeasurementButton.setEnabled(true);
    exportDataButton.setEnabled(false);
    exportEvidenceButton.setEnabled(false);
    progressLabel.setText("Running measurement...", juce::dontSendNotification);
    clearResultsSummary();
    {
        std::lock_guard<std::mutex> lock(pendingExportMutex);
        pendingExportDatasets.clear();
        pendingEvidenceSummary.clear();
        hasPendingEvidencePack = false;
    }

    std::thread([this, config, outDir, primaryState, serialState, serialPath, additionalSerialPaths, additionalSerialStates]() {
        try {

            // Build parameter name list
            std::vector<juce::String> paramNames;
            for (const auto& bucket : config.parameterBuckets)
                paramNames.push_back(bucket.paramName);

            // Build run grid
            juce::MessageManager::callAsync(
                [this]() { progressLabel.setText("Building run grid...", juce::dontSendNotification); });
            auto runs = buildRunGrid(config, paramNames);

            // Warn if too many runs and estimate time
            if (runs.size() > 100000) {
                std::cerr << "[Measurement] WARNING: " << runs.size()
                          << " runs is very large. This may take a long time." << std::endl;

                // Estimate time: assume ~0.1 seconds per run (very rough estimate)
                double estimatedSeconds = runs.size() * 0.1;
                double estimatedMinutes = estimatedSeconds / 60.0;
                double estimatedHours = estimatedMinutes / 60.0;

                juce::String timeEstimate;
                if (estimatedHours >= 1.0) {
                    timeEstimate = juce::String(estimatedHours, 1) + " hours";
                } else if (estimatedMinutes >= 1.0) {
                    timeEstimate = juce::String(estimatedMinutes, 1) + " minutes";
                } else {
                    timeEstimate = juce::String(estimatedSeconds, 1) + " seconds";
                }

                juce::MessageManager::callAsync([this, runs, timeEstimate]() {
                    progressLabel.setText("WARNING: " + juce::String(runs.size()) + " runs (~" + timeEstimate +
                                              ") - Consider reducing parameters/buckets!",
                                          juce::dontSendNotification);
                });
                std::this_thread::sleep_for(std::chrono::seconds(3)); // Give user time to see warning
            }

            WaveformSummary waveformSummary {
                "WAVE SHAPE  Run Raw with a sine signal",
                "WHAT IT DOES  unavailable",
                "SIMILAR TO  unavailable"
            };
            TemporalSummary temporalSummary {
                "TEMPORAL RESPONSE  Run Raw with a sine signal",
                "ATTACK HANDLING  unavailable",
                "AFTER-EFFECTS  unavailable"
            };

            MeasurementDataset sineRmsResult;
            MeasurementDataset sineThdResult;
            MeasurementDataset sineTransferResult;
            MeasurementDataset sweepLinearResult;
            MeasurementDataset sweepHarmonicFingerprintResult;
            MeasurementDataset noiseLinearResult;
            MeasurementDataset interactionResult;
            MeasurementDataset timingResult;
            MeasurementDataset residualResult;
            MeasurementDataset boundaryResult;
            MeasurementDataset stereoResult;
            MeasurementDataset summingResult;
            MeasurementDataset aliasResult;
            std::vector<PendingExportDataset> completedExports;

            bool haveSineRaw = false;
            bool haveSineRms = false;
            bool haveSineThd = false;
            bool haveSineTransfer = false;
            bool haveSweepLinear = false;
            bool haveSweepHarmonicFingerprint = false;
            bool haveNoiseLinear = false;
            bool haveInteraction = false;
            bool haveTiming = false;
            bool haveResidual = false;
            bool haveBoundary = false;
            bool haveStereo = false;
            bool haveSumming = false;
            bool haveAlias = false;

            const auto analyserSelected = [&config](const juce::String& name)
            {
                return std::find_if(config.analyzers.begin(), config.analyzers.end(), [&](const juce::String& item)
                {
                    return item.equalsIgnoreCase(name);
                }) != config.analyzers.end();
            };

            std::vector<juce::String> passSignals;
            if (config.signalType.equalsIgnoreCase("all"))
            {
                const bool needsSine = analyserSelected("RawCsv") || analyserSelected("RmsPeak")
                                    || analyserSelected("TransferCurve") || analyserSelected("Thd");
                const bool needsLinear = analyserSelected("LinearResponse");
                if (needsSine) passSignals.push_back("sine");
                if (needsLinear) { passSignals.push_back("sweep"); passSignals.push_back("noise"); }
                if (analyserSelected("Interaction")) passSignals.push_back("interaction");
                if (analyserSelected("Timing")) passSignals.push_back("timing");
                if (analyserSelected("Residual")) passSignals.push_back("residual");
                if (analyserSelected("Boundary")) passSignals.push_back("boundary");
                if (analyserSelected("Stereo")) passSignals.push_back("stereo");
                if (analyserSelected("Summing")) passSignals.push_back("summing");
                if (analyserSelected("Alias")) passSignals.push_back("alias");
            }
            else
                passSignals = { config.signalType };

            if (passSignals.empty())
                throw std::runtime_error("No analysers are enabled for the selected measurement mode.");

            const int totalPassRuns = (int)(runs.size() * passSignals.size());
            int completedPassRuns = 0;
            auto overallStartTime = std::chrono::steady_clock::now();

            for (const auto& passSignal : passSignals)
            {
                Config passConfig = config;
                passConfig.signalType = passSignal;
                if (passSignal.equalsIgnoreCase("timing"))
                    passConfig.seconds = 0.5; // Impulse, latency and settling need a short capture only.
                else if (passSignal.equalsIgnoreCase("residual"))
                    passConfig.seconds = 1.0; // 65k broadband samples are enough for aligned residual evidence.
                else if (passSignal.equalsIgnoreCase("boundary"))
                    passConfig.seconds = 10.0; // Includes five seconds of infrasonic and five seconds of ultrasonic stimulus.
                else if (passSignal.equalsIgnoreCase("stereo"))
                    passConfig.seconds = 4.0;
                else if (passSignal.equalsIgnoreCase("summing"))
                    passConfig.seconds = 3.0;
                else if (passSignal.equalsIgnoreCase("alias"))
                    passConfig.seconds = 5.0; // Ten stepped tones, 0.5 seconds per frequency.

                // Only create analysers that are valid for this signal. The user's
                // checkbox choices are preserved; this simply avoids meaningless
                // warnings and empty analysers during All mode.
                std::vector<juce::String> validAnalyzers;
                for (const auto& analyzerName : config.analyzers)
                {
                    const bool dedicatedPassAnalyzer = analyzerName.equalsIgnoreCase("Interaction")
                                                    || analyzerName.equalsIgnoreCase("Timing")
                                                    || analyzerName.equalsIgnoreCase("Residual")
                                                    || analyzerName.equalsIgnoreCase("Boundary")
                                                    || analyzerName.equalsIgnoreCase("Stereo")
                                                    || analyzerName.equalsIgnoreCase("Summing")
                                                    || analyzerName.equalsIgnoreCase("Alias");
                    if (dedicatedPassAnalyzer)
                        continue;
                    if (analyzerName.equalsIgnoreCase("Thd") && !passSignal.equalsIgnoreCase("sine"))
                        continue;
                    if (analyzerName.equalsIgnoreCase("LinearResponse") && passSignal.equalsIgnoreCase("sine"))
                        continue;
                    // In All mode, raw sample pairs are only needed for the sine
                    // behaviour/temporal pass. Keeping raw noise and sweep data in
                    // memory would multiply RAM use without improving the summary.
                    // Individual Noise or Sweep mode still permits Raw exports.
                    if (config.signalType.equalsIgnoreCase("all")
                        && analyzerName.equalsIgnoreCase("RawCsv")
                        && !passSignal.equalsIgnoreCase("sine"))
                        continue;
                    validAnalyzers.push_back(analyzerName);
                }
                passConfig.analyzers = std::move(validAnalyzers);

                juce::MessageManager::callAsync([this, passSignal]() {
                    progressLabel.setText("Analysing " + passSignal + "...", juce::dontSendNotification);
                });

                // Start every signal pass with a fresh plugin instance so state from
                // the sine pass cannot leak into sweep or noise measurements.
                juce::File passPluginFile(passConfig.pluginPath);
                juce::String passError;
                auto passPlugin = loadPluginInstance(passPluginFile, passConfig.sampleRate,
                                                     passConfig.blockSize, passError);
                if (passPlugin == nullptr)
                    throw std::runtime_error(("Failed to create plugin instance for " + passSignal +
                                              " pass: " + passError).toStdString());
                if (primaryState != nullptr && primaryState->getSize() > 0)
                    passPlugin->setStateInformation(primaryState->getData(), static_cast<int>(primaryState->getSize()));

                std::vector<std::unique_ptr<juce::AudioPluginInstance>> passSerialOwned;
                std::vector<juce::AudioPluginInstance*> passSerialPlugins;
                if (serialPath.isNotEmpty())
                {
                    juce::String serialError;
                    auto serial = loadPluginInstance(juce::File(serialPath), passConfig.sampleRate, passConfig.blockSize, serialError);
                    if (serial == nullptr)
                        throw std::runtime_error(("Failed to create Plugin 2 for " + passSignal + ": " + serialError).toStdString());
                    if (serialState != nullptr && serialState->getSize() > 0)
                        serial->setStateInformation(serialState->getData(), static_cast<int>(serialState->getSize()));
                    passSerialPlugins.push_back(serial.get());
                    passSerialOwned.push_back(std::move(serial));
                }
                for (size_t i = 0; i < additionalSerialPaths.size(); ++i)
                {
                    juce::String extraError;
                    auto extra = loadPluginInstance(juce::File(additionalSerialPaths[i]), passConfig.sampleRate, passConfig.blockSize, extraError);
                    if (extra == nullptr)
                        throw std::runtime_error(("Failed to create Plugin " + juce::String((int)i + 3) + " for " + passSignal + ": " + extraError).toStdString());
                    if (i < additionalSerialStates.size() && additionalSerialStates[i] != nullptr && additionalSerialStates[i]->getSize() > 0)
                        extra->setStateInformation(additionalSerialStates[i]->getData(), static_cast<int>(additionalSerialStates[i]->getSize()));
                    passSerialPlugins.push_back(extra.get());
                    passSerialOwned.push_back(std::move(extra));
                }

                auto passAnalyzers = createAnalyzers(passConfig, juce::File{}, paramNames);
                const int64_t totalSamples = (int64_t)(passConfig.seconds * passConfig.sampleRate);

                runMeasurementGrid(
                    *passPlugin, passConfig.sampleRate, passConfig.blockSize, totalSamples,
                    runs, passAnalyzers, passConfig, juce::File{},
                    [this, &completedPassRuns, totalPassRuns, overallStartTime, passSignal](int runIndex) {
                        const int absoluteRun = completedPassRuns + runIndex + 1;
                        const double overallProgress = totalPassRuns > 0
                            ? (double)absoluteRun / (double)totalPassRuns : 0.0;

                        auto currentTime = std::chrono::steady_clock::now();
                        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                            currentTime - overallStartTime).count();

                        juce::String statusText = passSignal.toUpperCase() + "  " +
                            juce::String(absoluteRun) + " / " + juce::String(totalPassRuns);

                        if (absoluteRun > 1 && elapsed > 0)
                        {
                            const double runsPerSecond = (double)absoluteRun / (double)elapsed;
                            const int remainingRuns = totalPassRuns - absoluteRun;
                            const int secondsRemaining = runsPerSecond > 0.0
                                ? (int)(remainingRuns / runsPerSecond) : 0;
                            if (secondsRemaining >= 60)
                                statusText += " (~" + juce::String(secondsRemaining / 60) + "m remaining)";
                            else
                                statusText += " (~" + juce::String(secondsRemaining) + "s remaining)";
                        }

                        juce::MessageManager::callAsync([this, statusText, overallProgress]() {
                            progressLabel.setText(statusText, juce::dontSendNotification);
                            progress = overallProgress;
                            progressBar.repaint();
                        });
                    }, passSerialPlugins, [this]() { return measurementCancelRequested.load(); });

                if (measurementCancelRequested.load()) break;
                completedPassRuns += (int)runs.size();

                for (auto& analyzer : passAnalyzers)
                    analyzer->finish(juce::File{});

                // Retain only the datasets needed by the current summary. All CSV
                // exports still remain on disk for every pass.
                for (const auto& analyzer : passAnalyzers)
                {
                    if (passSignal.equalsIgnoreCase("sine"))
                    {
                        if (const auto* raw = dynamic_cast<const RawCsvAnalyzer*>(analyzer.get()))
                        {
                            const auto& rawResult = raw->getResult();
                            haveSineRaw = !rawResult.isEmpty();
                            if (haveSineRaw)
                            {
                                waveformSummary = makeWaveformSummary(rawResult);
                                temporalSummary = makeTemporalSummary(rawResult);
                            }
                        }
                        else if (const auto* rms = dynamic_cast<const RmsPeakAnalyzer*>(analyzer.get()))
                        {
                            sineRmsResult = rms->getResult();
                            haveSineRms = !sineRmsResult.isEmpty();
                        }
                        else if (const auto* thd = dynamic_cast<const ThdAnalyzer*>(analyzer.get()))
                        {
                            sineThdResult = thd->getResult();
                            haveSineThd = !sineThdResult.isEmpty();
                        }
                        else if (const auto* transfer = dynamic_cast<const TransferCurveAnalyzer*>(analyzer.get()))
                        {
                            sineTransferResult = transfer->getResult();
                            haveSineTransfer = !sineTransferResult.isEmpty();
                        }
                    }
                    else if (passSignal.equalsIgnoreCase("sweep"))
                    {
                        if (const auto* linear = dynamic_cast<const LinearResponseAnalyzer*>(analyzer.get()))
                        {
                            sweepLinearResult = linear->getResult();
                            haveSweepLinear = !sweepLinearResult.isEmpty();
                            sweepHarmonicFingerprintResult = linear->getHarmonicFingerprintResult();
                            haveSweepHarmonicFingerprint = !sweepHarmonicFingerprintResult.isEmpty();
                        }
                    }
                    else if (passSignal.equalsIgnoreCase("noise"))
                    {
                        if (const auto* linear = dynamic_cast<const LinearResponseAnalyzer*>(analyzer.get()))
                        {
                            noiseLinearResult = linear->getResult();
                            haveNoiseLinear = !noiseLinearResult.isEmpty();
                        }
                    }
                    else if (passSignal.equalsIgnoreCase("interaction"))
                    {
                        if (const auto* interaction = dynamic_cast<const InteractionAnalyzer*>(analyzer.get()))
                        {
                            interactionResult = interaction->getResult();
                            haveInteraction = !interactionResult.isEmpty();
                        }
                    }
                    else if (passSignal.equalsIgnoreCase("timing"))
                    {
                        if (const auto* timing = dynamic_cast<const TimingAnalyzer*>(analyzer.get()))
                        {
                            timingResult = timing->getResult();
                            haveTiming = !timingResult.isEmpty();
                        }
                    }
                    else if (passSignal.equalsIgnoreCase("residual"))
                    {
                        if (const auto* residual = dynamic_cast<const ResidualAnalyzer*>(analyzer.get()))
                        {
                            residualResult = residual->getResult();
                            haveResidual = !residualResult.isEmpty();
                        }
                    }
                    else if (passSignal.equalsIgnoreCase("boundary"))
                    {
                        if (const auto* boundary = dynamic_cast<const BoundaryAnalyzer*>(analyzer.get()))
                        {
                            boundaryResult = boundary->getResult();
                            haveBoundary = !boundaryResult.isEmpty();
                        }
                    }
                    else if (passSignal.equalsIgnoreCase("stereo"))
                    {
                        if (const auto* stereo = dynamic_cast<const StereoAnalyzer*>(analyzer.get()))
                        { stereoResult = stereo->getResult(); haveStereo = !stereoResult.isEmpty(); }
                    }
                    else if (passSignal.equalsIgnoreCase("summing"))
                    {
                        if (const auto* summing = dynamic_cast<const SummingAnalyzer*>(analyzer.get()))
                        { summingResult = summing->getResult(); haveSumming = !summingResult.isEmpty(); }
                    }
                    else if (passSignal.equalsIgnoreCase("alias"))
                    {
                        if (const auto* alias = dynamic_cast<const AliasAnalyzer*>(analyzer.get()))
                        { aliasResult = alias->getResult(); haveAlias = !aliasResult.isEmpty(); }
                    }
                }

                // Move every completed in-memory dataset into the pending export
                // package. No CSV is written until Export Data is pressed.
                for (auto& analyzer : passAnalyzers)
                {
                    if (auto* raw = dynamic_cast<RawCsvAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ juce::String("grid_raw_") + passSignal.toLowerCase() + ".csv", raw->takeResult() });
                    else if (auto* rms = dynamic_cast<RmsPeakAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ juce::String("grid_rms_peak_") + passSignal.toLowerCase() + ".csv", rms->takeResult() });
                    else if (auto* transfer = dynamic_cast<TransferCurveAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ juce::String("grid_transfer_curves_") + passSignal.toLowerCase() + ".csv", transfer->takeResult() });
                    else if (auto* linear = dynamic_cast<LinearResponseAnalyzer*>(analyzer.get()))
                    {
                        auto harmonicFingerprint = linear->takeHarmonicFingerprintResult();
                        if (!harmonicFingerprint.isEmpty())
                            completedExports.emplace_back(PendingExportDataset{
                                juce::String("grid_harmonic_fingerprint_") + passSignal.toLowerCase() + ".csv",
                                std::move(harmonicFingerprint) });
                        completedExports.emplace_back(PendingExportDataset{
                            juce::String("grid_linear_response_") + passSignal.toLowerCase() + ".csv",
                            linear->takeResult() });
                    }
                    else if (auto* interaction = dynamic_cast<InteractionAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ "grid_interaction_dna.csv", interaction->takeResult() });
                    else if (auto* timing = dynamic_cast<TimingAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ "grid_timing_dna.csv", timing->takeResult() });
                    else if (auto* residual = dynamic_cast<ResidualAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ "grid_residual_dna.csv", residual->takeResult() });
                    else if (auto* boundary = dynamic_cast<BoundaryAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ "grid_boundary_dna.csv", boundary->takeResult() });
                    else if (auto* stereo = dynamic_cast<StereoAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ "grid_stereo_dna.csv", stereo->takeResult() });
                    else if (auto* summing = dynamic_cast<SummingAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ "grid_summing_dna.csv", summing->takeResult() });
                    else if (auto* alias = dynamic_cast<AliasAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ "grid_alias_stress.csv", alias->takeResult() });
                    else if (auto* thd = dynamic_cast<ThdAnalyzer*>(analyzer.get()))
                        completedExports.emplace_back(PendingExportDataset{ juce::String("grid_thd_") + passSignal.toLowerCase() + ".csv", thd->takeResult() });
                }
            }

            if (measurementCancelRequested.load()) {
                juce::MessageManager::callAsync([this]() {
                    progressLabel.setText("Measurement stopped", juce::dontSendNotification);
                    runMeasurementButton.setEnabled(true); stopMeasurementButton.setEnabled(false);
                    progress = 0.0; progressBar.repaint();
                });
                return;
            }

            LevelSummary levelSummary {
                "Run RMS/Peak to see level behaviour",
                "Run RMS/Peak to see transient behaviour",
                "Run RMS/Peak to see dynamic behaviour"
            };

            IdentitySummary identitySummary {
                "PRIMARY ROLE  waiting for complete measurement",
                "WHY  waiting for evidence",
                "WATCH OUT  waiting for evidence"
            };

            HarmonicSummary harmonicSummary {
                "HARMONIC DNA  Run THD with a sine signal",
                "BALANCE  waiting for harmonic measurement",
                "WHEN PUSHED  waiting for harmonic measurement"
            };

            BehaviourSummary behaviourSummary {
                "BEHAVIOUR  Run matching levels to measure response",
                "CHANGES WITH DRIVE  unavailable"
            };

            OperatingRangeSummary operatingRangeSummary {
                "OPERATING RANGE  Run -18, -12, -6 and -3 dBFS",
                "CHARACTER STARTS  unavailable",
                "UNDER STRESS  unavailable"
            };

            const MeasurementDataset* rmsResult = nullptr;
            const MeasurementDataset* thdResult = nullptr;
            const MeasurementDataset* linearResult = nullptr;
            const MeasurementDataset* transferResult = nullptr;

            NonlinearitySummary nonlinearitySummary {
                "NONLINEARITY  Run Transfer Curve to classify distortion",
                "CURVE BEHAVIOUR  unavailable"
            };

            ToneSummary toneSummary {
                "Run Linear Response with noise or sweep to detect the tone shape",
                "Run Linear Response with noise or sweep to inspect the low end",
                "Run Linear Response with noise or sweep to inspect the high end",
                "Run Linear Response with noise or sweep to find the strongest feature"
            };

            if (haveSineRms)
            {
                levelSummary = makeLevelSummary(sineRmsResult);
                rmsResult = &sineRmsResult;
            }

            if (haveSineThd)
            {
                harmonicSummary = makeHarmonicSummary(sineThdResult);
                thdResult = &sineThdResult;
            }

            if (haveSineTransfer)
            {
                transferResult = &sineTransferResult;
                nonlinearitySummary = makeNonlinearitySummary(sineTransferResult);
            }

            // Sweep is the preferred deterministic tone pass. Noise remains a
            // broadband fallback and is still exported for future texture rules.
            if (haveSweepLinear)
            {
                toneSummary = makeToneSummary(sweepLinearResult);
                linearResult = &sweepLinearResult;
            }
            else if (haveNoiseLinear)
            {
                toneSummary = makeToneSummary(noiseLinearResult);
                linearResult = &noiseLinearResult;
            }

            behaviourSummary = makeBehaviourSummary(rmsResult, thdResult, linearResult);
            operatingRangeSummary = makeOperatingRangeSummary(rmsResult, thdResult, linearResult);
            identitySummary = makeIdentitySummary(levelSummary, harmonicSummary, toneSummary, behaviourSummary,
                                                  operatingRangeSummary, nonlinearitySummary, waveformSummary,
                                                  temporalSummary);
            const auto fingerprintEvidence = haveSweepHarmonicFingerprint
                ? makeHarmonicFingerprintEvidence(sweepHarmonicFingerprintResult)
                : HarmonicFingerprintEvidence{};
            const auto producerSummary = makeProducerVocabularySummary(identitySummary, levelSummary,
                                                                             harmonicSummary, toneSummary,
                                                                             behaviourSummary, temporalSummary,
                                                                             fingerprintEvidence);

            const auto compactSummary = juce::String("PluginDNA Summary\n\n")
                + identitySummary.role + "\n"
                + producerSummary.character + "\n"
                + identitySummary.watch + "\n"
                + levelSummary.level + "\n"
                + producerSummary.behaviour;

            {
                std::lock_guard<std::mutex> lock(pendingExportMutex);
                pendingExportDatasets = std::move(completedExports);
                pendingEvidenceConfig = config;
                pendingEvidenceSummary = compactSummary;
                hasPendingEvidencePack = true;
            }

            juce::MessageManager::callAsync([this, identitySummary, producerSummary, levelSummary]() {
                progressLabel.setText(
                    "Measurement complete",
                    juce::dontSendNotification);
                identityRoleLabel.setText(identitySummary.role, juce::dontSendNotification);
                identityWhyLabel.setText(producerSummary.character, juce::dontSendNotification);
                identityWatchLabel.setText(identitySummary.watch, juce::dontSendNotification);
                levelSummaryLabel.setText(levelSummary.level, juce::dontSendNotification);
                peakSummaryLabel.setText({}, juce::dontSendNotification);
                dynamicsSummaryLabel.setText({}, juce::dontSendNotification);
                harmonicSummaryLabel.setText(producerSummary.behaviour, juce::dontSendNotification);
                harmonicBalanceLabel.setText({}, juce::dontSendNotification);
                harmonicGrowthLabel.setText({}, juce::dontSendNotification);
                toneBassLabel.setText({}, juce::dontSendNotification);
                toneMidLabel.setText({}, juce::dontSendNotification);
                toneTrebleLabel.setText({}, juce::dontSendNotification);
                toneLargestLabel.setText({}, juce::dontSendNotification);
                behaviourSummaryLabel.setText({}, juce::dontSendNotification);
                behaviourChangesLabel.setText({}, juce::dontSendNotification);
                operatingRangeLabel.setText({}, juce::dontSendNotification);
                characterStartsLabel.setText({}, juce::dontSendNotification);
                underStressLabel.setText({}, juce::dontSendNotification);
                nonlinearityLabel.setText({}, juce::dontSendNotification);
                curveBehaviourLabel.setText({}, juce::dontSendNotification);
                waveformShapeLabel.setText({}, juce::dontSendNotification);
                waveformKneeLabel.setText({}, juce::dontSendNotification);
                waveformStyleLabel.setText({}, juce::dontSendNotification);
                temporalResponseLabel.setText({}, juce::dontSendNotification);
                temporalAttackLabel.setText({}, juce::dontSendNotification);
                temporalAfterLabel.setText({}, juce::dontSendNotification);
                progress = 1.0;
                progressBar.repaint();
                runMeasurementButton.setEnabled(true);
                stopMeasurementButton.setEnabled(false);
                copySummaryButton.setEnabled(true);
                exportDataButton.setEnabled(true);
                exportEvidenceButton.setEnabled(true);
            });
        } catch (const std::exception& e) {
            std::cerr << "[Measurement] Exception: " << e.what() << std::endl;
            juce::MessageManager::callAsync([this, e]() {
                showError("Error: " + juce::String(e.what()));
                runMeasurementButton.setEnabled(true); stopMeasurementButton.setEnabled(false);
            });
        } catch (...) {
            std::cerr << "[Measurement] Unknown exception occurred" << std::endl;
            juce::MessageManager::callAsync([this]() {
                showError("Unknown error occurred during measurement");
                runMeasurementButton.setEnabled(true); stopMeasurementButton.setEnabled(false);
            });
        }
    }).detach();
}

Config MainComponent::buildConfigFromUI() {
    Config config;
    measurementConfig->fillConfig(config);
    saveCurrentParameterPanelState();

    auto appendBuckets = [&config](const std::vector<ParameterBucketConfig>& buckets,
                                   const std::map<juce::String, juce::AudioProcessorParameter*>& map,
                                   const juce::String& prefix) {
        for (auto bucket : buckets) {
            const auto rawName = bucket.paramName;
            const auto it = map.find(rawName);
            if (it != map.end() && it->second != nullptr)
                bucket.pluginDefaultValue = it->second->getDefaultValue();
            bucket.includePluginDefault = true;
            bucket.paramName = prefix + rawName;
            config.parameterBuckets.push_back(bucket);
        }
    };
    appendBuckets(primarySavedBuckets, primaryParameterMap, "P1::");
    appendBuckets(serialSavedBuckets, serialParameterMap, "P2::");
    for (const auto& stage : additionalSerialStages) {
        if (stage->pluginInstance == nullptr) break;
        appendBuckets(stage->savedBuckets, stage->parameterMap, "P" + juce::String(stage->stageNumber) + "::");
    }
    return config;
}

void MainComponent::showError(const juce::String& message) {
    juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon, "Error", message);
    progressLabel.setText(message, juce::dontSendNotification);
}

juce::String MainComponent::buildSummaryText() const
{
    juce::StringArray lines;
    lines.add("PluginDNA Summary");
    lines.add({});

    for (const auto* label : { &identityRoleLabel, &identityWhyLabel, &identityWatchLabel,
                               &levelSummaryLabel, &peakSummaryLabel, &dynamicsSummaryLabel,
                               &harmonicSummaryLabel, &harmonicBalanceLabel, &harmonicGrowthLabel,
                               &toneBassLabel, &toneMidLabel, &toneTrebleLabel, &toneLargestLabel,
                               &behaviourSummaryLabel, &behaviourChangesLabel,
                               &operatingRangeLabel, &characterStartsLabel, &underStressLabel,
                               &nonlinearityLabel, &curveBehaviourLabel,
                               &waveformShapeLabel, &waveformKneeLabel, &waveformStyleLabel,
                               &temporalResponseLabel, &temporalAttackLabel, &temporalAfterLabel })
    {
        const auto text = label->getText().trim();
        if (text.isNotEmpty())
            lines.add(text);
    }

    return lines.joinIntoString("\n");
}


void MainComponent::exportEvidencePack(const juce::File& overrideOutputDirectory)
{
    const bool exportingFullPackage = overrideOutputDirectory.getFullPathName().isNotEmpty();
    const auto outputPath = exportingFullPackage
                          ? overrideOutputDirectory.getFullPathName()
                          : outputPathEditor.getText().trim();
    if (outputPath.isEmpty())
    {
        showError("Please specify an output path");
        return;
    }

    juce::File outDir(outputPath);
    if (!outDir.exists() && !outDir.createDirectory())
    {
        showError("Could not create output directory");
        return;
    }

    std::vector<PendingExportDataset> datasets;
    juce::StringArray fullPackageFiles;
    Config config;
    juce::String summary;
    {
        std::lock_guard<std::mutex> lock(pendingExportMutex);
        if (!hasPendingEvidencePack || pendingExportDatasets.empty())
        {
            showError("Run a measurement before exporting evidence");
            return;
        }
        config = pendingEvidenceConfig;
        summary = pendingEvidenceSummary;
        for (const auto& item : pendingExportDatasets)
        {
            fullPackageFiles.addIfNotAlreadyThere(item.filename);
            if (!item.filename.startsWithIgnoreCase("grid_raw_"))
                datasets.push_back(item);
        }
    }

    exportEvidenceButton.setEnabled(false);
    exportDataButton.setEnabled(false);
    runMeasurementButton.setEnabled(false);
    progressLabel.setText("Building compact evidence pack...", juce::dontSendNotification);

    const auto pluginPath = pluginPathEditor.getText().trim();
    const auto pluginInfo = pluginInfoLabel.getText();
    const auto serialPluginPath = serialPluginPathEditor.getText().trim();
    const auto serialPluginInfo = serialPluginInfoLabel.getText();
    std::vector<juce::String> additionalPluginPaths;
    std::vector<juce::String> additionalPluginInfos;
    for (const auto& stage : additionalSerialStages)
    {
        if (stage->pluginInstance == nullptr) break;
        additionalPluginPaths.push_back(stage->pathEditor->getText().trim());
        additionalPluginInfos.push_back(stage->infoLabel->getText());
    }
    std::thread([this, outDir, config, summary, datasets = std::move(datasets), fullPackageFiles, pluginPath, pluginInfo, serialPluginPath, serialPluginInfo, additionalPluginPaths, additionalPluginInfos, exportingFullPackage]() mutable {
        auto* root = new juce::DynamicObject();
        root->setProperty("schema", "PluginDNA Evidence Pack");
        root->setProperty("schemaVersion", 24);
        root->setProperty("createdUtc", juce::Time::getCurrentTime().toISO8601(true));

        const juce::File pluginFile(pluginPath);
        auto* plugin = new juce::DynamicObject();
        plugin->setProperty("name", pluginFile.getFileNameWithoutExtension());
        plugin->setProperty("path", pluginPath);
        plugin->setProperty("info", pluginInfo);
        root->setProperty("plugin", juce::var(plugin));
        juce::Array<juce::var> chain;
        {
            auto* stage1 = new juce::DynamicObject();
            stage1->setProperty("position", 1); stage1->setProperty("name", pluginFile.getFileNameWithoutExtension());
            stage1->setProperty("path", pluginPath); stage1->setProperty("info", pluginInfo);
            chain.add(juce::var(stage1));
        }
        if (serialPluginPath.isNotEmpty())
        {
            auto* serial = new juce::DynamicObject();
            serial->setProperty("name", juce::File(serialPluginPath).getFileNameWithoutExtension());
            serial->setProperty("path", serialPluginPath); serial->setProperty("info", serialPluginInfo);
            serial->setProperty("position", 2);
            root->setProperty("serialPlugin", juce::var(serial)); // compatibility with schema 23 readers
            auto* stage2 = new juce::DynamicObject();
            stage2->setProperty("position", 2); stage2->setProperty("name", juce::File(serialPluginPath).getFileNameWithoutExtension());
            stage2->setProperty("path", serialPluginPath); stage2->setProperty("info", serialPluginInfo);
            chain.add(juce::var(stage2));
        }
        for (size_t i = 0; i < additionalPluginPaths.size(); ++i)
        {
            auto* extra = new juce::DynamicObject();
            extra->setProperty("position", (int)i + 3);
            extra->setProperty("name", juce::File(additionalPluginPaths[i]).getFileNameWithoutExtension());
            extra->setProperty("path", additionalPluginPaths[i]);
            if (i < additionalPluginInfos.size()) extra->setProperty("info", additionalPluginInfos[i]);
            chain.add(juce::var(extra));
        }
        root->setProperty("pluginChain", chain);
        const int chainCount = 1 + (serialPluginPath.isNotEmpty() ? 1 : 0) + (int)additionalPluginPaths.size();
        root->setProperty("systemType", chainCount > 1 ? "serial plugin chain" : "single plugin");

        auto* test = new juce::DynamicObject();
        test->setProperty("sampleRate", config.sampleRate);
        test->setProperty("durationSeconds", config.seconds);
        test->setProperty("blockSize", config.blockSize);
        test->setProperty("analysisMode", config.signalType);
        test->setProperty("sineFrequencyHz", config.sineFrequency);
        test->setProperty("sweepStartHz", config.sweepStartHz);
        test->setProperty("sweepEndHz", config.sweepEndHz);
        juce::Array<juce::var> aliasFrequencies;
        for (const auto frequency : AliasAnalyzer::testFrequencies()) aliasFrequencies.add(frequency);
        test->setProperty("aliasStressFrequenciesHz", aliasFrequencies);
        test->setProperty("aliasDetectionFloorDb", -90.0);
        juce::Array<juce::var> gainBuckets;
        for (const auto gain : config.inputGainBucketsDb)
            gainBuckets.add(gain);
        test->setProperty("inputGainBucketsDb", gainBuckets);
        juce::Array<juce::var> enabledAnalyzers;
        for (const auto& analyzer : config.analyzers) enabledAnalyzers.add(analyzer);
        test->setProperty("enabledAnalyzers", enabledAnalyzers);
        root->setProperty("testConfiguration", juce::var(test));

        juce::Array<juce::var> parameters;
        for (const auto& bucket : config.parameterBuckets)
        {
            auto* parameter = new juce::DynamicObject();
            parameter->setProperty("name", bucket.paramName);
            parameter->setProperty("strategy", bucket.strategy);
            const auto values = generatedBucketValues(bucket);
            parameter->setProperty("hostRangeMin", bucket.min);
            parameter->setProperty("hostRangeMax", bucket.max);
            parameter->setProperty("pluginDefaultValue", bucket.pluginDefaultValue);
            parameter->setProperty("defaultIncludedInScan", bucket.includePluginDefault && !bucket.strategy.equalsIgnoreCase("Enumerated"));
            parameter->setProperty("testedValues", values);
            if (bucket.strategy == "Enumerated")
            {
                juce::Array<juce::var> states;
                for (size_t i = 0; i < bucket.values.size(); ++i)
                {
                    auto* state = new juce::DynamicObject();
                    state->setProperty("value", bucket.values[i]);
                    state->setProperty("label", i < bucket.valueLabels.size() ? bucket.valueLabels[i]
                                                                                 : juce::String("State ") + juce::String(static_cast<int>(i) + 1));
                    states.add(juce::var(state));
                }
                parameter->setProperty("states", states);
                parameter->setProperty("parameterType", "enumerated mode selector");
            }
            else
            {
                parameter->setProperty("parameterType", "continuous");
            }
            parameter->setProperty("referencePolicy", "actual plugin default");
            parameter->setProperty("referenceValue", bucket.pluginDefaultValue);
            parameters.add(juce::var(parameter));
        }
        root->setProperty("parameters", parameters);
        root->setProperty("parameterScanEnabled", !config.parameterBuckets.empty());
        root->setProperty("scanType", config.parameterBuckets.size() == 0 ? "baseline current-state test"
                                      : config.parameterBuckets.size() == 1 ? "parameter profile"
                                      : config.parameterBuckets.size() == 2 ? "parameter coupling"
                                      : "multi-parameter research grid");
        root->setProperty("producerSummary", summary);

        juce::Array<juce::var> evidence;
        for (const auto& item : datasets)
        {
            if (item.filename.startsWithIgnoreCase("grid_raw_"))
                continue;
            evidence.add(summariseDatasetForEvidence(item.filename, item.dataset, config));
        }
        root->setProperty("evidence", evidence);

        // Compact Evidence v2: one directly readable fingerprint per tested operating point.
        // This is assembled from the derived analyser evidence above; no extra DSP pass is run.
        const auto operatingPointDNA = buildOperatingPointDNA(evidence);
        root->setProperty("operatingPointDNA", operatingPointDNA);
        root->setProperty("compactThresholds", buildCompactThresholds(operatingPointDNA));

        // Phase 2: Feature Extraction v2. Extracts a rich objective feature
        // set from the evidence already in memory. No additional test passes.
        {
            juce::Array<juce::var> analyserFeatures;

            const auto col = [](const MeasurementDataset& data, const juce::String& name)
            {
                return findColumn(data, name.toStdString());
            };
            const auto valueAt = [](const std::vector<double>& row, int c)
            {
                return c >= 0 && static_cast<std::size_t>(c) < row.size()
                    ? row[static_cast<std::size_t>(c)]
                    : std::numeric_limits<double>::quiet_NaN();
            };
            const auto collect = [&](const MeasurementDataset& data, const juce::String& name,
                                     const std::function<bool(const std::vector<double>&)>& include)
            {
                std::vector<double> values;
                const int c = col(data, name);
                if (c < 0) return values;
                values.reserve(data.rows.size());
                for (const auto& row : data.rows)
                {
                    const double v = valueAt(row, c);
                    if (std::isfinite(v) && (!include || include(row))) values.push_back(v);
                }
                return values;
            };
            const auto mean = [](const std::vector<double>& values)
            {
                if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
                return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
            };
            const auto standardDeviation = [&](const std::vector<double>& values)
            {
                if (values.size() < 2) return 0.0;
                const double m = mean(values);
                double sum = 0.0;
                for (const double v : values) sum += (v - m) * (v - m);
                return std::sqrt(sum / static_cast<double>(values.size()));
            };
            const auto percentile = [](std::vector<double> values, double q)
            {
                if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
                std::sort(values.begin(), values.end());
                const double pos = juce::jlimit(0.0, 1.0, q) * static_cast<double>(values.size() - 1);
                const auto lo = static_cast<std::size_t>(std::floor(pos));
                const auto hi = static_cast<std::size_t>(std::ceil(pos));
                const double t = pos - static_cast<double>(lo);
                return values[lo] * (1.0 - t) + values[hi] * t;
            };
            const auto meanColumn = [&](const MeasurementDataset& data, const juce::String& name,
                                        const std::function<bool(const std::vector<double>&)>& include)
            {
                return mean(collect(data, name, include));
            };
            const auto allRows = [](const std::vector<double>&) { return true; };
            const auto addNumber = [](juce::DynamicObject* object, const juce::String& key, double value)
            {
                if (std::isfinite(value)) object->setProperty(key, value);
            };
            const auto addStats = [&](juce::DynamicObject* object, const juce::String& prefix,
                                      const std::vector<double>& values)
            {
                if (values.empty()) return;
                addNumber(object, prefix + "Mean", mean(values));
                addNumber(object, prefix + "StdDev", standardDeviation(values));
                addNumber(object, prefix + "Min", *std::min_element(values.begin(), values.end()));
                addNumber(object, prefix + "Max", *std::max_element(values.begin(), values.end()));
                addNumber(object, prefix + "P10", percentile(values, 0.10));
                addNumber(object, prefix + "P90", percentile(values, 0.90));
            };

            for (const auto& exportItem : datasets)
            {
                const auto& data = exportItem.dataset;
                auto* feature = new juce::DynamicObject();
                feature->setProperty("source", exportItem.filename);
                bool useful = false;

                if (exportItem.filename.containsIgnoreCase("grid_rms_peak"))
                {
                    feature->setProperty("family", "levelDynamics");
                    std::vector<double> levelChanges, peakChanges, crestInValues, crestOutValues, crestChanges;
                    std::vector<double> inputRmsDb, outputRmsDb, inputPeakDb, outputPeakDb;
                    const int ri = col(data, "rmsInL"), ro = col(data, "rmsOutL");
                    const int pi = col(data, "peakInL"), po = col(data, "peakOutL");
                    for (const auto& row : data.rows)
                    {
                        const double inRms = valueAt(row, ri), outRms = valueAt(row, ro);
                        const double inPeak = valueAt(row, pi), outPeak = valueAt(row, po);
                        if (!(inRms > 0.0 && outRms > 0.0 && inPeak > 0.0 && outPeak > 0.0)) continue;
                        const double inRmsDb = amplitudeToDb(inRms), outRmsDb = amplitudeToDb(outRms);
                        const double inPeakDb = amplitudeToDb(inPeak), outPeakDb = amplitudeToDb(outPeak);
                        const double crestIn = inPeakDb - inRmsDb;
                        const double crestOut = outPeakDb - outRmsDb;
                        inputRmsDb.push_back(inRmsDb); outputRmsDb.push_back(outRmsDb);
                        inputPeakDb.push_back(inPeakDb); outputPeakDb.push_back(outPeakDb);
                        levelChanges.push_back(outRmsDb - inRmsDb);
                        peakChanges.push_back(outPeakDb - inPeakDb);
                        crestInValues.push_back(crestIn); crestOutValues.push_back(crestOut);
                        crestChanges.push_back(crestOut - crestIn);
                    }
                    addStats(feature, "levelChangeDb", levelChanges);
                    addStats(feature, "peakChangeDb", peakChanges);
                    addStats(feature, "crestFactorInDb", crestInValues);
                    addStats(feature, "crestFactorOutDb", crestOutValues);
                    addStats(feature, "crestFactorChangeDb", crestChanges);
                    if (!inputRmsDb.empty())
                    {
                        const double xMean = mean(inputRmsDb), yMean = mean(outputRmsDb);
                        double cov = 0.0, var = 0.0;
                        for (std::size_t i = 0; i < inputRmsDb.size(); ++i)
                        { cov += (inputRmsDb[i]-xMean)*(outputRmsDb[i]-yMean); var += (inputRmsDb[i]-xMean)*(inputRmsDb[i]-xMean); }
                        if (var > 1.0e-12) addNumber(feature, "rmsTransferSlope", cov / var);
                    }
                    addNumber(feature, "peakToRmsChangeCorrelation", [&]()
                    {
                        if (levelChanges.size() < 2) return std::numeric_limits<double>::quiet_NaN();
                        const double a = mean(levelChanges), b = mean(peakChanges);
                        double num=0.0, da=0.0, db=0.0;
                        for (std::size_t i=0;i<levelChanges.size();++i)
                        { const double xa=levelChanges[i]-a, xb=peakChanges[i]-b; num+=xa*xb; da+=xa*xa; db+=xb*xb; }
                        return da>0.0 && db>0.0 ? num/std::sqrt(da*db) : 0.0;
                    }());
                    useful = !levelChanges.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_alias_stress"))
                {
                    feature->setProperty("family", "aliasStress");
                    const auto aliasRatios = collect(data, "aliasPowerRatio", allRows);
                    const auto aliasDb = collect(data, "aliasPowerDb", allRows);
                    const auto detected = collect(data, "detectedFoldCount", allRows);
                    const auto integrity = collect(data, "ladderIntegrity", allRows);
                    addStats(feature, "aliasPowerRatio", aliasRatios);
                    addStats(feature, "aliasPowerDb", aliasDb);
                    addStats(feature, "detectedFoldCount", detected);
                    addStats(feature, "ladderIntegrity", integrity);

                    const int fcol = col(data, "fundamentalHz");
                    const int acol = col(data, "aliasPowerDb");
                    double firstDetected = std::numeric_limits<double>::quiet_NaN();
                    double highestClean = std::numeric_limits<double>::quiet_NaN();
                    double strongestDb = -300.0;
                    double strongestAt = 0.0;
                    for (const auto& row : data.rows)
                    {
                        const double f = valueAt(row, fcol);
                        const double a = valueAt(row, acol);
                        if (!std::isfinite(f) || !std::isfinite(a)) continue;
                        if (a >= -90.0 && !std::isfinite(firstDetected)) firstDetected = f;
                        if (a < -90.0) highestClean = std::isfinite(highestClean) ? std::max(highestClean, f) : f;
                        if (a > strongestDb) { strongestDb = a; strongestAt = f; }
                    }
                    addNumber(feature, "firstDetectedAliasFundamentalHz", firstDetected);
                    addNumber(feature, "highestCleanFundamentalHz", highestClean);
                    addNumber(feature, "strongestAliasPowerDb", strongestDb);
                    addNumber(feature, "strongestAliasAtFundamentalHz", strongestAt);
                    useful = !aliasRatios.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_thd_sine"))
                {
                    feature->setProperty("family", "harmonicLadder");
                    std::vector<double> thdDbValues, oddEvenDbValues, ladderSlopeValues, ladderSmoothnessValues;
                    std::array<std::vector<double>, 11> harmonics;
                    for (const auto& row : data.rows)
                    {
                        const double thd = valueAt(row, col(data, "thd"));
                        if (thd > 0.0) thdDbValues.push_back(amplitudeToDb(thd));
                        std::vector<std::pair<int,double>> active;
                        double oddPower=0.0, evenPower=0.0;
                        for (int h=2; h<=10; ++h)
                        {
                            const double v=valueAt(row,col(data,"h"+juce::String(h)));
                            if (!(v>0.0)) continue;
                            harmonics[static_cast<std::size_t>(h)].push_back(amplitudeToDb(v));
                            active.push_back({h,amplitudeToDb(v)});
                            if (h%2==0) evenPower+=v*v; else oddPower+=v*v;
                        }
                        if (oddPower>0.0 || evenPower>0.0)
                            oddEvenDbValues.push_back(10.0*std::log10(std::max(oddPower,1.0e-30)/std::max(evenPower,1.0e-30)));
                        if (active.size()>=3)
                        {
                            double sx=0,sy=0,sxx=0,sxy=0;
                            for (const auto& hv:active){ const double x=static_cast<double>(hv.first), y=hv.second; sx+=x;sy+=y;sxx+=x*x;sxy+=x*y; }
                            const double n=static_cast<double>(active.size()), den=n*sxx-sx*sx;
                            if (std::abs(den)>1.0e-12)
                            {
                                const double slope=(n*sxy-sx*sy)/den, intercept=(sy-slope*sx)/n;
                                ladderSlopeValues.push_back(slope);
                                double err=0.0; for(const auto& hv:active){const double d=hv.second-(intercept+slope*hv.first);err+=d*d;}
                                ladderSmoothnessValues.push_back(std::sqrt(err/n));
                            }
                        }
                    }
                    addStats(feature,"thdDb",thdDbValues);
                    addStats(feature,"oddVsEvenDb",oddEvenDbValues);
                    addStats(feature,"ladderSlopeDbPerOrder",ladderSlopeValues);
                    addStats(feature,"ladderDeviationDb",ladderSmoothnessValues);
                    for(int h=2;h<=10;++h) addStats(feature,"h"+juce::String(h)+"Db",harmonics[static_cast<std::size_t>(h)]);
                    const int gc=col(data,"inputGainDb");
                    std::vector<std::pair<double,double>> growth;
                    for(const double gain:config.inputGainBucketsDb)
                    {
                        const double v=meanColumn(data,"thd",[&](const std::vector<double>& row){return std::abs(valueAt(row,gc)-gain)<0.01;});
                        if(v>0.0) growth.push_back({gain,amplitudeToDb(v)});
                    }
                    if(growth.size()>=2)
                    {
                        addNumber(feature,"thdGrowthDbAcrossInputRange",growth.back().second-growth.front().second);
                        addNumber(feature,"thdGrowthDbPerInputDb",(growth.back().second-growth.front().second)/std::max(1.0,growth.back().first-growth.front().first));
                        double maxCurvature=0.0, kneeGain=growth.front().first;
                        for(std::size_t i=1;i+1<growth.size();++i)
                        {
                            const double curvature=std::abs((growth[i+1].second-growth[i].second)-(growth[i].second-growth[i-1].second));
                            if(curvature>maxCurvature){maxCurvature=curvature;kneeGain=growth[i].first;}
                        }
                        addNumber(feature,"thdGrowthCurvatureDb",maxCurvature);
                        addNumber(feature,"estimatedThdKneeInputDb",kneeGain);
                    }
                    // Stability across input level: low values mean the harmonic fingerprint
                    // retains its shape as the processor is driven.
                    std::vector<double> harmonicLevelVariation;
                    for (int h=2; h<=10; ++h)
                    {
                        std::vector<double> byGain;
                        for (const double gain : config.inputGainBucketsDb)
                        {
                            const double v = meanColumn(data, "h"+juce::String(h), [&](const std::vector<double>& row)
                            { return std::abs(valueAt(row,gc)-gain)<0.01; });
                            if (v > 0.0) byGain.push_back(amplitudeToDb(v));
                        }
                        if (byGain.size() >= 2) harmonicLevelVariation.push_back(standardDeviation(byGain));
                    }
                    if (!harmonicLevelVariation.empty())
                    {
                        const double variation = mean(harmonicLevelVariation);
                        addNumber(feature,"harmonicStabilityAcrossLevelDb",variation);
                        addNumber(feature,"harmonicStabilityAcrossLevel",1.0/(1.0+variation));
                    }
                    if (!oddEvenDbValues.empty())
                        addNumber(feature,"evenOddPowerRatioDb",-mean(oddEvenDbValues));
                    useful=!thdDbValues.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_transfer_curves"))
                {
                    feature->setProperty("family", "transferCurve");
                    const int xc=col(data,"x"), yc=col(data,"meanY");
                    std::vector<std::pair<double,double>> points;
                    for(const auto& row:data.rows){const double x=valueAt(row,xc),y=valueAt(row,yc);if(std::isfinite(x)&&std::isfinite(y))points.push_back({x,y});}
                    std::sort(points.begin(),points.end());
                    if(points.size()>=4)
                    {
                        double maxDev=0.0,maxSlope=-1.0e30,minSlope=1.0e30,positiveDev=0.0,negativeDev=0.0;
                        int posCount=0,negCount=0,reversals=0,inflections=0; double prevSlope=0.0,prevDelta=0.0; bool haveSlope=false;
                        for(std::size_t i=1;i<points.size();++i)
                        {
                            const double dx=points[i].first-points[i-1].first; if(std::abs(dx)<1.0e-12)continue;
                            const double slope=(points[i].second-points[i-1].second)/dx;
                            maxSlope=std::max(maxSlope,slope);minSlope=std::min(minSlope,slope);
                            if(haveSlope){if(prevSlope>0.02&&slope<-0.02)++reversals; const double delta=slope-prevSlope;if((prevDelta>0&&delta<0)||(prevDelta<0&&delta>0))++inflections;prevDelta=delta;}
                            else{prevDelta=0.0;haveSlope=true;} prevSlope=slope;
                        }
                        for(const auto& pt:points)
                        {
                            const double dev=pt.second-pt.first;maxDev=std::max(maxDev,std::abs(dev));
                            if(pt.first>0){positiveDev+=dev;++posCount;}else if(pt.first<0){negativeDev+=dev;++negCount;}
                        }
                        const double symmetryError=(posCount&&negCount)?std::abs(positiveDev/posCount+negativeDev/negCount):0.0;
                        std::vector<double> localSlopes;
                        for(std::size_t i=1;i<points.size();++i){const double dx=points[i].first-points[i-1].first;if(std::abs(dx)>1.0e-12)localSlopes.push_back((points[i].second-points[i-1].second)/dx);}
                        const double meanLocalSlope = mean(localSlopes);
                        double localSlopeAtZero = std::numeric_limits<double>::quiet_NaN(), closestZero = std::numeric_limits<double>::max();
                        for(std::size_t i=0;i<localSlopes.size();++i){const double midpoint=0.5*(points[i].first+points[i+1].first);if(std::abs(midpoint)<closestZero){closestZero=std::abs(midpoint);localSlopeAtZero=localSlopes[i];}}
                        addNumber(feature,"maxStaticDeviation",maxDev);addNumber(feature,"maximumLocalSlope",maxSlope);addNumber(feature,"minimumLocalSlope",minSlope);
                        addStats(feature,"localSlope",localSlopes);addNumber(feature,"localSlopeAtZero",localSlopeAtZero);
                        const double symmetryIndex=1.0/(1.0+symmetryError);
                        addNumber(feature,"curveSymmetryError",symmetryError);addNumber(feature,"curveAsymmetryScore",symmetryError/std::max(maxDev,1.0e-12));addNumber(feature,"symmetryIndex",symmetryIndex);
                        double deviationPower=0.0,gainErrorPower=0.0;
                        for(const auto& pt:points){const double nonlinearDeviation=pt.second-meanLocalSlope*pt.first;deviationPower+=nonlinearDeviation*nonlinearDeviation;const double gainError=(meanLocalSlope-1.0)*pt.first;gainErrorPower+=gainError*gainError;}
                        const double nonlinearRms=std::sqrt(deviationPower/points.size()),gainErrorRms=std::sqrt(gainErrorPower/points.size());
                        addNumber(feature,"saturationEfficiency",nonlinearRms/std::max(1.0e-15,nonlinearRms+gainErrorRms));
                        feature->setProperty("foldbackReversalCount",reversals);feature->setProperty("inflectionCount",inflections);
                        addNumber(feature,"outputAtPositiveLimit",points.back().second);addNumber(feature,"outputAtNegativeLimit",points.front().second);
                        addNumber(feature,"limitCompressionRatio",std::abs(points.back().second-points.front().second)/std::max(1.0e-12,std::abs(points.back().first-points.front().first)));
                        useful=true;
                    }
                }
                else if (exportItem.filename.containsIgnoreCase("linear_response"))
                {
                    feature->setProperty("family", "frequencyResponse");
                    const int fc=col(data,"freqHz"), mc=col(data,"magDb");
                    const auto band=[&](double lo,double hi){return meanColumn(data,"magDb",[&](const std::vector<double>& row){const double f=valueAt(row,fc);return f>=lo&&f<=hi;});};
                    const double sub=band(20,60),bass=band(60,200),lowMid=band(200,500),mid=band(500,2000),presence=band(2000,6000),high=band(6000,12000),air=band(12000,20000);
                    addNumber(feature,"subDb",sub);addNumber(feature,"bassDb",bass);addNumber(feature,"lowMidDb",lowMid);addNumber(feature,"midDb",mid);
                    addNumber(feature,"presenceDb",presence);addNumber(feature,"highDb",high);addNumber(feature,"airDb",air);
                    if(std::isfinite(mid)){addNumber(feature,"subVsMidDb",sub-mid);addNumber(feature,"bassVsMidDb",bass-mid);addNumber(feature,"lowMidVsMidDb",lowMid-mid);addNumber(feature,"presenceVsMidDb",presence-mid);addNumber(feature,"highVsMidDb",high-mid);addNumber(feature,"airVsMidDb",air-mid);}
                    std::vector<double> logF,mags;double weighted=0.0,totalWeight=0.0,maxMag=-1.0e30,minMag=1.0e30,maxFreq=0.0,minFreq=0.0;
                    for(const auto& row:data.rows){const double f=valueAt(row,fc),m=valueAt(row,mc);if(!(f>0.0&&std::isfinite(m)))continue;logF.push_back(std::log2(f));mags.push_back(m);const double w=std::pow(10.0,m/20.0);weighted+=f*w;totalWeight+=w;if(m>maxMag){maxMag=m;maxFreq=f;}if(m<minMag){minMag=m;minFreq=f;}}
                    if(logF.size()>=2)
                    {
                        const double xm=mean(logF),ym=mean(mags);double cov=0,var=0;for(std::size_t i=0;i<logF.size();++i){cov+=(logF[i]-xm)*(mags[i]-ym);var+=(logF[i]-xm)*(logF[i]-xm);}if(var>0)addNumber(feature,"overallTiltDbPerOctave",cov/var);
                        addNumber(feature,"responseRippleDb",maxMag-minMag);addNumber(feature,"strongestBoostDb",maxMag);addNumber(feature,"strongestBoostFrequencyHz",maxFreq);addNumber(feature,"deepestCutDb",minMag);addNumber(feature,"deepestCutFrequencyHz",minFreq);
                        if(totalWeight>0)addNumber(feature,"responseEnergyCentroidHz",weighted/totalWeight);
                    }
                    useful=!mags.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("harmonic_fingerprint"))
                {
                    feature->setProperty("family", "frequencyDependentHarmonics");
                    const int fc=col(data,"fundamentalHz"),hc=col(data,"harmonicOrder"),mc=col(data,"magDbRelativeToInput");
                    std::array<std::vector<double>,8> byOrder;std::vector<double> low,mid,high;
                    for(const auto& row:data.rows){const int h=static_cast<int>(std::llround(valueAt(row,hc)));const double f=valueAt(row,fc),m=valueAt(row,mc);if(h>=2&&h<=7&&std::isfinite(m)){byOrder[static_cast<std::size_t>(h)].push_back(m);if(f<200)low.push_back(m);else if(f<3000)mid.push_back(m);else high.push_back(m);}}
                    for(int h=2;h<=7;++h)addStats(feature,"h"+juce::String(h)+"Db",byOrder[static_cast<std::size_t>(h)]);
                    addStats(feature,"lowBandHarmonicDb",low);addStats(feature,"midBandHarmonicDb",mid);addStats(feature,"highBandHarmonicDb",high);
                    if(!byOrder[2].empty()&&!byOrder[3].empty())addNumber(feature,"h2VsH3Db",mean(byOrder[2])-mean(byOrder[3]));
                    double lowOrder=0,highOrder=0;for(int h=2;h<=7;++h){const double m=mean(byOrder[static_cast<std::size_t>(h)]);if(!std::isfinite(m))continue;const double a=std::pow(10.0,m/20.0);if(h<=3)lowOrder+=a*a;else highOrder+=a*a;}if(lowOrder>0)addNumber(feature,"highOrderVsLowOrderDb",10.0*std::log10(std::max(highOrder,1.0e-30)/lowOrder));
                    useful=!low.empty()||!mid.empty()||!high.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_residual_dna"))
                {
                    feature->setProperty("family", "residual");
                    const auto overall=collect(data,"residualDbRelativeToOutput",allRows),low=collect(data,"lowResidualDb",allRows),mid=collect(data,"midResidualDb",allRows),high=collect(data,"highResidualDb",allRows),ultra=collect(data,"ultrasonicResidualDb",allRows);
                    addStats(feature,"overallResidualDb",overall);addStats(feature,"lowResidualDb",low);addStats(feature,"midResidualDb",mid);addStats(feature,"highResidualDb",high);addStats(feature,"ultrasonicResidualDb",ultra);
                    addNumber(feature,"residualLowVsMidDb",mean(low)-mean(mid));addNumber(feature,"residualHighVsMidDb",mean(high)-mean(mid));addNumber(feature,"residualUltrasonicVsHighDb",mean(ultra)-mean(high));
                    addStats(feature,"alignmentLagSamples",collect(data,"lagSamples",allRows));addStats(feature,"alignmentGainDb",collect(data,"gainDb",allRows));
                    addStats(feature,"residualCentroidHz",collect(data,"residualSpectralCentroidHz",allRows));
                    addStats(feature,"residualBandwidthHz",collect(data,"residualEffectiveBandwidthHz",allRows));
                    addStats(feature,"residualEntropyBits",collect(data,"residualSpectralEntropyBits",allRows));
                    addStats(feature,"residualEntropyNormalized",collect(data,"residualSpectralEntropyNormalized",allRows));
                    if(!overall.empty())addNumber(feature,"residualDynamicRangeDb",*std::max_element(overall.begin(),overall.end())-*std::min_element(overall.begin(),overall.end()));
                    useful=!overall.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_interaction_dna"))
                {
                    feature->setProperty("family", "interaction");
                    const int pc=col(data,"productClass"),gc=col(data,"inputGainDb");
                    const auto overall=collect(data,"magDbRelativeToInput",allRows);addStats(feature,"productDb",overall);
                    for(int cls=1;cls<=4;++cls){const auto vals=collect(data,"magDbRelativeToInput",[&](const std::vector<double>& row){return std::llround(valueAt(row,pc))==cls;});addStats(feature,"productClass"+juce::String(cls)+"Db",vals);}
                    std::vector<std::pair<double,double>> growth;for(const double gain:config.inputGainBucketsDb){const double v=meanColumn(data,"magDbRelativeToInput",[&](const std::vector<double>& row){return std::abs(valueAt(row,gc)-gain)<0.01;});if(std::isfinite(v))growth.push_back({gain,v});}
                    if(growth.size()>=2){addNumber(feature,"interactionGrowthDb",growth.back().second-growth.front().second);addNumber(feature,"interactionGrowthDbPerInputDb",(growth.back().second-growth.front().second)/std::max(1.0,growth.back().first-growth.front().first));}
                    if(!overall.empty()){addNumber(feature,"productDensityAboveMinus80Db",static_cast<double>(std::count_if(overall.begin(),overall.end(),[](double v){return v>-80.0;}))/overall.size());addNumber(feature,"worstProductDb",*std::max_element(overall.begin(),overall.end()));}
                    useful=!overall.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_timing_dna"))
                {
                    feature->setProperty("family", "timingPhase");
                    const int rc=col(data,"recordType"),fc=col(data,"frequencyHz"),pc=col(data,"phaseDegrees"),mc=col(data,"magnitudeDb");
                    const auto summary=[&](const std::vector<double>& row){return std::llround(valueAt(row,rc))==2;};
                    addStats(feature,"latencySamples",collect(data,"latencySamples",summary));addStats(feature,"latencyMs",collect(data,"latencyMs",summary));
                    addStats(feature,"preRingingDb",collect(data,"preRingingDb",summary));addStats(feature,"postRingingDb",collect(data,"postRingingDb",summary));
                    addStats(feature,"settlingMs",collect(data,"settlingMs",summary));addStats(feature,"overshootDb",collect(data,"overshootDb",summary));
                    std::vector<std::pair<double,double>> phase;std::vector<double> lowPhase,midPhase,highPhase,lowMag,midMag,highMag;
                    for(const auto& row:data.rows){if(std::llround(valueAt(row,rc))!=1)continue;const double f=valueAt(row,fc),pdeg=valueAt(row,pc),mag=valueAt(row,mc);if(!(f>0&&std::isfinite(pdeg)))continue;phase.push_back({f,pdeg});if(f<200){lowPhase.push_back(pdeg);lowMag.push_back(mag);}else if(f<3000){midPhase.push_back(pdeg);midMag.push_back(mag);}else{highPhase.push_back(pdeg);highMag.push_back(mag);}}
                    addStats(feature,"lowPhaseDegrees",lowPhase);addStats(feature,"midPhaseDegrees",midPhase);addStats(feature,"highPhaseDegrees",highPhase);addStats(feature,"lowMagnitudeDb",lowMag);addStats(feature,"midMagnitudeDb",midMag);addStats(feature,"highMagnitudeDb",highMag);
                    if(phase.size()>=2){std::sort(phase.begin(),phase.end());std::vector<double> groupDelayMs;for(std::size_t i=1;i<phase.size();++i){const double df=phase[i].first-phase[i-1].first;if(df<=0)continue;double dp=phase[i].second-phase[i-1].second;while(dp>180)dp-=360;while(dp<-180)dp+=360;groupDelayMs.push_back(-dp/360.0/df*1000.0);}addStats(feature,"groupDelayMs",groupDelayMs);}
                    useful=true;
                }
                else if (exportItem.filename.containsIgnoreCase("grid_stereo_dna"))
                {
                    feature->setProperty("family", "stereo");
                    addStats(feature,"leftOnlyCrosstalkDb",collect(data,"leftOnlyLeakDb",allRows));
                    addStats(feature,"rightOnlyCrosstalkDb",collect(data,"rightOnlyLeakDb",allRows));
                    addStats(feature,"midChannelMismatchDb",collect(data,"midChannelMismatchDb",allRows));
                    addStats(feature,"sideChannelMismatchDb",collect(data,"sideChannelMismatchDb",allRows));
                    addStats(feature,"midToSideLeakDb",collect(data,"midToSideLeakDb",allRows));
                    addStats(feature,"sideToMidLeakDb",collect(data,"sideToMidLeakDb",allRows));
                    addStats(feature,"midCorrelation",collect(data,"midCorrelation",allRows));
                    addStats(feature,"sideCorrelation",collect(data,"sideCorrelation",allRows));
                    useful=!data.rows.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_summing_dna"))
                {
                    feature->setProperty("family", "summing");
                    addStats(feature,"separateProcessedRms",collect(data,"separateRms",allRows));
                    addStats(feature,"summedProcessedRms",collect(data,"summedRms",allRows));
                    addStats(feature,"interactionResidualRms",collect(data,"interactionResidualRms",allRows));
                    addStats(feature,"interactionResidualDbRelative",collect(data,"interactionResidualDbRelative",allRows));
                    addStats(feature,"summedCrestDb",collect(data,"summedCrestDb",allRows));
                    addStats(feature,"nonAdditivityPercent",collect(data,"nonAdditivityPercent",allRows));
                    useful=!data.rows.empty();
                }
                else if (exportItem.filename.containsIgnoreCase("grid_boundary_dna"))
                {
                    feature->setProperty("family", "boundaries");
                    const int cc=col(data,"componentClass"),tc=col(data,"testId"),gc=col(data,"inputGainDb"),fc=col(data,"frequencyHz");
                    addStats(feature,"dcOffset",collect(data,"dcOffset",allRows));addStats(feature,"transferDb",collect(data,"transferDb",allRows));addStats(feature,"phaseDegrees",collect(data,"phaseDegrees",allRows));addStats(feature,"outputDbRelativeToInputRms",collect(data,"outputDbRelativeToInputRms",allRows));
                    for(int cls=0;cls<=8;++cls){const auto vals=collect(data,"outputDbRelativeToInputRms",[&](const std::vector<double>& row){return std::llround(valueAt(row,cc))==cls;});if(!vals.empty())addStats(feature,"componentClass"+juce::String(cls)+"Db",vals);}
                    std::vector<std::pair<double,double>> worstByGain;for(const double gain:config.inputGainBucketsDb){const auto vals=collect(data,"outputDbRelativeToInputRms",[&](const std::vector<double>& row){return std::abs(valueAt(row,gc)-gain)<0.01;});if(!vals.empty())worstByGain.push_back({gain,*std::max_element(vals.begin(),vals.end())});}
                    if(worstByGain.size()>=2){addNumber(feature,"boundaryGrowthDb",worstByGain.back().second-worstByGain.front().second);double onset=worstByGain.back().first;for(const auto& p:worstByGain)if(p.second>-80.0){onset=p.first;break;}addNumber(feature,"firstBoundaryProductAboveMinus80InputDb",onset);}
                    addNumber(feature,"highestMeasuredFrequencyHz",meanColumn(data,"frequencyHz",[&](const std::vector<double>& row){return valueAt(row,fc)>0;}));
                    useful=true;
                }

                if (useful) analyserFeatures.add(juce::var(feature)); else delete feature;
            }

            auto* featureDna = new juce::DynamicObject();
            featureDna->setProperty("extractorVersion", 7);
            featureDna->setProperty("policy", "compact evidence v2 plus complete static feature audit, operating-point DNA, alias thresholds, stereo, summing, residual entropy, harmonic stability and transfer efficiency extraction");
            featureDna->setProperty("analysers", analyserFeatures);
            root->setProperty("featureDNA", juce::var(featureDna));
        }

        // Phase 2: Behaviour Engine v1. Converts numeric evidence into concise,
        // measurable signal behaviours. No subjective sound adjectives.
        {
            juce::Array<juce::var> behaviours;
            const auto addBehaviour = [&behaviours](const juce::String& id,
                                                    const juce::String& statement,
                                                    double confidence,
                                                    const juce::StringArray& sources,
                                                    const juce::String& detail)
            {
                auto* item = new juce::DynamicObject();
                item->setProperty("id", id);
                item->setProperty("statement", statement);
                item->setProperty("confidence", juce::jlimit(0.0, 1.0, confidence));
                juce::Array<juce::var> sourceArray;
                for (const auto& source : sources) sourceArray.add(source);
                item->setProperty("evidenceSources", sourceArray);
                item->setProperty("detail", detail);
                behaviours.add(juce::var(item));
            };
            const auto meanColumn = [](const MeasurementDataset& data, const juce::String& name,
                                       const std::function<bool(const std::vector<double>&)>& include)
            {
                const int c = findColumn(data, name.toStdString());
                if (c < 0) return std::numeric_limits<double>::quiet_NaN();
                double sum = 0.0; int count = 0;
                for (const auto& row : data.rows)
                {
                    if (static_cast<std::size_t>(c) >= row.size() || !std::isfinite(row[static_cast<std::size_t>(c)])) continue;
                    if (include && !include(row)) continue;
                    sum += row[static_cast<std::size_t>(c)]; ++count;
                }
                return count > 0 ? sum / static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
            };
            const auto allRows = [](const std::vector<double>&) { return true; };

            for (const auto& exportItem : datasets)
            {
                const auto& data = exportItem.dataset;
                if (exportItem.filename.containsIgnoreCase("grid_thd_sine"))
                {
                    const double h2 = meanColumn(data, "h2", allRows);
                    const double h3 = meanColumn(data, "h3", allRows);
                    if (std::isfinite(h2) && std::isfinite(h3) && h2 > 0.0 && h3 > 0.0)
                    {
                        const double differenceDb = amplitudeToDb(h2) - amplitudeToDb(h3);
                        if (differenceDb > 4.0)
                            addBehaviour("builds_even_harmonics", "Builds mainly even-order harmonics",
                                         juce::jlimit(0.60, 0.96, 0.68 + differenceDb / 50.0),
                                         {"THD sine"}, "H2 exceeds H3 by " + juce::String(differenceDb, 1) + " dB.");
                        else if (differenceDb < -4.0)
                            addBehaviour("builds_odd_harmonics", "Builds mainly odd-order harmonics",
                                         juce::jlimit(0.60, 0.96, 0.68 + std::abs(differenceDb) / 50.0),
                                         {"THD sine"}, "H3 exceeds H2 by " + juce::String(std::abs(differenceDb), 1) + " dB.");
                        else
                            addBehaviour("builds_mixed_harmonics", "Builds a mixed even/odd harmonic ladder",
                                         0.70, {"THD sine"}, "H2 and H3 remain within 4 dB.");
                    }
                }
                else if (exportItem.filename.containsIgnoreCase("linear_response_sweep"))
                {
                    const int fc = findColumn(data, "freqHz");
                    const auto band = [&](double lo, double hi)
                    {
                        return meanColumn(data, "magDb", [=](const std::vector<double>& row)
                        {
                            if (fc < 0 || static_cast<std::size_t>(fc) >= row.size()) return false;
                            const double f = row[static_cast<std::size_t>(fc)];
                            return f >= lo && f <= hi;
                        });
                    };
                    const double low = band(40.0, 150.0), mid = band(500.0, 2000.0), high = band(10000.0, 20000.0);
                    if (std::isfinite(low) && std::isfinite(mid) && std::isfinite(high))
                    {
                        const double lowDelta = low - mid, highDelta = high - mid;
                        if (lowDelta > 0.5)
                            addBehaviour("pushes_lows", "Pushes lows relative to the midrange", 0.82,
                                         {"Sweep linear response"}, "Low band is " + juce::String(lowDelta, 2) + " dB above the midrange.");
                        else if (lowDelta < -0.5)
                            addBehaviour("reduces_lows", "Reduces lows relative to the midrange", 0.82,
                                         {"Sweep linear response"}, "Low band is " + juce::String(std::abs(lowDelta), 2) + " dB below the midrange.");
                        if (highDelta < -0.5)
                            addBehaviour("rolls_highs", "Rolls extreme highs relative to the midrange", 0.86,
                                         {"Sweep linear response"}, "High band is " + juce::String(std::abs(highDelta), 2) + " dB below the midrange.");
                        else if (highDelta > 0.5)
                            addBehaviour("pushes_highs", "Pushes extreme highs relative to the midrange", 0.86,
                                         {"Sweep linear response"}, "High band is " + juce::String(highDelta, 2) + " dB above the midrange.");
                        if (std::abs(lowDelta) <= 0.5 && std::abs(highDelta) <= 0.5)
                            addBehaviour("keeps_frequency_balance", "Keeps broad frequency balance essentially flat", 0.80,
                                         {"Sweep linear response"}, "Low and high bands remain within 0.5 dB of the midrange.");
                    }
                }
                else if (exportItem.filename.containsIgnoreCase("grid_residual_dna"))
                {
                    const double residual = meanColumn(data, "residualDbRelativeToOutput", allRows);
                    const double low = meanColumn(data, "lowResidualDb", allRows);
                    const double mid = meanColumn(data, "midResidualDb", allRows);
                    const double high = meanColumn(data, "highResidualDb", allRows);
                    if (std::isfinite(residual))
                    {
                        if (residual < -50.0)
                            addBehaviour("adds_little_fingerprint", "Adds little residual fingerprint", 0.88,
                                         {"Residual DNA"}, "Aligned residual averages " + juce::String(residual, 1) + " dB relative to output.");
                        else if (residual > -25.0)
                            addBehaviour("adds_strong_fingerprint", "Adds a strong residual fingerprint", 0.84,
                                         {"Residual DNA"}, "Aligned residual averages " + juce::String(residual, 1) + " dB relative to output.");
                    }
                    if (std::isfinite(low) && std::isfinite(mid) && std::isfinite(high))
                    {
                        if (low > mid + 3.0 && low > high + 3.0)
                            addBehaviour("adds_low_texture", "Concentrates its added fingerprint in the low band", 0.84,
                                         {"Residual DNA"}, "Low-band residual exceeds mid and high bands by more than 3 dB.");
                        else if (high > low + 3.0 && high > mid + 3.0)
                            addBehaviour("adds_high_texture", "Concentrates its added fingerprint in the high band", 0.84,
                                         {"Residual DNA"}, "High-band residual exceeds low and mid bands by more than 3 dB.");
                    }
                }
                else if (exportItem.filename.containsIgnoreCase("grid_timing_dna"))
                {
                    const int rc = findColumn(data, "recordType");
                    const auto summaryRows = [rc](const std::vector<double>& row)
                    {
                        return rc < 0 || (static_cast<std::size_t>(rc) < row.size() && std::llround(row[static_cast<std::size_t>(rc)]) == 2);
                    };
                    const double settling = meanColumn(data, "settlingMs", summaryRows);
                    const double latency = meanColumn(data, "latencySamples", summaryRows);
                    const double post = meanColumn(data, "postRingingDb", summaryRows);
                    if (std::isfinite(settling) && settling <= 0.10)
                        addBehaviour("settles_quickly", "Settles quickly after an impulse", 0.90,
                                     {"Timing DNA"}, "Mean settling time is " + juce::String(settling, 3) + " ms.");
                    if (std::isfinite(latency) && latency <= 1.0)
                        addBehaviour("keeps_timing_immediate", "Keeps timing essentially immediate", 0.90,
                                     {"Timing DNA"}, "Measured latency averages " + juce::String(latency, 2) + " samples.");
                    if (std::isfinite(post) && post > -35.0)
                        addBehaviour("adds_ringing", "Adds measurable post-impulse ringing", 0.72,
                                     {"Timing DNA"}, "Post-ringing reaches " + juce::String(post, 1) + " dB relative to the main peak.");
                }
                else if (exportItem.filename.containsIgnoreCase("grid_interaction_dna"))
                {
                    const double interaction = meanColumn(data, "magDbRelativeToInput", allRows);
                    if (std::isfinite(interaction) && interaction < -60.0)
                        addBehaviour("keeps_signals_independent", "Keeps simultaneous tones largely independent", 0.86,
                                     {"Interaction DNA"}, "Generated interaction products average " + juce::String(interaction, 1) + " dB relative to input.");
                    else if (std::isfinite(interaction) && interaction > -40.0)
                        addBehaviour("creates_interaction_products", "Creates substantial products between simultaneous tones", 0.84,
                                     {"Interaction DNA"}, "Generated interaction products average " + juce::String(interaction, 1) + " dB relative to input.");
                }
            }

            auto* behaviourDna = new juce::DynamicObject();
            behaviourDna->setProperty("dictionaryVersion", 1);
            behaviourDna->setProperty("languagePolicy", "signal behaviour only; no subjective sound adjectives");
            behaviourDna->setProperty("behaviours", behaviours);
            behaviourDna->setProperty("status", behaviours.isEmpty()
                ? "insufficient derived evidence"
                : "measured behaviours; thresholds require calibration across the reference set");
            root->setProperty("behaviourDNA", juce::var(behaviourDna));
        }

        // Phase 2: compact Parameter DNA. This deliberately describes measured
        // movement first; semantic producer labels are calibrated separately.
        juce::Array<juce::var> parameterDna;
        for (const auto& bucket : config.parameterBuckets)
        {
            std::vector<double> values;
            const auto generatedValues = generatedBucketValues(bucket);
            for (const auto& value : generatedValues)
                values.push_back(static_cast<double>(value));
            if (values.empty())
                continue;

            const double defaultValue = static_cast<double>(bucket.pluginDefaultValue);
            std::size_t baselineIndex = 0;
            for (std::size_t i = 1; i < values.size(); ++i)
                if (std::abs(values[i] - defaultValue) < std::abs(values[baselineIndex] - defaultValue))
                    baselineIndex = i;

            std::vector<double> activity(values.size(), 0.0);
            std::vector<int> activityCounts(values.size(), 0);
            std::map<juce::String, double> sourceActivity;

            const auto isStructuralColumn = [&](const juce::String& name)
            {
                if (name == bucket.paramName || name.equalsIgnoreCase("runId")
                    || name.equalsIgnoreCase("inputGainDb") || name.equalsIgnoreCase("binIndex")
                    || name.equalsIgnoreCase("centreSample") || name.equalsIgnoreCase("count")
                    || name.equalsIgnoreCase("frequencyHz") || name.equalsIgnoreCase("fundamentalHz")
                    || name.equalsIgnoreCase("harmonicOrder") || name.equalsIgnoreCase("recordType")
                    || name.equalsIgnoreCase("testId") || name.equalsIgnoreCase("productHz")
                    || name.equalsIgnoreCase("productClass"))
                    return true;
                for (const auto& other : config.parameterBuckets)
                    if (name == other.paramName)
                        return true;
                return false;
            };

            for (const auto& item : datasets)
            {
                const auto& dataset = item.dataset;
                int parameterColumn = -1;
                for (std::size_t c = 0; c < dataset.columns.size(); ++c)
                    if (juce::String(dataset.columns[c]) == bucket.paramName)
                        parameterColumn = static_cast<int>(c);
                if (parameterColumn < 0)
                    continue;

                double sourceTotal = 0.0;
                int sourceCount = 0;
                for (std::size_t metricColumn = 0; metricColumn < dataset.columns.size(); ++metricColumn)
                {
                    const auto metricName = juce::String(dataset.columns[metricColumn]);
                    if (isStructuralColumn(metricName))
                        continue;

                    double globalMin = std::numeric_limits<double>::infinity();
                    double globalMax = -std::numeric_limits<double>::infinity();
                    for (const auto& row : dataset.rows)
                    {
                        if (metricColumn >= row.size() || !std::isfinite(row[metricColumn]))
                            continue;
                        bool otherParametersAtDefault = true;
                        for (const auto& other : config.parameterBuckets)
                        {
                            if (other.paramName == bucket.paramName)
                                continue;
                            const int otherColumn = findColumn(dataset, other.paramName.toStdString());
                            if (otherColumn >= 0 && static_cast<std::size_t>(otherColumn) < row.size()
                                && std::abs(row[static_cast<std::size_t>(otherColumn)]
                                            - static_cast<double>(other.pluginDefaultValue)) > 1.0e-4)
                            {
                                otherParametersAtDefault = false;
                                break;
                            }
                        }
                        if (!otherParametersAtDefault)
                            continue;
                        globalMin = std::min(globalMin, row[metricColumn]);
                        globalMax = std::max(globalMax, row[metricColumn]);
                    }
                    const double range = globalMax - globalMin;
                    if (!std::isfinite(range) || range < 1.0e-12)
                        continue;

                    std::vector<double> sums(values.size(), 0.0);
                    std::vector<int> counts(values.size(), 0);
                    for (const auto& row : dataset.rows)
                    {
                        if (static_cast<std::size_t>(parameterColumn) >= row.size()
                            || metricColumn >= row.size() || !std::isfinite(row[metricColumn]))
                            continue;
                        bool otherParametersAtDefault = true;
                        for (const auto& other : config.parameterBuckets)
                        {
                            if (other.paramName == bucket.paramName)
                                continue;
                            const int otherColumn = findColumn(dataset, other.paramName.toStdString());
                            if (otherColumn >= 0 && static_cast<std::size_t>(otherColumn) < row.size()
                                && std::abs(row[static_cast<std::size_t>(otherColumn)]
                                            - static_cast<double>(other.pluginDefaultValue)) > 1.0e-4)
                            {
                                otherParametersAtDefault = false;
                                break;
                            }
                        }
                        if (!otherParametersAtDefault)
                            continue;
                        const double parameterValue = row[static_cast<std::size_t>(parameterColumn)];
                        std::size_t nearest = 0;
                        for (std::size_t i = 1; i < values.size(); ++i)
                            if (std::abs(values[i] - parameterValue) < std::abs(values[nearest] - parameterValue))
                                nearest = i;
                        if (std::abs(values[nearest] - parameterValue) > 1.0e-4)
                            continue;
                        sums[nearest] += row[metricColumn];
                        ++counts[nearest];
                    }
                    if (counts[baselineIndex] == 0)
                        continue;
                    const double baselineMean = sums[baselineIndex] / static_cast<double>(counts[baselineIndex]);
                    for (std::size_t i = 0; i < values.size(); ++i)
                    {
                        if (counts[i] == 0)
                            continue;
                        const double mean = sums[i] / static_cast<double>(counts[i]);
                        const double delta = juce::jlimit(0.0, 1.0, std::abs(mean - baselineMean) / range);
                        activity[i] += delta;
                        ++activityCounts[i];
                        sourceTotal += delta;
                        ++sourceCount;
                    }
                }
                if (sourceCount > 0)
                    sourceActivity[item.filename] = sourceTotal / static_cast<double>(sourceCount);
            }

            double maximumActivity = 0.0;
            for (std::size_t i = 0; i < activity.size(); ++i)
            {
                if (activityCounts[i] > 0)
                    activity[i] /= static_cast<double>(activityCounts[i]);
                maximumActivity = std::max(maximumActivity, activity[i]);
            }

            const bool defaultInside = baselineIndex > 0 && baselineIndex + 1 < values.size();
            const double endpointDifference = std::abs(activity.front() - activity.back());
            const bool symmetric = values.size() >= 3 && maximumActivity > 0.01
                                && endpointDifference <= std::max(0.03, maximumActivity * 0.25);
            const bool centreNeutral = defaultInside && maximumActivity > 0.02
                                    && activity[baselineIndex] <= maximumActivity * 0.15;

            bool monotonicFromDefault = !defaultInside;
            if (monotonicFromDefault && values.size() > 2)
            {
                if (baselineIndex == 0)
                    for (std::size_t i = 1; i < activity.size(); ++i)
                        monotonicFromDefault = monotonicFromDefault && activity[i] + 0.02 >= activity[i - 1];
                else
                    for (std::size_t i = activity.size() - 1; i > 0; --i)
                        monotonicFromDefault = monotonicFromDefault && activity[i - 1] + 0.02 >= activity[i];
            }

            bool endLoaded = false;
            if (!defaultInside && activity.size() >= 4)
            {
                std::vector<double> ordered = activity;
                if (baselineIndex != 0)
                    std::reverse(ordered.begin(), ordered.end());
                double earlier = 0.0;
                for (std::size_t i = 1; i + 1 < ordered.size(); ++i)
                    earlier += std::max(0.0, ordered[i] - ordered[i - 1]);
                earlier /= static_cast<double>(ordered.size() - 2);
                const double finalStep = std::max(0.0, ordered.back() - ordered[ordered.size() - 2]);
                endLoaded = finalStep > std::max(0.03, earlier * 1.75);
            }

            auto* dna = new juce::DynamicObject();
            dna->setProperty("name", bucket.paramName);
            dna->setProperty("pluginDefaultValue", defaultValue);
            dna->setProperty("baselineValue", values[baselineIndex]);
            dna->setProperty("baselineIsExactDefault", std::abs(values[baselineIndex] - defaultValue) < 1.0e-5);
            dna->setProperty("strategy", bucket.strategy);
            dna->setProperty("profileCondition", "all other selected parameters held at their plugin defaults");
            dna->setProperty("shape", centreNeutral && symmetric ? "bipolar / centre-neutral"
                                      : monotonicFromDefault ? (endLoaded ? "monotonic / end-loaded" : "monotonic")
                                      : symmetric ? "symmetric" : "complex / non-monotonic");
            dna->setProperty("symmetryDetected", symmetric);
            dna->setProperty("centreNeutralDetected", centreNeutral);
            dna->setProperty("endLoadedDetected", endLoaded);
            dna->setProperty("maximumActivity", maximumActivity);
            double sensitivitySum=0.0, maxStep=0.0; int sensitivityCount=0;
            for(std::size_t i=1;i<activity.size();++i){const double dv=std::abs(values[i]-values[i-1]);if(dv<=0.0)continue;const double step=std::abs(activity[i]-activity[i-1])/dv;sensitivitySum+=step;maxStep=std::max(maxStep,step);++sensitivityCount;}
            dna->setProperty("meanSensitivity",sensitivityCount? sensitivitySum/sensitivityCount : 0.0);
            dna->setProperty("maximumSensitivity",maxStep);
            juce::Array<juce::var> deadValues, usefulValues, criticalValues;
            for(std::size_t i=0;i<values.size();++i){if(activity[i]<0.02)deadValues.add(values[i]);else usefulValues.add(values[i]);if(activity[i]>std::max(0.15,maximumActivity*0.75))criticalValues.add(values[i]);}
            dna->setProperty("deadZoneValues",deadValues);dna->setProperty("usefulRangeValues",usefulValues);dna->setProperty("criticalRangeValues",criticalValues);
            dna->setProperty("usefulRangeCoverage",values.empty()?0.0:static_cast<double>(usefulValues.size())/values.size());
            dna->setProperty("controlClassCandidate", centreNeutral && symmetric ? "setup control"
                                                   : monotonicFromDefault ? "amount control"
                                                   : "character control");

            juce::Array<juce::var> positions;
            for (std::size_t i = 0; i < values.size(); ++i)
            {
                auto* position = new juce::DynamicObject();
                position->setProperty("value", values[i]);
                position->setProperty("isDefault", std::abs(values[i] - defaultValue) < 1.0e-5);
                position->setProperty("deltaActivity", activity[i]);
                positions.add(juce::var(position));
            }
            dna->setProperty("positions", positions);

            std::vector<std::pair<juce::String, double>> rankedSources(sourceActivity.begin(), sourceActivity.end());
            std::sort(rankedSources.begin(), rankedSources.end(), [](const auto& a, const auto& b)
            {
                return a.second > b.second;
            });
            juce::Array<juce::var> dominantSources;
            for (std::size_t i = 0; i < std::min<std::size_t>(3, rankedSources.size()); ++i)
            {
                auto* source = new juce::DynamicObject();
                source->setProperty("source", rankedSources[i].first);
                source->setProperty("relativeMovement", rankedSources[i].second);
                dominantSources.add(juce::var(source));
            }
            dna->setProperty("dominantEvidenceSources", dominantSources);

            // Phase 3.1: universal DSP vocabulary. Preserve the manufacturer's
            // parameter name, then describe the measured dimensions it changes.
            struct InfluenceAxis
            {
                juce::String id;
                juce::String name;
                juce::String category;
                juce::String generalisation;
                juce::String explanation;
                double score = 0.0;
                juce::StringArray evidence;
            };

            std::vector<InfluenceAxis> axes {
                { "transfer_symmetry", "Transfer symmetry", "Transfer",
                  "Changes how the nonlinear stage treats positive and negative waveform excursions",
                  "Changes waveform symmetry and therefore the balance of even and odd harmonics.", 0.0, {} },
                { "transfer_shape", "Transfer shape", "Transfer",
                  "Changes the shape or strength of the nonlinear transfer curve",
                  "Changes how progressively peaks are bent, compressed, expanded or saturated.", 0.0, {} },
                { "harmonic_balance", "Harmonic balance", "Harmonics",
                  "Changes the balance and distribution of generated harmonics",
                  "Changes which harmonic orders dominate and how the harmonic ladder develops.", 0.0, {} },
                { "spectral_tilt", "Spectral tilt", "Frequency",
                  "Changes the balance between low, mid and high frequencies",
                  "Changes broad frequency weighting, such as low lift, mid emphasis or high-frequency roll-off.", 0.0, {} },
                { "crest_behaviour", "Crest behaviour", "Dynamics",
                  "Changes the relationship between average level and peaks",
                  "Changes peak density or crest factor without necessarily changing average level by the same amount.", 0.0, {} },
                { "residual_fingerprint", "Residual fingerprint", "Residual",
                  "Changes the added signal left after the aligned dry signal is removed",
                  "Changes the level, spectrum, density or tonality of the processor's added fingerprint.", 0.0, {} },
                { "timing_response", "Timing response", "Timing",
                  "Changes how energy arrives and settles in time",
                  "Changes onset, ringing, overshoot, impulse width or settling behaviour.", 0.0, {} },
                { "phase_response", "Phase response", "Phase",
                  "Changes frequency-dependent timing alignment",
                  "Changes phase rotation or group delay across the spectrum.", 0.0, {} },
                { "signal_interaction", "Signal interaction", "Interaction",
                  "Changes how simultaneous tones interact",
                  "Changes intermodulation products and how complex signals generate new components.", 0.0, {} },
                { "boundary_behaviour", "Boundary behaviour", "Limits",
                  "Changes behaviour near DC, infrasonic, ultrasonic or level limits",
                  "Changes stability, foldback, DC offset or other edge-of-range behaviour.", 0.0, {} }
            };

            const auto addAxisEvidence = [&axes](const juce::String& id, double amount, const juce::String& source)
            {
                for (auto& axis : axes)
                {
                    if (axis.id != id) continue;
                    axis.score += juce::jlimit(0.0, 1.0, amount);
                    if (!axis.evidence.contains(source, true)) axis.evidence.add(source);
                    return;
                }
            };

            for (const auto& [sourceName, movement] : sourceActivity)
            {
                const auto source = sourceName.toLowerCase();
                if (source.contains("transfer"))
                {
                    addAxisEvidence("transfer_shape", movement, sourceName);
                    addAxisEvidence("transfer_symmetry", movement * 0.85, sourceName);
                }
                if (source.contains("thd") || source.contains("harmonic"))
                {
                    addAxisEvidence("harmonic_balance", movement, sourceName);
                    addAxisEvidence("transfer_symmetry", movement * 0.55, sourceName);
                }
                if (source.contains("linear") || source.contains("sweep"))
                    addAxisEvidence("spectral_tilt", movement, sourceName);
                if (source.contains("rms") || source.contains("peak"))
                    addAxisEvidence("crest_behaviour", movement, sourceName);
                if (source.contains("residual"))
                    addAxisEvidence("residual_fingerprint", movement, sourceName);
                if (source.contains("timing"))
                {
                    addAxisEvidence("timing_response", movement, sourceName);
                    addAxisEvidence("phase_response", movement * 0.65, sourceName);
                }
                if (source.contains("interaction"))
                    addAxisEvidence("signal_interaction", movement, sourceName);
                if (source.contains("boundary"))
                    addAxisEvidence("boundary_behaviour", movement, sourceName);
            }

            double topAxis = 0.0;
            for (const auto& axis : axes) topAxis = std::max(topAxis, axis.score);
            std::stable_sort(axes.begin(), axes.end(), [](const InfluenceAxis& a, const InfluenceAxis& b)
            {
                return a.score > b.score;
            });

            juce::Array<juce::var> influenceMatrix;
            juce::Array<juce::var> whatItChanges;
            juce::Array<juce::var> littleEffectOn;
            for (const auto& axis : axes)
            {
                auto* influence = new juce::DynamicObject();
                influence->setProperty("id", axis.id);
                influence->setProperty("name", axis.name);
                influence->setProperty("category", axis.category);
                const double normalised = topAxis > 1.0e-12 ? juce::jlimit(0.0, 1.0, axis.score / topAxis) : 0.0;
                influence->setProperty("influence", normalised);
                influence->setProperty("influencePercent", juce::roundToInt(normalised * 100.0));
                influence->setProperty("generalisation", axis.generalisation);
                influence->setProperty("engineeringEffect", axis.explanation);
                juce::Array<juce::var> axisEvidence;
                for (const auto& e : axis.evidence) axisEvidence.add(e);
                influence->setProperty("evidenceSources", axisEvidence);
                influenceMatrix.add(juce::var(influence));

                if (normalised >= 0.45 && whatItChanges.size() < 3)
                    whatItChanges.add(axis.name);
                if (normalised <= 0.12 && littleEffectOn.size() < 4)
                    littleEffectOn.add(axis.name);
            }

            dna->setProperty("whatItChanges", whatItChanges);
            dna->setProperty("measuredInfluence", influenceMatrix);
            dna->setProperty("littleEffectOn", littleEffectOn);

            if (!axes.empty() && axes.front().score > 1.0e-12)
            {
                dna->setProperty("primaryMeasuredConcept", axes.front().name);
                dna->setProperty("generalisation", axes.front().generalisation);
                dna->setProperty("engineeringEffect", axes.front().explanation);
                dna->setProperty("conceptConfidence", juce::jlimit(0.0, 1.0,
                    0.45 + 0.12 * static_cast<double>(axes.front().evidence.size())
                    + 0.25 * (topAxis > 0.0 ? axes.front().score / topAxis : 0.0)));
            }

            dna->setProperty("interpretationStatus",
                "universal measured vocabulary; manufacturer parameter name preserved");
            parameterDna.add(juce::var(dna));
        }
        root->setProperty("parameterDNA", parameterDna);

        // Phase 2: Coupling DNA. Only two selected parameters are interpreted as
        // a deliberate pair. Individual profiles above remain conditioned on the
        // other control's real plugin default.
        if (config.parameterBuckets.size() == 2)
        {
            const auto& first = config.parameterBuckets[0];
            const auto& second = config.parameterBuckets[1];
            std::vector<double> firstValues;
            std::vector<double> secondValues;
            for (const auto& value : generatedBucketValues(first)) firstValues.push_back(static_cast<double>(value));
            for (const auto& value : generatedBucketValues(second)) secondValues.push_back(static_cast<double>(value));

            const auto nearestIndex = [](const std::vector<double>& values, double target)
            {
                std::size_t index = 0;
                for (std::size_t i = 1; i < values.size(); ++i)
                    if (std::abs(values[i] - target) < std::abs(values[index] - target)) index = i;
                return index;
            };

            if (!firstValues.empty() && !secondValues.empty())
            {
                const auto firstDefaultIndex = nearestIndex(firstValues, first.pluginDefaultValue);
                const auto secondDefaultIndex = nearestIndex(secondValues, second.pluginDefaultValue);
                std::vector<std::vector<double>> coupling(firstValues.size(), std::vector<double>(secondValues.size(), 0.0));
                std::vector<std::vector<int>> couplingCounts(firstValues.size(), std::vector<int>(secondValues.size(), 0));
                std::map<juce::String, double> sourceCoupling;

                for (const auto& item : datasets)
                {
                    const auto& dataset = item.dataset;
                    const int firstColumn = findColumn(dataset, first.paramName.toStdString());
                    const int secondColumn = findColumn(dataset, second.paramName.toStdString());
                    if (firstColumn < 0 || secondColumn < 0)
                        continue;

                    double sourceTotal = 0.0;
                    int sourceCount = 0;
                    for (std::size_t metricColumn = 0; metricColumn < dataset.columns.size(); ++metricColumn)
                    {
                        const auto metricName = juce::String(dataset.columns[metricColumn]);
                        if (metricName == first.paramName || metricName == second.paramName
                            || metricName.equalsIgnoreCase("runId") || metricName.equalsIgnoreCase("inputGainDb")
                            || metricName.equalsIgnoreCase("binIndex") || metricName.equalsIgnoreCase("centreSample")
                            || metricName.equalsIgnoreCase("count") || metricName.equalsIgnoreCase("frequencyHz")
                            || metricName.equalsIgnoreCase("fundamentalHz") || metricName.equalsIgnoreCase("harmonicOrder")
                            || metricName.equalsIgnoreCase("recordType") || metricName.equalsIgnoreCase("testId")
                            || metricName.equalsIgnoreCase("productHz") || metricName.equalsIgnoreCase("productClass"))
                            continue;

                        double globalMin = std::numeric_limits<double>::infinity();
                        double globalMax = -std::numeric_limits<double>::infinity();
                        std::vector<std::vector<double>> sums(firstValues.size(), std::vector<double>(secondValues.size(), 0.0));
                        std::vector<std::vector<int>> counts(firstValues.size(), std::vector<int>(secondValues.size(), 0));

                        for (const auto& row : dataset.rows)
                        {
                            if (static_cast<std::size_t>(firstColumn) >= row.size()
                                || static_cast<std::size_t>(secondColumn) >= row.size()
                                || metricColumn >= row.size() || !std::isfinite(row[metricColumn]))
                                continue;
                            const auto i = nearestIndex(firstValues, row[static_cast<std::size_t>(firstColumn)]);
                            const auto j = nearestIndex(secondValues, row[static_cast<std::size_t>(secondColumn)]);
                            if (std::abs(firstValues[i] - row[static_cast<std::size_t>(firstColumn)]) > 1.0e-4
                                || std::abs(secondValues[j] - row[static_cast<std::size_t>(secondColumn)]) > 1.0e-4)
                                continue;
                            sums[i][j] += row[metricColumn];
                            ++counts[i][j];
                            globalMin = std::min(globalMin, row[metricColumn]);
                            globalMax = std::max(globalMax, row[metricColumn]);
                        }

                        const double range = globalMax - globalMin;
                        if (!std::isfinite(range) || range < 1.0e-12
                            || counts[firstDefaultIndex][secondDefaultIndex] == 0)
                            continue;
                        const double baseline = sums[firstDefaultIndex][secondDefaultIndex]
                                              / static_cast<double>(counts[firstDefaultIndex][secondDefaultIndex]);

                        for (std::size_t i = 0; i < firstValues.size(); ++i)
                        {
                            if (counts[i][secondDefaultIndex] == 0) continue;
                            const double firstOnly = sums[i][secondDefaultIndex]
                                                   / static_cast<double>(counts[i][secondDefaultIndex]);
                            for (std::size_t j = 0; j < secondValues.size(); ++j)
                            {
                                if (counts[firstDefaultIndex][j] == 0 || counts[i][j] == 0) continue;
                                const double secondOnly = sums[firstDefaultIndex][j]
                                                        / static_cast<double>(counts[firstDefaultIndex][j]);
                                const double joint = sums[i][j] / static_cast<double>(counts[i][j]);
                                const double expectedAdditive = firstOnly + secondOnly - baseline;
                                const double synergy = juce::jlimit(0.0, 1.0, std::abs(joint - expectedAdditive) / range);
                                coupling[i][j] += synergy;
                                ++couplingCounts[i][j];
                                sourceTotal += synergy;
                                ++sourceCount;
                            }
                        }
                    }
                    if (sourceCount > 0)
                        sourceCoupling[item.filename] = sourceTotal / static_cast<double>(sourceCount);
                }

                double overall = 0.0;
                int overallCount = 0;
                double maximum = 0.0;
                std::size_t maximumI = firstDefaultIndex;
                std::size_t maximumJ = secondDefaultIndex;
                juce::Array<juce::var> cells;
                for (std::size_t i = 0; i < firstValues.size(); ++i)
                {
                    for (std::size_t j = 0; j < secondValues.size(); ++j)
                    {
                        if (couplingCounts[i][j] > 0)
                            coupling[i][j] /= static_cast<double>(couplingCounts[i][j]);
                        if (i == firstDefaultIndex && j == secondDefaultIndex)
                            continue;
                        overall += coupling[i][j];
                        ++overallCount;
                        if (coupling[i][j] > maximum)
                        {
                            maximum = coupling[i][j];
                            maximumI = i;
                            maximumJ = j;
                        }
                        auto* cell = new juce::DynamicObject();
                        cell->setProperty(first.paramName, firstValues[i]);
                        cell->setProperty(second.paramName, secondValues[j]);
                        cell->setProperty("nonAdditiveCoupling", coupling[i][j]);
                        cells.add(juce::var(cell));
                    }
                }
                if (overallCount > 0) overall /= static_cast<double>(overallCount);

                auto* couplingDna = new juce::DynamicObject();
                couplingDna->setProperty("parameters", juce::StringArray{first.paramName, second.paramName});
                couplingDna->setProperty("reference", "both parameters at their actual plugin defaults");
                couplingDna->setProperty("overallCoupling", overall);
                couplingDna->setProperty("maximumCoupling", maximum);
                couplingDna->setProperty("strength", maximum >= 0.12 ? "strong"
                                                  : maximum >= 0.04 ? "moderate" : "weak");
                auto* strongest = new juce::DynamicObject();
                strongest->setProperty(first.paramName, firstValues[maximumI]);
                strongest->setProperty(second.paramName, secondValues[maximumJ]);
                couplingDna->setProperty("strongestInteractionAt", juce::var(strongest));
                couplingDna->setProperty("grid", cells);

                std::vector<std::pair<juce::String, double>> rankedSources(sourceCoupling.begin(), sourceCoupling.end());
                std::sort(rankedSources.begin(), rankedSources.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
                juce::Array<juce::var> dominantSources;
                for (std::size_t i = 0; i < std::min<std::size_t>(3, rankedSources.size()); ++i)
                {
                    auto* source = new juce::DynamicObject();
                    source->setProperty("source", rankedSources[i].first);
                    source->setProperty("relativeCoupling", rankedSources[i].second);
                    dominantSources.add(juce::var(source));
                }
                couplingDna->setProperty("dominantEvidenceSources", dominantSources);
                couplingDna->setProperty("interpretationStatus", "non-additive interaction only; musical meaning pending calibration");
                root->setProperty("couplingDNA", juce::var(couplingDna));
            }
        }

        root->setProperty("rawDataIncluded", false);
        root->setProperty("rawDataPackagedAlongside", exportingFullPackage);
        if (exportingFullPackage)
        {
            root->setProperty("rawDataLocation", ".");
            root->setProperty("note", "Compact evidence JSON. Complete analyser CSV rows are stored alongside this file in the same full-evidence package.");
        }
        else
        {
            root->setProperty("note", "Compact Evidence v2. Per-run raw columns and redundant dataset bookkeeping are omitted; operatingPointDNA retains the derived fingerprint for every tested state. Use Export Raw Data for every analyser CSV and lossless evidence.");
        }

        const auto json = juce::JSON::toString(juce::var(root), true);
        auto baseName = pluginFile.getFileNameWithoutExtension();
        if (baseName.isEmpty())
            baseName = "PluginDNA";
        if (!config.parameterBuckets.empty())
        {
            juce::StringArray parameterNames;
            for (const auto& bucket : config.parameterBuckets)
                parameterNames.add(bucket.paramName);
            baseName += "_" + parameterNames.joinIntoString("_");
        }
        baseName = baseName.retainCharacters("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_ ")
                           .trim().replaceCharacter(' ', '_');
        const auto filename = baseName + "_evidence_"
                            + juce::Time::getCurrentTime().formatted("%Y%m%d_%H%M%S") + ".json";
        const auto file = outDir.getChildFile(filename);
        bool success = file.replaceWithText(json);
        const auto bytes = file.existsAsFile() ? file.getSize() : 0;

        if (success && exportingFullPackage)
        {
            auto* manifest = new juce::DynamicObject();
            manifest->setProperty("schema", "PluginDNA Full Evidence Manifest");
            manifest->setProperty("schemaVersion", 24);
            manifest->setProperty("createdUtc", juce::Time::getCurrentTime().toISO8601(true));
            manifest->setProperty("evidenceFile", filename);
            manifest->setProperty("plugin", juce::File(pluginPath).getFileNameWithoutExtension());
            if (serialPluginPath.isNotEmpty())
                manifest->setProperty("serialPlugin", juce::File(serialPluginPath).getFileNameWithoutExtension());
            juce::Array<juce::var> manifestChain;
            manifestChain.add(juce::File(pluginPath).getFileNameWithoutExtension());
            if (serialPluginPath.isNotEmpty()) manifestChain.add(juce::File(serialPluginPath).getFileNameWithoutExtension());
            for (const auto& path : additionalPluginPaths) manifestChain.add(juce::File(path).getFileNameWithoutExtension());
            manifest->setProperty("pluginChain", manifestChain);
            manifest->setProperty("sampleRate", config.sampleRate);
            manifest->setProperty("blockSize", config.blockSize);
            manifest->setProperty("durationSeconds", config.seconds);

            juce::Array<juce::var> analyserNames;
            for (const auto& analyser : config.analyzers)
                analyserNames.add(analyser);
            manifest->setProperty("analysers", analyserNames);

            juce::Array<juce::var> files;
            for (const auto& item : fullPackageFiles)
                files.add(item);
            files.add(filename);
            files.add("manifest.json");
            manifest->setProperty("files", files);

            const auto manifestJson = juce::JSON::toString(juce::var(manifest), true);
            success = outDir.getChildFile("manifest.json").replaceWithText(manifestJson);
        }

        juce::MessageManager::callAsync([this, success, filename, bytes]() {
            runMeasurementButton.setEnabled(true);
            exportDataButton.setEnabled(true);
            exportEvidenceButton.setEnabled(true);
            if (success)
                progressLabel.setText("Exported evidence pack: " + filename + " ("
                                      + juce::File::descriptionOfSizeInBytes(bytes) + ")",
                                      juce::dontSendNotification);
            else
                showError("Failed to export evidence pack");
        });
    }).detach();
}

void MainComponent::exportMeasurementData()
{
    const auto outputPath = outputPathEditor.getText().trim();
    if (outputPath.isEmpty())
    {
        showError("Please specify an output path");
        return;
    }

    juce::File outDir(outputPath);
    if (!outDir.exists() && !outDir.createDirectory())
    {
        showError("Could not create output directory");
        return;
    }

    exportDataButton.setEnabled(false);
    runMeasurementButton.setEnabled(false);
    progressLabel.setText("Building full evidence package...", juce::dontSendNotification);

    const auto pluginPath = pluginPathEditor.getText().trim();
    auto baseName = juce::File(pluginPath).getFileNameWithoutExtension();
    if (baseName.isEmpty())
        baseName = "PluginDNA";
    baseName = baseName.retainCharacters("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_ ")
                       .trim().replaceCharacter(' ', '_');
    const auto packageName = baseName + "_full_evidence_"
                           + juce::Time::getCurrentTime().formatted("%Y%m%d_%H%M%S");
    const auto packageDir = outDir.getChildFile(packageName);
    if (!packageDir.createDirectory())
    {
        runMeasurementButton.setEnabled(true);
        exportDataButton.setEnabled(true);
        showError("Could not create full evidence package directory");
        return;
    }

    std::thread([this, packageDir]() {
        bool success = true;
        juce::String failedFile;
        int exportedCount = 0;
        int totalCount = 0;

        {
            std::lock_guard<std::mutex> lock(pendingExportMutex);
            totalCount = static_cast<int>(pendingExportDatasets.size());

            for (const auto& item : pendingExportDatasets)
            {
                const auto file = packageDir.getChildFile(item.filename);
                std::ofstream out(file.getFullPathName().toStdString(),
                                  std::ios::out | std::ios::trunc);

                if (!out.is_open())
                {
                    success = false;
                    failedFile = item.filename;
                    break;
                }

                out << std::setprecision(10);
                for (std::size_t column = 0; column < item.dataset.columns.size(); ++column)
                {
                    if (column > 0)
                        out << ',';
                    out << item.dataset.columns[column];
                }
                out << '\n';

                for (const auto& row : item.dataset.rows)
                {
                    for (std::size_t column = 0; column < row.size(); ++column)
                    {
                        if (column > 0)
                            out << ',';
                        out << row[column];
                    }
                    out << '\n';
                }

                if (!out.good())
                {
                    success = false;
                    failedFile = item.filename;
                    break;
                }

                ++exportedCount;
                juce::MessageManager::callAsync([this, exportedCount, totalCount]() {
                    progressLabel.setText(
                        "Exporting data " + juce::String(exportedCount) + " / " + juce::String(totalCount),
                        juce::dontSendNotification);
                });
            }
        }

        juce::MessageManager::callAsync([this, success, failedFile, exportedCount, packageDir]() {
            if (success)
            {
                progressLabel.setText(
                    "Exported " + juce::String(exportedCount) + " raw data files; building compact evidence...",
                    juce::dontSendNotification);
                exportEvidencePack(packageDir);
            }
            else
            {
                runMeasurementButton.setEnabled(true);
                exportDataButton.setEnabled(true);
                exportEvidenceButton.setEnabled(hasPendingEvidencePack);
                showError("Failed to export " + failedFile);
            }
        });
    }).detach();
}

void MainComponent::clearResultsSummary()
{
    copySummaryButton.setEnabled(false);
    exportDataButton.setEnabled(false);
    exportEvidenceButton.setEnabled(false);
    identityRoleLabel.setText("PRIMARY ROLE  waiting for complete measurement", juce::dontSendNotification);
    identityWhyLabel.setText("WHY  waiting for evidence", juce::dontSendNotification);
    identityWatchLabel.setText("WATCH OUT  waiting for evidence", juce::dontSendNotification);
    levelSummaryLabel.setText("LEVEL MATCH  waiting for measurement", juce::dontSendNotification);
    peakSummaryLabel.setText("TRANSIENT IMPACT  waiting for measurement", juce::dontSendNotification);
    dynamicsSummaryLabel.setText("DYNAMIC RANGE  waiting for measurement", juce::dontSendNotification);
    harmonicSummaryLabel.setText("HARMONIC DNA  waiting for sine measurement", juce::dontSendNotification);
    harmonicBalanceLabel.setText("BALANCE  waiting for harmonic measurement", juce::dontSendNotification);
    harmonicGrowthLabel.setText("WHEN PUSHED  waiting for harmonic measurement", juce::dontSendNotification);
    toneBassLabel.setText("TONE SHAPE  waiting for linear response measurement", juce::dontSendNotification);
    toneMidLabel.setText("LOW END  waiting for linear response measurement", juce::dontSendNotification);
    toneTrebleLabel.setText("HIGH END  waiting for linear response measurement", juce::dontSendNotification);
    toneLargestLabel.setText("STRONGEST FEATURE  waiting for linear response measurement", juce::dontSendNotification);
    behaviourSummaryLabel.setText("BEHAVIOUR  waiting for matched input levels", juce::dontSendNotification);
    behaviourChangesLabel.setText("CHANGES WITH DRIVE  waiting for measurement", juce::dontSendNotification);
    operatingRangeLabel.setText("OPERATING RANGE  waiting for stress levels", juce::dontSendNotification);
    characterStartsLabel.setText("CHARACTER STARTS  waiting for measurement", juce::dontSendNotification);
    underStressLabel.setText("UNDER STRESS  waiting for measurement", juce::dontSendNotification);
    nonlinearityLabel.setText("NONLINEARITY  waiting for transfer curve", juce::dontSendNotification);
    curveBehaviourLabel.setText("CURVE BEHAVIOUR  waiting for measurement", juce::dontSendNotification);
    waveformShapeLabel.setText("WAVE SHAPE  waiting for raw sine measurement", juce::dontSendNotification);
    waveformKneeLabel.setText("WHAT IT DOES  waiting for measurement", juce::dontSendNotification);
    waveformStyleLabel.setText("SIMILAR TO  waiting for measurement", juce::dontSendNotification);
    temporalResponseLabel.setText("TEMPORAL RESPONSE  waiting for raw sine onset", juce::dontSendNotification);
    temporalAttackLabel.setText("ATTACK HANDLING  waiting for measurement", juce::dontSendNotification);
    temporalAfterLabel.setText("AFTER-EFFECTS  waiting for measurement", juce::dontSendNotification);
}