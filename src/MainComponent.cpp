#include "MainComponent.h"
#include "RmsPeakAnalyzer.h"
#include "ThdAnalyzer.h"
#include "LinearResponseAnalyzer.h"
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

MainComponent::MainComponent()
    : pluginPathLabel("Plugin Path:", "Plugin Path:"), pluginInfoLabel("", "No plugin loaded"),
      parametersLabel("Parameters:", "Parameters:"), outputPathLabel("Output Path:", "Output Path:"),
      progressLabel("", "Ready"), progressBar(progress) {
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

    addAndMakeVisible(pluginInfoLabel);

    // Parameter list
    addAndMakeVisible(parametersLabel);
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

    // Progress
    addAndMakeVisible(progressLabel);
    addAndMakeVisible(progressBar);
    progressBar.setPercentageDisplay(false);

    setSize(1200, 1200);
}

MainComponent::~MainComponent() {}

void MainComponent::paint(juce::Graphics& g) {
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));
}

void MainComponent::resized() {
    auto bounds = getLocalBounds().reduced(10);

    // Plugin path section (top)
    auto pluginSection = bounds.removeFromTop(80);
    pluginPathLabel.setBounds(pluginSection.removeFromTop(25));
    auto pluginRow = pluginSection.removeFromTop(30);
    pluginPathEditor.setBounds(pluginRow.removeFromLeft(600));
    pluginRow.removeFromLeft(10);
    browseButton.setBounds(pluginRow.removeFromLeft(100));
    pluginRow.removeFromLeft(10);
    loadPluginButton.setBounds(pluginRow.removeFromLeft(120));
    pluginInfoLabel.setBounds(pluginSection);
    bounds.removeFromTop(10);

    // Main content area (side by side)
    auto contentArea = bounds.removeFromTop(410);

    // Left: Parameter list and config
    auto leftPanel = contentArea.removeFromLeft(500);
    parametersLabel.setBounds(leftPanel.removeFromTop(25));
    auto buttonRow = leftPanel.removeFromTop(30);
    selectAllButton.setBounds(buttonRow.removeFromLeft(100));
    buttonRow.removeFromLeft(10);
    deselectAllButton.setBounds(buttonRow.removeFromLeft(100));
    parameterListBox.setBounds(leftPanel.removeFromTop(150));
    leftPanel.removeFromTop(10);
    parameterConfigViewport.setBounds(leftPanel);

    // Right: Measurement config
    measurementConfig->setBounds(contentArea);

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
    if (button == &browseButton) {
        auto chooser = std::make_shared<juce::FileChooser>("Select VST3 Plugin", juce::File(), "*.vst3");
        auto chooserFlags = juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles;
        chooser->launchAsync(chooserFlags, [this, chooser](const juce::FileChooser& fc) {
            if (fc.getResults().size() > 0) {
                pluginPathEditor.setText(fc.getResult().getFullPathName(), juce::dontSendNotification);
            }
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
    } else if (button == &runMeasurementButton) {
        runMeasurement();
    } else if (button == &copySummaryButton) {
        juce::SystemClipboard::copyTextToClipboard(buildSummaryText());
        progressLabel.setText("Summary copied to clipboard", juce::dontSendNotification);
    }
}

void MainComponent::textEditorTextChanged(juce::TextEditor& editor) {
    // Handle text changes if needed
}

int MainComponent::getNumRows() {
    return (int)availableParameters.size();
}

void MainComponent::paintListBoxItem(int rowNumber, juce::Graphics& g, int width, int height, bool rowIsSelected) {
    if (rowNumber < 0 || rowNumber >= (int)availableParameters.size())
        return;

    // Background - use subtle alternating colors, ignore rowIsSelected since we disabled selection
    g.fillAll(rowNumber % 2 == 0 ? juce::Colours::white : juce::Colours::lightgrey.withAlpha(0.3f));

    // Draw checkbox
    const int checkboxSize = 18;
    const int checkboxX = 5;
    const int checkboxY = (height - checkboxSize) / 2;
    juce::Rectangle<float> checkboxBounds((float)checkboxX, (float)checkboxY, (float)checkboxSize, (float)checkboxSize);

    bool isChecked = rowNumber < (int)selectedParameters.size() && selectedParameters[rowNumber];

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
    g.drawText(availableParameters[rowNumber], checkboxX + checkboxSize + 10, 0, width - checkboxX - checkboxSize - 10,
               height, juce::Justification::centredLeft);
}

void MainComponent::listBoxItemClicked(int row, const juce::MouseEvent& e) {
    // Ensure selectedParameters is the right size
    if (selectedParameters.size() != availableParameters.size()) {
        selectedParameters.resize(availableParameters.size(), false);
    }

    // Only toggle if clicking within the checkbox area or the row itself
    // This prevents accidental toggles when clicking elsewhere
    const int checkboxArea = 30; // Approximate checkbox area width

    // Toggle if clicking in the left part of the row (where checkbox is) or anywhere on the row
    if (row >= 0 && row < (int)selectedParameters.size() && row < (int)availableParameters.size()) {
        // Always toggle on row click - the checkbox visual will update
        selectedParameters[row] = !selectedParameters[row];
        parameterListBox.updateContent();
        parameterListBox.repaintRow(row); // Repaint the specific row that changed
        parameterListBox.repaint();       // Also repaint the whole component
        updateParameterList();
    }
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
    scanPluginParameters();
    runMeasurementButton.setEnabled(true);
    progressLabel.setText("Plugin loaded successfully", juce::dontSendNotification);
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
    parameterListBox.updateContent();
    updateParameterList();
}

void MainComponent::updateParameterList() {
    // Remove old config components
    parameterConfigComponents.clear();
    parameterConfigContainer.removeAllChildren();

    // Create config components for selected parameters
    int y = 10;
    for (size_t i = 0; i < availableParameters.size(); ++i) {
        if (i < selectedParameters.size() && selectedParameters[i]) {
            auto* comp = new ParameterConfigComponent(availableParameters[i]);
            comp->setBounds(10, y, 480, 150);
            parameterConfigContainer.addAndMakeVisible(comp);
            parameterConfigComponents.push_back(std::unique_ptr<ParameterConfigComponent>(comp));
            y += 160;
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

    // Count selected parameters
    int selectedCount = 0;
    for (bool selected : selectedParameters) {
        if (selected)
            selectedCount++;
    }

    if (selectedCount == 0) {
        showError("Please select at least one parameter");
        return;
    }

    juce::String outputPath = outputPathEditor.getText();
    if (outputPath.isEmpty()) {
        showError("Please specify an output path");
        return;
    }

    juce::File outDir(outputPath);
    if (!outDir.exists()) {
        outDir.createDirectory();
    }

    if (!outDir.isDirectory()) {
        showError("Output path is not a directory");
        return;
    }

    // Build config from UI
    Config config = buildConfigFromUI();
    config.pluginPath = pluginPathEditor.getText();

    // Run measurement in background thread
    runMeasurementButton.setEnabled(false);
    progressLabel.setText("Running measurement...", juce::dontSendNotification);
    clearResultsSummary();

    std::thread([this, config, outDir]() {
        try {

            // Create a separate plugin instance for the measurement thread
            // This is necessary because JUCE plugins should not be accessed from multiple threads
            juce::File pluginFile(config.pluginPath);
            juce::String errorMessage;
            auto measurementPlugin = loadPluginInstance(pluginFile, config.sampleRate, config.blockSize, errorMessage);

            if (measurementPlugin == nullptr) {
                std::cerr << "[Measurement] Failed to create plugin instance: " << errorMessage << std::endl;
                juce::MessageManager::callAsync([this, errorMessage]() {
                    showError("Failed to create plugin instance for measurement: " + errorMessage);
                    runMeasurementButton.setEnabled(true);
                });
                return;
            }

            // Build parameter name list
            std::vector<juce::String> paramNames;
            for (size_t i = 0; i < availableParameters.size(); ++i) {
                if (i < selectedParameters.size() && selectedParameters[i]) {
                    paramNames.push_back(availableParameters[i]);
                }
            }

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

            // Create analyzers
            auto analyzers = createAnalyzers(config, outDir, paramNames);

            // Run measurements
            int64_t totalSamples = (int64_t)(config.seconds * config.sampleRate);
            juce::MessageManager::callAsync([this, runs]() {
                progressLabel.setText("Running " + juce::String(runs.size()) + " measurements...",
                                      juce::dontSendNotification);
            });

            // Pass progress callback to update UI with time estimate
            auto startTime = std::chrono::steady_clock::now();
            runMeasurementGrid(
                *measurementPlugin, config.sampleRate, config.blockSize, totalSamples, runs, analyzers, config, outDir,
                [this, runs, startTime](int runIndex) {
                    double progress = (double)(runIndex + 1) / (double)runs.size();
                    auto currentTime = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

                    juce::String statusText = "Run " + juce::String(runIndex + 1) + " / " + juce::String(runs.size());

                    // Estimate remaining time
                    if (runIndex > 0 && elapsed > 0) {
                        double runsPerSecond = (double)(runIndex + 1) / (double)elapsed;
                        int remainingRuns = runs.size() - (runIndex + 1);
                        int estimatedSecondsRemaining = (int)(remainingRuns / runsPerSecond);

                        int hours = estimatedSecondsRemaining / 3600;
                        int minutes = (estimatedSecondsRemaining % 3600) / 60;
                        int seconds = estimatedSecondsRemaining % 60;

                        if (hours > 0) {
                            statusText += " (~" + juce::String(hours) + "h " + juce::String(minutes) + "m remaining)";
                        } else if (minutes > 0) {
                            statusText += " (~" + juce::String(minutes) + "m " + juce::String(seconds) + "s remaining)";
                        } else {
                            statusText += " (~" + juce::String(seconds) + "s remaining)";
                        }
                    }

                    juce::MessageManager::callAsync([this, statusText, progress]() {
                        progressLabel.setText(statusText, juce::dontSendNotification);
                        this->progress = progress;
                        progressBar.repaint();
                    });
                });

            // Finish analyzers
            // Finish the original analysers exactly as before.
            for (auto& analyzer : analyzers) {
                analyzer->finish(outDir);
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

            for (const auto& analyzer : analyzers)
            {
                if (const auto* rawAnalyzer = dynamic_cast<const RawCsvAnalyzer*>(analyzer.get()))
                {
                    waveformSummary = makeWaveformSummary(rawAnalyzer->getResult());
                    temporalSummary = makeTemporalSummary(rawAnalyzer->getResult());
                }

                if (const auto* rmsAnalyzer = dynamic_cast<const RmsPeakAnalyzer*>(analyzer.get()))
                {
                    levelSummary = makeLevelSummary(rmsAnalyzer->getResult());
                    rmsResult = &rmsAnalyzer->getResult();
                }

                if (const auto* thdAnalyzer = dynamic_cast<const ThdAnalyzer*>(analyzer.get()))
                {
                    harmonicSummary = makeHarmonicSummary(thdAnalyzer->getResult());
                    thdResult = &thdAnalyzer->getResult();
                }

                if (const auto* linearAnalyzer = dynamic_cast<const LinearResponseAnalyzer*>(analyzer.get()))
                {
                    toneSummary = makeToneSummary(linearAnalyzer->getResult());
                    linearResult = &linearAnalyzer->getResult();
                }

                if (const auto* transferAnalyzer = dynamic_cast<const TransferCurveAnalyzer*>(analyzer.get()))
                {
                    transferResult = &transferAnalyzer->getResult();
                    nonlinearitySummary = makeNonlinearitySummary(transferAnalyzer->getResult());
                }
            }

            behaviourSummary = makeBehaviourSummary(rmsResult, thdResult, linearResult);
            operatingRangeSummary = makeOperatingRangeSummary(rmsResult, thdResult, linearResult);
            identitySummary = makeIdentitySummary(levelSummary, harmonicSummary, toneSummary, behaviourSummary,
                                                  operatingRangeSummary, nonlinearitySummary, waveformSummary,
                                                  temporalSummary);

            juce::MessageManager::callAsync([this, identitySummary, levelSummary, harmonicSummary, toneSummary, behaviourSummary, operatingRangeSummary, nonlinearitySummary, waveformSummary, temporalSummary]() {
                progressLabel.setText(
                    "Measurement complete",
                    juce::dontSendNotification);
                identityRoleLabel.setText(identitySummary.role, juce::dontSendNotification);
                identityWhyLabel.setText(identitySummary.why, juce::dontSendNotification);
                identityWatchLabel.setText(identitySummary.watch, juce::dontSendNotification);
                levelSummaryLabel.setText(levelSummary.level, juce::dontSendNotification);
                peakSummaryLabel.setText(levelSummary.peaks, juce::dontSendNotification);
                dynamicsSummaryLabel.setText(levelSummary.dynamics, juce::dontSendNotification);
                harmonicSummaryLabel.setText(harmonicSummary.identity, juce::dontSendNotification);
                harmonicBalanceLabel.setText(harmonicSummary.balance, juce::dontSendNotification);
                harmonicGrowthLabel.setText(harmonicSummary.growth, juce::dontSendNotification);
                toneBassLabel.setText(toneSummary.bass, juce::dontSendNotification);
                toneMidLabel.setText(toneSummary.midrange, juce::dontSendNotification);
                toneTrebleLabel.setText(toneSummary.treble, juce::dontSendNotification);
                toneLargestLabel.setText(toneSummary.largest, juce::dontSendNotification);
                behaviourSummaryLabel.setText(behaviourSummary.headline, juce::dontSendNotification);
                behaviourChangesLabel.setText(behaviourSummary.changes, juce::dontSendNotification);
                operatingRangeLabel.setText(operatingRangeSummary.range, juce::dontSendNotification);
                characterStartsLabel.setText(operatingRangeSummary.onset, juce::dontSendNotification);
                underStressLabel.setText(operatingRangeSummary.stress, juce::dontSendNotification);
                nonlinearityLabel.setText(nonlinearitySummary.type, juce::dontSendNotification);
                curveBehaviourLabel.setText(nonlinearitySummary.detail, juce::dontSendNotification);
                waveformShapeLabel.setText(waveformSummary.shape, juce::dontSendNotification);
                waveformKneeLabel.setText(waveformSummary.detail, juce::dontSendNotification);
                waveformStyleLabel.setText(waveformSummary.style, juce::dontSendNotification);
                temporalResponseLabel.setText(temporalSummary.response, juce::dontSendNotification);
                temporalAttackLabel.setText(temporalSummary.attack, juce::dontSendNotification);
                temporalAfterLabel.setText(temporalSummary.after, juce::dontSendNotification);
                progress = 1.0;
                progressBar.repaint();
                runMeasurementButton.setEnabled(true);
                copySummaryButton.setEnabled(true);
            });
        } catch (const std::exception& e) {
            std::cerr << "[Measurement] Exception: " << e.what() << std::endl;
            juce::MessageManager::callAsync([this, e]() {
                showError("Error: " + juce::String(e.what()));
                runMeasurementButton.setEnabled(true);
            });
        } catch (...) {
            std::cerr << "[Measurement] Unknown exception occurred" << std::endl;
            juce::MessageManager::callAsync([this]() {
                showError("Unknown error occurred during measurement");
                runMeasurementButton.setEnabled(true);
            });
        }
    }).detach();
}

Config MainComponent::buildConfigFromUI() {
    Config config;

    // Fill measurement config
    measurementConfig->fillConfig(config);

    // Add parameter buckets from selected parameters
    for (size_t i = 0; i < availableParameters.size(); ++i) {
        if (i < selectedParameters.size() && selectedParameters[i]) {
            // Find corresponding config component
            for (const auto& comp : parameterConfigComponents) {
                if (comp->getConfig().paramName == availableParameters[i]) {
                    config.parameterBuckets.push_back(comp->getConfig());
                    break;
                }
            }
        }
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
    lines.add("PluginDNA Signal Summary");
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

void MainComponent::clearResultsSummary()
{
    copySummaryButton.setEnabled(false);
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
