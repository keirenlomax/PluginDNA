#pragma once

#include "Analyzer.h"
#include "MeasurementResult.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

class BoundaryAnalyzer final : public Analyzer
{
public:
    BoundaryAnalyzer(double sampleRateIn, double durationIn, std::vector<juce::String> paramsIn)
        : sampleRate(sampleRateIn), duration(durationIn), params(std::move(paramsIn))
    { result.analyserName = "BoundaryDNA"; }

    void processBlock(const BlockContext& ctx) override
    {
        auto& run = runs[ctx.runId];
        if (!run.ready)
        {
            run.ready = true; run.inputGainDb = ctx.inputGainDb; run.values = ctx.paramNamedValues;
            init(run);
        }
        const int64_t split = (int64_t)(duration * sampleRate * 0.5);
        for (int i = 0; i < ctx.numSamples; ++i)
        {
            const bool infra = ctx.firstSample + i < split;
            auto& bank = infra ? run.infra : run.ultra;
            const double in = ctx.inL ? ctx.inL[i] : 0.0;
            const double out = ctx.outL ? ctx.outL[i] : 0.0;
            bank.inEnergy += in*in; bank.outMean += out; bank.n++;
            for (auto& d : bank.detectors) d.process(in, out);
        }
    }

    void finish(const juce::File&) override { build(); }
    const MeasurementDataset& getResult() const noexcept { return result; }
    MeasurementDataset takeResult() noexcept { return std::move(result); }

private:
    struct Detector
    {
        double hz=0, c=1,s=0,pc=1,ps=0, inR=0,inI=0,outR=0,outI=0; int64_t n=0; int kind=0;
        void setup(double f,int k,double sr){hz=f;kind=k; const double w=2*juce::MathConstants<double>::pi*f/sr;c=std::cos(w);s=std::sin(w);}
        void process(double in,double out){inR+=in*pc;inI-=in*ps;outR+=out*pc;outI-=out*ps;double nc=pc*c-ps*s;ps=ps*c+pc*s;pc=nc;if((n&4095)==4095){double q=std::hypot(pc,ps);if(q){pc/=q;ps/=q;}}n++;}
        double inMag()const{return n?2.0*std::hypot(inR,inI)/n:0;} double outMag()const{return n?2.0*std::hypot(outR,outI)/n:0;}
        double phase()const{return (std::atan2(outI,outR)-std::atan2(inI,inR))*180.0/juce::MathConstants<double>::pi;}
    };
    struct Bank{std::vector<Detector> detectors;double inEnergy=0,outMean=0;int64_t n=0;};
    struct Run{bool ready=false;float inputGainDb=0;std::map<juce::String,float> values;Bank infra,ultra;};
    void add(Bank& b,double f,int k){if(f<sampleRate*0.49){Detector d;d.setup(f,k,sampleRate);b.detectors.push_back(d);}}
    void init(Run& r)
    {
        for(double f:{0.5,1.0,2.0,5.0,10.0,20.0}) add(r.infra,f,1);
        for(double f:{22000.0,26000.0,32000.0,40000.0}) add(r.ultra,f,2);
        for(double f:{1000.0,4000.0,6000.0,8000.0,10000.0,14000.0,18000.0}) add(r.ultra,f,3);
    }
    static double db(double x){return 20*std::log10(std::max(x,1e-12));}
    void rowsFor(int id,int test,const Bank& b,const Run& r)
    {
        const double rms=b.n?std::sqrt(b.inEnergy/b.n):0;
        for(const auto& d:b.detectors)
        {
            const double in=d.inMag(), out=d.outMag();
            std::vector<double> row{(double)id,(double)test,d.hz,(double)d.kind,
                db(out/std::max(in,1e-12)),d.phase(),db(out/std::max(rms,1e-12)),
                b.n?b.outMean/b.n:0.0,r.inputGainDb};
            for(const auto& p:params){auto it=r.values.find(p);row.push_back(it!=r.values.end()?it->second:0.0);} result.rows.push_back(std::move(row));
        }
    }
    void build()
    {
        result.columns={"runId","testId","frequencyHz","componentClass","transferDb","phaseDegrees","outputDbRelativeToInputRms","dcOffset","inputGainDb"};
        for(const auto& p:params)result.columns.push_back(p.toStdString());
        for(const auto& [id,r]:runs){rowsFor(id,1,r.infra,r);rowsFor(id,2,r.ultra,r);}
    }
    double sampleRate,duration;std::vector<juce::String> params;std::map<int,Run> runs;MeasurementDataset result;
};
inline std::unique_ptr<Analyzer> createBoundaryAnalyzer(double sr,double dur,const std::vector<juce::String>& p){return std::make_unique<BoundaryAnalyzer>(sr,dur,p);}
