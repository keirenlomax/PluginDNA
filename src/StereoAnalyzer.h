#pragma once
#include "Analyzer.h"
#include "MeasurementResult.h"
#include <cmath>
#include <map>

class StereoAnalyzer final : public Analyzer
{
public:
    explicit StereoAnalyzer(std::vector<juce::String> names) : parameterNames(std::move(names)) { result.analyserName = "StereoDNA"; }
    void processBlock(const BlockContext& ctx) override
    {
        auto& r = runs[ctx.runId];
        if (!r.init) { r.init=true; r.inputGainDb=ctx.inputGainDb; r.params=ctx.paramNamedValues; }
        for (int i=0;i<ctx.numSamples;++i)
        {
            const auto section = static_cast<int>((ctx.firstSample+i) / std::max<int64_t>(1, sectionSamples));
            auto& s = r.sections[juce::jlimit(0,3,section)];
            const double l=ctx.outL?ctx.outL[i]:0.0, rr=ctx.outR?ctx.outR[i]:0.0;
            s.l2+=l*l; s.r2+=rr*rr; s.lr+=l*rr; s.m2+=0.25*(l+rr)*(l+rr); s.s2+=0.25*(l-rr)*(l-rr); ++s.n;
        }
    }
    void setSectionSamples(int64_t n) { sectionSamples=n; }
    void finish(const juce::File&) override
    {
        result.columns={"runId","inputGainDb","leftOnlyLeakDb","rightOnlyLeakDb","midChannelMismatchDb","sideChannelMismatchDb","midToSideLeakDb","sideToMidLeakDb","midCorrelation","sideCorrelation"};
        for (const auto& n:parameterNames) result.columns.push_back(n.toStdString());
        for (const auto& [id,r]:runs)
        {
            auto rms=[](double e,int64_t n){return std::sqrt(e/std::max<int64_t>(1,n));};
            auto db=[](double x){return 20.0*std::log10(std::max(x,1e-12));};
            const auto& L=r.sections[0]; const auto& R=r.sections[1]; const auto& M=r.sections[2]; const auto& S=r.sections[3];
            const double lLeak=db(rms(L.r2,L.n)/std::max(rms(L.l2,L.n),1e-12));
            const double rLeak=db(rms(R.l2,R.n)/std::max(rms(R.r2,R.n),1e-12));
            const double midMismatch=std::abs(db(rms(M.l2,M.n)/std::max(rms(M.r2,M.n),1e-12)));
            const double sideMismatch=std::abs(db(rms(S.l2,S.n)/std::max(rms(S.r2,S.n),1e-12)));
            const double midLeak=db(rms(M.s2,M.n)/std::max(rms(M.m2,M.n),1e-12));
            const double sideLeak=db(rms(S.m2,S.n)/std::max(rms(S.s2,S.n),1e-12));
            auto corr=[](const Stats& x){return x.lr/std::sqrt(std::max(x.l2*x.r2,1e-24));};
            std::vector<double> row{(double)id,r.inputGainDb,lLeak,rLeak,midMismatch,sideMismatch,midLeak,sideLeak,corr(M),corr(S)};
            for (const auto& n:parameterNames){auto it=r.params.find(n);row.push_back(it==r.params.end()?0.0:it->second);} result.rows.push_back(std::move(row));
        }
    }
    const MeasurementDataset& getResult() const noexcept { return result; }
    MeasurementDataset takeResult() noexcept { return std::move(result); }
private:
    struct Stats{double l2=0,r2=0,lr=0,m2=0,s2=0;int64_t n=0;};
    struct Run{bool init=false;double inputGainDb=0;std::map<juce::String,float> params;Stats sections[4];};
    int64_t sectionSamples=1; std::vector<juce::String> parameterNames; std::map<int,Run> runs; MeasurementDataset result;
};
inline std::unique_ptr<Analyzer> createStereoAnalyzer(int64_t sectionSamples,const std::vector<juce::String>& names){auto p=std::make_unique<StereoAnalyzer>(names);p->setSectionSamples(sectionSamples);return p;}
