#pragma once
#include "Analyzer.h"
#include "MeasurementResult.h"
#include <cmath>
#include <map>
#include <vector>

class SummingAnalyzer final : public Analyzer
{
public:
    explicit SummingAnalyzer(std::vector<juce::String> names) : parameterNames(std::move(names)) { result.analyserName="SummingDNA"; }
    void setSectionSamples(int64_t n){sectionSamples=n;}
    void processBlock(const BlockContext& ctx) override
    {
        auto& r=runs[ctx.runId]; if(!r.init){r.init=true;r.inputGainDb=ctx.inputGainDb;r.params=ctx.paramNamedValues;}
        for(int i=0;i<ctx.numSamples;++i){const int section=juce::jlimit(0,2,(int)((ctx.firstSample+i)/std::max<int64_t>(1,sectionSamples)));const float y=ctx.outL?ctx.outL[i]:0.0f;r.y[section].push_back(y);}
    }
    void finish(const juce::File&) override
    {
        result.columns={"runId","inputGainDb","separateRms","summedRms","interactionResidualRms","interactionResidualDbRelative","summedCrestDb","nonAdditivityPercent"};
        for(const auto& n:parameterNames)result.columns.push_back(n.toStdString());
        for(const auto& [id,r]:runs){const size_t n=std::min({r.y[0].size(),r.y[1].size(),r.y[2].size()});double sep2=0,sum2=0,res2=0,pk=0;for(size_t i=0;i<n;++i){double sep=r.y[0][i]+r.y[1][i],sum=r.y[2][i],res=sum-sep;sep2+=sep*sep;sum2+=sum*sum;res2+=res*res;pk=std::max(pk,std::abs(sum));}auto rr=[n](double e){return std::sqrt(e/std::max<size_t>(1,n));};auto db=[](double x){return 20*std::log10(std::max(x,1e-12));};double sr=rr(sep2),mr=rr(sum2),er=rr(res2);std::vector<double> row{(double)id,r.inputGainDb,sr,mr,er,db(er/std::max(mr,1e-12)),db(pk/std::max(mr,1e-12)),100.0*er/std::max(mr,1e-12)};for(const auto& name:parameterNames){auto it=r.params.find(name);row.push_back(it==r.params.end()?0.0:it->second);}result.rows.push_back(std::move(row));}
    }
    const MeasurementDataset& getResult()const noexcept{return result;} MeasurementDataset takeResult()noexcept{return std::move(result);}
private:
    struct Run{bool init=false;double inputGainDb=0;std::map<juce::String,float> params;std::vector<float> y[3];};
    int64_t sectionSamples=1;std::vector<juce::String> parameterNames;std::map<int,Run> runs;MeasurementDataset result;
};
inline std::unique_ptr<Analyzer> createSummingAnalyzer(int64_t sectionSamples,const std::vector<juce::String>& names){auto p=std::make_unique<SummingAnalyzer>(names);p->setSectionSamples(sectionSamples);return p;}
