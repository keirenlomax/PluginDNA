#pragma once
#include "Analyzer.h"
#include "MeasurementResult.h"
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

class V20Analyzer final : public Analyzer {
public:
 enum class Mode { Envelope, Hysteresis, Silence, TruePeak };
 V20Analyzer(Mode m,double sr,double dur,std::vector<juce::String> pn):mode(m),sampleRate(sr),duration(dur),params(std::move(pn)){
   result.analyserName = mode==Mode::Envelope?"DynamicEnvelopeDNA":mode==Mode::Hysteresis?"HysteresisDNA":mode==Mode::Silence?"SelfNoiseDNA":"TruePeakDNA";
 }
 void processBlock(const BlockContext& c) override {
   auto& r=runs[c.runId]; if(!r.init){r.init=true;r.gain=c.inputGainDb;r.p=c.paramNamedValues;}
   for(int i=0;i<c.numSamples;++i){r.in.push_back(c.inL?c.inL[i]:0);r.out.push_back(c.outL?c.outL[i]:0);}
 }
 void finish(const juce::File&) override { build(); }
 const MeasurementDataset& getResult()const{return result;} MeasurementDataset takeResult(){return std::move(result);}
private:
 struct R{bool init=false;float gain=0;std::map<juce::String,float>p;std::vector<float>in,out;};
 static double db(double x){return 20*std::log10(std::max(x,1e-12));}
 static double rms(const std::vector<float>&v,size_t a,size_t b){b=std::min(b,v.size());if(b<=a)return 0;double e=0;for(size_t i=a;i<b;++i)e+=(double)v[i]*v[i];return std::sqrt(e/(b-a));}
 void addParams(std::vector<double>& row,const R&r){row.push_back(r.gain);for(auto&n:params){auto it=r.p.find(n);row.push_back(it==r.p.end()?0:it->second);}}
 void build(){
   if(mode==Mode::Envelope) envelope(); else if(mode==Mode::Hysteresis) hysteresis(); else if(mode==Mode::Silence) silence(); else truepeak();
 }
 void envelope(){
   result.columns={"runId","attackMs","releaseMs","attackOvershootDb","releaseUndershootDb","lowGainDb","highGainDb","inputGainDb"};
   for(auto&n:params)result.columns.push_back(n.toStdString());
   for(auto&kv:runs){auto&r=kv.second;size_t q=r.out.size()/4;if(q<8)continue;
     double lo=rms(r.out,q/2,q), hi=rms(r.out,2*q-q/2,2*q); double targetA=lo+0.9*(hi-lo),targetR=hi-0.9*(hi-lo);
     size_t ai=q,ri=3*q;for(size_t i=q;i<2*q;i++){if(std::abs(r.out[i])>=targetA){ai=i;break;}}for(size_t i=3*q;i<r.out.size();i++){if(std::abs(r.out[i])<=std::max(targetR,1e-12)){ri=i;break;}}
     double pk=*std::max_element(r.out.begin()+q,r.out.begin()+2*q);double mn=*std::min_element(r.out.begin()+3*q,r.out.end());
     std::vector<double>row={(double)kv.first,1000.0*(ai-q)/sampleRate,1000.0*(ri-3*q)/sampleRate,db(std::abs(pk)/std::max(hi,1e-12)),db(std::abs(mn)/std::max(lo,1e-12)),db(lo),db(hi)};addParams(row,r);result.rows.push_back(row);
   }}
 void hysteresis(){
   result.columns={"runId","levelDb","upGainDb","downGainDb","hysteresisDb","inputGainDb"};for(auto&n:params)result.columns.push_back(n.toStdString());
   static const double levels[]={-30,-24,-18,-12,-6};for(auto&kv:runs){auto&r=kv.second;size_t seg=r.out.size()/11;if(seg<8)continue;
    for(int j=0;j<5;j++){int up=j,down=10-j;double inU=rms(r.in,up*seg+(seg/4),(up+1)*seg);double outU=rms(r.out,up*seg+(seg/4),(up+1)*seg);double inD=rms(r.in,down*seg+(seg/4),std::min((down+1)*seg,r.out.size()));double outD=rms(r.out,down*seg+(seg/4),std::min((down+1)*seg,r.out.size()));
      double gu=db(outU/std::max(inU,1e-12)),gd=db(outD/std::max(inD,1e-12));std::vector<double>row={(double)kv.first,levels[j],gu,gd,gd-gu};addParams(row,r);result.rows.push_back(row);}}}
 void silence(){
   result.columns={"runId","outputRmsDbFS","outputPeakDbFS","dcOffset","noiseVariationDb","inputGainDb"};for(auto&n:params)result.columns.push_back(n.toStdString());
   for(auto&kv:runs){auto&r=kv.second;if(r.out.empty())continue;double mean=std::accumulate(r.out.begin(),r.out.end(),0.0)/r.out.size(), rr=rms(r.out,0,r.out.size()),pk=0;for(float x:r.out)pk=std::max(pk,std::abs((double)x));size_t h=r.out.size()/2;double a=rms(r.out,0,h),b=rms(r.out,h,r.out.size());std::vector<double>row={(double)kv.first,db(rr),db(pk),mean,db(std::max(b,1e-12))-db(std::max(a,1e-12))};addParams(row,r);result.rows.push_back(row);}}
 void truepeak(){
   result.columns={"runId","samplePeakDbFS","truePeakDbTP","intersampleOvershootDb","inputGainDb"};for(auto&n:params)result.columns.push_back(n.toStdString());
   for(auto&kv:runs){auto&r=kv.second;if(r.out.size()<16)continue;double sp=0,tp=0;for(float x:r.out)sp=std::max(sp,std::abs((double)x));
     // 8x windowed-sinc reconstruction, 8 taps each side.
     for(size_t i=8;i+8<r.out.size();++i)for(int p=0;p<8;p++){double frac=p/8.0,y=0;for(int k=-7;k<=8;k++){double x=k-frac;double sinc=std::abs(x)<1e-12?1.0:std::sin(juce::MathConstants<double>::pi*x)/(juce::MathConstants<double>::pi*x);double w=0.5+0.5*std::cos(juce::MathConstants<double>::pi*x/8.0);y+=r.out[i+k]*sinc*w;}tp=std::max(tp,std::abs(y));}
     std::vector<double>row={(double)kv.first,db(sp),db(tp),db(tp/std::max(sp,1e-12))};addParams(row,r);result.rows.push_back(row);}}
 Mode mode;double sampleRate,duration;std::vector<juce::String>params;std::map<int,R>runs;MeasurementDataset result;
};
inline std::unique_ptr<Analyzer> createV20Analyzer(V20Analyzer::Mode m,double sr,double dur,const std::vector<juce::String>&p){return std::make_unique<V20Analyzer>(m,sr,dur,p);}
