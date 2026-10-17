#ifndef IMSRG_SCALAR_PANDYA_PLAN_HH
#define IMSRG_SCALAR_PANDYA_PLAN_HH
#include "PandyaHelpers.hh"
#include "MemoryBudget.hh"
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <mutex>

namespace Commutator { namespace detail {
// Geometry only. Matrix values and occupations are read from each current call.
// The exact structural key also permits reuse across identical ModelSpace copies.
inline std::vector<int64_t> ScalarPandyaKey(ModelSpace &ms)
{
  std::vector<int64_t> key{ms.GetEmax(),static_cast<int64_t>(ms.Orbits.size()),ms.nTwoBodyChannels};
  for (const auto &o:ms.Orbits) for (int v:{o.n,o.l,o.j2,o.tz2}) key.push_back(v);
  key.push_back(ms.Kets.size());
  for (const auto &k:ms.Kets) { key.push_back(k.p); key.push_back(k.q); key.push_back(k.phase_prefactor); }
  for (size_t ch=0;ch<ms.GetNumberTwoBodyChannels();++ch)
  {
    const auto &s=ms.GetTwoBodyChannel(ch);
    const auto &c=ms.GetTwoBodyChannel_CC(ch);
    for (int v:{s.J,s.parity,s.Tz,c.J,c.parity,c.Tz}) key.push_back(v);
    key.push_back(s.KetList.size()); for (auto i:s.KetList) key.push_back(i);
    key.push_back(c.KetList.size()); for (auto i:c.KetList) key.push_back(i);
    key.push_back(c.GetKetIndex_hh().size()); for (auto i:c.GetKetIndex_hh()) key.push_back(i);
    key.push_back(c.GetKetIndex_ph().size()); for (auto i:c.GetKetIndex_ph()) key.push_back(i);
  }
  return key;
}

class ScalarPandyaPlan
{
  // Keep the summation order and phase multiplication separate from coefficients.
  // A term occupies 16 bytes; reservations count against the cache budget.
  struct Term { double coefficient; uint32_t offset; uint16_t channel; int16_t phase; };
  struct ForwardElement { double norm; uint32_t end; };
  struct InverseElement { double norm; uint32_t end1,end2; int exchange,reuse; };
  struct ForwardChannel
  {
    bool ready=false;
    std::vector<Term> terms;
    std::vector<ForwardElement> elements;
  };
  struct InverseRow
  {
    bool ready=false;
    std::vector<Term> terms;
    std::vector<InverseElement> elements;
  };
  std::vector<ForwardChannel> forward_;
  std::vector<std::vector<InverseRow>> inverse_;
  size_t bytes_=0;

  // Conservative reservations avoid capacity growth and bound peak construction
  // storage as well as retained storage. Uncached channels/rows use the old loops.
  static bool Fits(size_t elements,size_t max_terms,size_t element_bytes,size_t remaining)
  {
    const size_t per=element_bytes+max_terms*sizeof(Term);
    return per && elements<=remaining/per && elements<=UINT32_MAX/std::max(size_t(1),max_terms);
  }

  void BuildForward(ModelSpace &ms,size_t ch,size_t &remaining)
  {
    const auto &cc=ms.GetTwoBodyChannel_CC(ch);
    const auto ph=arma::join_cols(cc.GetKetIndex_hh(),cc.GetKetIndex_ph()).eval();
    const size_t nk=cc.GetNumberKets(),ne=2*ph.size()*nk,max_terms=2*ms.GetEmax()+2;
    if (!Fits(ne,max_terms,sizeof(ForwardElement),remaining)) return;
    auto &out=forward_[ch];
    out.elements.reserve(ne); out.terms.reserve(ne*max_terms);
    for (size_t r=0;r<2*ph.size();++r)
    {
      const auto &bra=cc.GetKet(ph[r%ph.size()]);
      const int a=r<ph.size()?bra.p:bra.q,b=r<ph.size()?bra.q:bra.p;
      const auto &oa=ms.GetOrbit(a), &ob=ms.GetOrbit(b);
      for (size_t k=0;k<nk;++k)
      {
        const auto &ket=cc.GetKet(k);
        const int c=ket.p,d=ket.q;
        const auto &oc=ms.GetOrbit(c), &od=ms.GetOrbit(d);
        const TwoBodyElement element(ms,c,b,a,d);
        int low=std::max(std::abs(oa.j2-od.j2),std::abs(oc.j2-ob.j2))/2;
        const int high=std::min(oa.j2+od.j2,oc.j2+ob.j2)/2;
        const int step=(a==d || c==b)?2:1;
        if (step==2) low+=low%2;
        if (element.pb==element.pk && element.tb==element.tk)
          for (int J=low;J<=high;J+=step)
          {
            const double sixj=ms.GetCachedSixJ(oa.j2,ob.j2,cc.J,oc.j2,od.j2,J);
            if (std::abs(sixj)<1e-8) continue;
            const auto source=ms.GetTwoBodyChannelIndex(J,element.pb,element.tb);
            const auto &channel=ms.GetTwoBodyChannel(source);
            const auto ib=channel.GetLocalIndex(element.bra),ik=channel.GetLocalIndex(element.ket);
            if (ib>=channel.GetNumberKets() || ik>=channel.GetNumberKets()) continue;
            int phase=1;
            if (element.flip_b) phase*=element.phase_b*ms.phase(J);
            if (element.flip_k) phase*=element.phase_k*ms.phase(J);
            out.terms.push_back({(2*J+1)*sixj,static_cast<uint32_t>(ib+channel.GetNumberKets()*ik),
                                 static_cast<uint16_t>(source),static_cast<int16_t>(phase)});
          }
        out.elements.push_back({element.norm,static_cast<uint32_t>(out.terms.size())});
      }
    }
    out.ready=true;
    const size_t used=out.terms.capacity()*sizeof(Term)+out.elements.capacity()*sizeof(ForwardElement);
    remaining-=used; bytes_+=used;
  }

  void BuildInverse(ModelSpace &ms,size_t ch,size_t row,size_t &remaining)
  {
    const auto &channel=ms.GetTwoBodyChannel(ch);
    const auto &bra=channel.GetKet(row);
    const int i=bra.p,j=bra.q,J=channel.J;
    const auto &oi=ms.GetOrbit(i), &oj=ms.GetOrbit(j);
    const size_t ne=channel.GetNumberKets()-row,max_terms=4*ms.GetEmax()+4;
    if (!Fits(ne,max_terms,sizeof(InverseElement),remaining)) return;
    auto &out=inverse_[ch][row];
    out.elements.reserve(ne); out.terms.reserve(ne*max_terms);
    for (size_t kidx=row;kidx<channel.GetNumberKets();++kidx)
    {
      const auto &ket=channel.GetKet(kidx);
      const int k=ket.p,l=ket.q;
      const auto &ok=ms.GetOrbit(k), &ol=ms.GetOrbit(l);
      const auto il=ms.GetKetIndex(std::min(i,l),std::max(i,l));
      const auto kj=ms.GetKetIndex(std::min(k,j),std::max(k,j));
      const auto ik=ms.GetKetIndex(std::min(i,k),std::max(i,k));
      const auto lj=ms.GetKetIndex(std::min(l,j),std::max(l,j));
      auto terms=[&](int a,int b,int parity,int tz,size_t pair1,size_t pair2,bool flip1,bool flip2,int low,int high)
      {
        for (int K=low;K<=high;++K)
        {
          const double sixj=ms.GetCachedSixJ(a,b,J,ok.j2,ol.j2,K);
          if (std::abs(sixj)<1e-8) continue;
          const auto source=ms.GetTwoBodyChannelIndex(K,parity,tz);
          const auto &cc=ms.GetTwoBodyChannel_CC(source);
          const size_t n=cc.GetNumberKets();
          const size_t r=cc.GetLocalIndex(pair1)+(flip1?n:0),c=cc.GetLocalIndex(pair2)+(flip2?n:0);
          out.terms.push_back({(2*K+1)*sixj,static_cast<uint32_t>(r+n*c),static_cast<uint16_t>(source),1});
        }
      };
      terms(oi.j2,oj.j2,(oi.l+ol.l)%2,std::abs(oi.tz2-ol.tz2)/2,il,kj,i>l,k>j,
            std::max(std::abs(oi.j2-ol.j2),std::abs(ok.j2-oj.j2))/2,
            std::min(oi.j2+ol.j2,ok.j2+oj.j2)/2);
      const uint32_t end1=out.terms.size();
      int reuse=0;
      if (k==l) reuse=1;
      else if (i==j) reuse=ms.phase((oi.j2+oj.j2+ok.j2+ol.j2)/2);
      else
        terms(oj.j2,oi.j2,(oi.l+ok.l)%2,std::abs(oi.tz2-ok.tz2)/2,ik,lj,i>k,l>j,
              std::max(std::abs(oj.j2-ol.j2),std::abs(ok.j2-oi.j2))/2,
              std::min(oj.j2+ol.j2,ok.j2+oi.j2)/2);
      const double norm=bra.delta_pq()==ket.delta_pq()?1+bra.delta_pq():PhysConst::SQRT2;
      out.elements.push_back({norm,end1,static_cast<uint32_t>(out.terms.size()),ms.phase((ok.j2+ol.j2)/2-J),reuse});
    }
    out.ready=true;
    const size_t used=out.terms.capacity()*sizeof(Term)+out.elements.capacity()*sizeof(InverseElement);
    remaining-=used; bytes_+=used;
  }
public:
  std::vector<int64_t> key;
  explicit ScalarPandyaPlan(ModelSpace &ms,std::vector<int64_t> signature,size_t limit)
      : forward_(ms.GetNumberTwoBodyChannels()),inverse_(ms.GetNumberTwoBodyChannels()),key(std::move(signature))
  {
    bytes_=key.capacity()*sizeof(int64_t)+forward_.capacity()*sizeof(ForwardChannel)+inverse_.capacity()*sizeof(std::vector<InverseRow>);
    for (size_t ch=0;ch<inverse_.size();++ch)
    { inverse_[ch].resize(ms.GetTwoBodyChannel(ch).GetNumberKets()); bytes_+=inverse_[ch].capacity()*sizeof(InverseRow); }
    if (bytes_>limit) return;
    size_t forward_budget=(limit-bytes_)/2,inverse_budget=limit-bytes_-forward_budget;
    for (size_t ch=0;ch<forward_.size();++ch) BuildForward(ms,ch,forward_budget);
    for (size_t ch=0;ch<inverse_.size();++ch)
      for (size_t r=0;r<inverse_[ch].size();++r) BuildInverse(ms,ch,r,inverse_budget);
  }
  size_t Bytes() const { return bytes_; }
  bool Forward(size_t ch,const Operator &x,const Operator &y,arma::mat &xb,arma::mat &yb) const
  {
    const auto &plan=forward_[ch];
    if (!plan.ready) return false;
    const auto &cc=x.modelspace->GetTwoBodyChannel_CC(ch);
    const auto ph=arma::join_cols(cc.GetKetIndex_hh(),cc.GetKetIndex_ph()).eval();
    const size_t nk=cc.GetNumberKets(),np=ph.size();
    xb.zeros(nk,2*np); yb.zeros(2*np,nk);
    std::vector<const double*> xp(forward_.size()),yp(forward_.size());
    for (const auto &m:x.TwoBody.MatEl) if (m.first[0]==m.first[1]) xp[m.first[0]]=m.second.memptr();
    for (const auto &m:y.TwoBody.MatEl) if (m.first[0]==m.first[1]) yp[m.first[0]]=m.second.memptr();
    const int hy=y.IsHermitian()?1:-1;
    size_t t=0,e=0;
    for (size_t r=0;r<2*np;++r)
    {
      const auto &bra=cc.GetKet(ph[r%np]);
      const double na=x.modelspace->GetOrbit(bra.p).occ,nb=x.modelspace->GetOrbit(bra.q).occ;
      const double occupation=r<np?na-nb:nb-na;
      for (size_t k=0;k<nk;++k,++e)
      {
        const auto &el=plan.elements[e]; double vx=0,vy=0;
        for (;t<el.end;++t)
        {
          const auto &term=plan.terms[t];
          const double a=term.phase*xp[term.channel][term.offset],b=term.phase*yp[term.channel][term.offset];
          vx-=term.coefficient*a; vy-=term.coefficient*b*hy;
        }
        xb(k,r)=vx*el.norm*occupation; yb(r,k)=vy*el.norm;
      }
    }
    return true;
  }
  bool Inverse(size_t ch,size_t row,const std::deque<arma::mat> &zbar,Operator &z) const
  {
    const auto &plan=inverse_[ch][row];
    if (!plan.ready) return false;
    auto &matrix=z.TwoBody.GetMatrix(ch,ch);
    const int hz=z.IsHermitian()?1:-1;
    const size_t start=z.IsHermitian()?0:1;
    size_t t=start?plan.elements[0].end2:0;
    for (size_t e=start;e<plan.elements.size();++e)
    {
      const auto &el=plan.elements[e]; double ij=0,ji=0;
      for (;t<el.end1;++t) { const auto &v=plan.terms[t]; ij-=v.coefficient*zbar[v.channel][v.offset]; }
      for (;t<el.end2;++t) { const auto &v=plan.terms[t]; ji-=v.coefficient*zbar[v.channel][v.offset]; }
      if (el.reuse) ji=el.reuse*ij;
      const double value=-(ij-el.exchange*ji)/el.norm;
      matrix(row,row+e)+=value;
      if (e) matrix(row+e,row)+=hz*value;
    }
    return true;
  }
};

inline std::shared_ptr<const ScalarPandyaPlan> PrepareScalarPandya(ModelSpace &ms,size_t limit=MemoryBudgetBytes("IMSRG_SCALAR_PANDYA_CACHE_MB",1024))
{
  if (!limit || ms.GetNumberTwoBodyChannels()>UINT16_MAX) return {};
  for (const auto &ch:ms.TwoBodyChannels)
    if (ch.GetNumberKets()>std::sqrt(static_cast<double>(UINT32_MAX))) return {};
  for (const auto &ch:ms.TwoBodyChannels_CC)
    if (ch.GetNumberKets()>std::sqrt(UINT32_MAX/2.0)) return {};
  std::vector<int64_t> key;
  try { key=ScalarPandyaKey(ms); }
  catch (const std::bad_alloc &) { return {}; }
  if (key.capacity()*sizeof(int64_t)>limit/2) return {};
  static std::mutex mutex;
  static std::shared_ptr<const ScalarPandyaPlan> cache;
  std::lock_guard<std::mutex> lock(mutex);
  if (cache && cache->key==key && cache->Bytes()<=limit) return cache;
  if (cache && cache.use_count()>1) return {};
  cache.reset();
  try { cache=std::make_shared<ScalarPandyaPlan>(ms,std::move(key),limit); }
  catch (const std::bad_alloc &) { return {}; }
  if (cache->Bytes()>limit) { cache.reset(); return {}; }
  return cache;
}
}}
#endif
