#ifndef IMSRG_CANONICAL_SIX_J_HH
#define IMSRG_CANONICAL_SIX_J_HH
#include "ModelSpace.hh"
#include "MemoryBudget.hh"
#include <algorithm>
#include <memory>
#include <mutex>
#include <new>

namespace Commutator { namespace detail {
// Angular coefficients only: independent of orbit ordering, occupations and
// operator values. Populate through GetSixJ to retain its exact convention.
class CanonicalSixJTable {
  size_t dim1_, dim2_;
  std::vector<double> values_;
public:
  explicit CanonicalSixJTable(ModelSpace &ms)
    : dim1_(ms.GetEmax()+1), dim2_(2*dim1_),
      values_(4*dim1_*dim1_*dim1_*dim1_*dim1_*dim1_,0.0)
  {
    const int jmax=2*dim1_-1;
    for(int a=1;a<=jmax;a+=2) for(int b=1;b<=jmax;b+=2)
      for(int c=1;c<=jmax;c+=2) for(int d=1;d<=jmax;d+=2)
        for(int J=std::max(std::abs(a-b),std::abs(c-d))/2;J<=std::min(a+b,c+d)/2;++J)
          for(int K=std::max(std::abs(a-d),std::abs(c-b))/2;K<=std::min(a+d,c+b)/2;++K)
            values_[JJSToIndex(a,b,J,c,d,K,dim1_,dim2_)]=ms.GetSixJ(a*.5,b*.5,J,c*.5,d*.5,K);
  }
  int Emax() const { return dim1_-1; }
  size_t Bytes() const { return values_.size()*sizeof(double); }
  double Get(ModelSpace &ms,int a,int b,int J,int c,int d,int K) const
  {
    const int jmax=dim2_-1;
    if(a>0 && b>0 && c>0 && d>0 && a<=jmax && b<=jmax && c<=jmax && d<=jmax
       && (a&1) && (b&1) && (c&1) && (d&1) && J>=0 && K>=0 && J<=jmax && K<=jmax)
      return values_[JJSToIndex(a,b,J,c,d,K,dim1_,dim2_)];
    return ms.GetSixJ(a*.5,b*.5,J,c*.5,d*.5,K);
  }
};

inline std::shared_ptr<const CanonicalSixJTable> PrepareCanonicalSixJs(ModelSpace &ms,
                                              size_t limit=MemoryBudgetBytes("IMSRG_3F2_SIXJ_CACHE_MB",512))
{
  const size_t dim=static_cast<size_t>(ms.GetEmax())+1;
  size_t elements=4;
  for(int k=0;k<6;++k) {
    if(dim==0 || elements>limit/sizeof(double)/dim) return {};
    elements*=dim;
  }
  static std::mutex mutex;
  static std::shared_ptr<const CanonicalSixJTable> cache;
  std::lock_guard<std::mutex> lock(mutex);
  if(cache && cache->Emax()>=ms.GetEmax())
    return cache->Bytes()<=limit ? cache : std::shared_ptr<const CanonicalSixJTable>{};
  // Do not retain two generations during concurrent calls or construction.
  if(cache && cache.use_count()>1)
    return cache->Bytes()<=limit ? cache : std::shared_ptr<const CanonicalSixJTable>{};
  cache.reset();
  try { cache=std::make_shared<CanonicalSixJTable>(ms); }
  catch(const std::bad_alloc &) { return {}; }
  return cache;
}

class CanonicalSixJ {
  ModelSpace &ms_;
  std::shared_ptr<const CanonicalSixJTable> table_;
public:
  explicit CanonicalSixJ(ModelSpace &ms,size_t limit=MemoryBudgetBytes("IMSRG_3F2_SIXJ_CACHE_MB",512))
      : ms_(ms),table_(PrepareCanonicalSixJs(ms,limit)) {}
  double Get(int a,int b,int J,int c,int d,int K) const {
    return table_?table_->Get(ms_,a,b,J,c,d,K):ms_.GetSixJ(a*.5,b*.5,J,c*.5,d*.5,K);
  }
  size_t Bytes() const { return table_?table_->Bytes():0; }
};
}}
#endif
