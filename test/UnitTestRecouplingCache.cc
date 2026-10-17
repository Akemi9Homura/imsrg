#include "CanonicalSixJ.hh"
#include "ScalarPandyaPlan.hh"
#include "Commutator.hh"
#include "UnitTest.hh"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace
{
void require(bool ok,const char *message)
{
  if (!ok) throw std::runtime_error(message);
}

void compare_flow_kernel(ModelSpace &ms,const char *budget)
{
  UnitTest ut(ms);
  ut.SetRandomSeed(20261017);
  for (int hx:{-1,1}) for (int hy:{-1,1})
  {
    auto x=ut.RandomOp(ms,0,0,0,2,hx);
    auto y=ut.RandomOp(ms,0,0,0,2,hy);
    auto expected=ut.RandomOp(ms,0,0,0,2,-hx*hy),actual=expected;
    setenv("IMSRG_SCALAR_PANDYA_CACHE_MB","0",1);
    Commutator::comm222_phss(x,y,expected);
    setenv("IMSRG_SCALAR_PANDYA_CACHE_MB",budget,1);
    Commutator::comm222_phss(x,y,actual);
    for (const auto &block:expected.TwoBody.MatEl)
      require(block.second.empty() || arma::abs(block.second-actual.TwoBody.GetMatrix(block.first[0],block.first[1])).max()<1e-10,
              "cached scalar commutator differs from uncached result");
  }
}

void check_sixjs(ModelSpace &ms)
{
  ms.PreCalculateSixJ();
  Commutator::detail::CanonicalSixJ cache(ms);
  const int jmax=2*ms.GetEmax()+1;
  for (int a=1;a<=jmax;a+=2) for (int b=1;b<=jmax;b+=2)
    for (int c=1;c<=jmax;c+=2) for (int d=1;d<=jmax;d+=2)
      for (int J=std::max(std::abs(a-b),std::abs(c-d))/2;J<=std::min(a+b,c+d)/2;++J)
        for (int K=std::max(std::abs(a-d),std::abs(c-b))/2;K<=std::min(a+d,c+b)/2;++K)
          require(cache.Get(a,b,J,c,d,K)==ms.GetSixJ(a*.5,b*.5,J,c*.5,d*.5,K),
                  "canonical SixJ cache changed a coefficient");
  Commutator::detail::CanonicalSixJ fallback(ms,0);
  require(fallback.Bytes()==0,"zero budget did not disable SixJ cache");
  require(fallback.Get(1,1,0,1,1,1)==ms.GetSixJ(.5,.5,0,.5,.5,1),"SixJ fallback differs");
}
}

int main()
{
  try
  {
    using namespace Commutator::detail;
    setenv("IMSRG_TEST_MEMORY_BUDGET","-1",1);
    require(MemoryBudgetBytes("IMSRG_TEST_MEMORY_BUDGET",7)==7*1024*1024,"negative budget accepted");
    setenv("IMSRG_TEST_MEMORY_BUDGET","18446744073709551616",1);
    require(MemoryBudgetBytes("IMSRG_TEST_MEMORY_BUDGET",7)==7*1024*1024,"overflowing budget accepted");

    ModelSpace small(2,"He6","He6");
    check_sixjs(small);
    compare_flow_kernel(small,"1024");
    auto pinned=PrepareScalarPandya(small);
    require(pinned && pinned->Bytes()<=1024ULL*1024*1024,"scalar cache budget exceeded");
    ModelSpace copied(small);
    require(PrepareScalarPandya(copied).get()==pinned.get(),"identical geometry did not reuse cache");
    // Occupation values may change without changing the selected hh/ph ket lists.
    copied.GetOrbit(*copied.holes.begin()).occ*=0.75;
    compare_flow_kernel(copied,"1024");
    require(PrepareScalarPandya(copied).get()==pinned.get(),"occupation values should be read live");
    small.SetReference("He4");
    require(!PrepareScalarPandya(small),"retained cache for another reference was reused");
    pinned.reset();
    compare_flow_kernel(small,"1024");

    ModelSpace large(4,"O18","sd-shell");
    // Growth must not invalidate views already in use by a smaller calculation.
    auto angular=PrepareCanonicalSixJs(small);
    large.PreCalculateSixJ();
    {
      CanonicalSixJ during_growth(large);
      require(during_growth.Get(9,9,0,9,9,1)==large.GetSixJ(4.5,4.5,0,4.5,4.5,1),
              "SixJ growth fallback differs");
    }
    angular.reset();
    check_sixjs(large);
    compare_flow_kernel(large,"1"); // Deliberately partial forward/inverse plans.
    auto partial=PrepareScalarPandya(large,1024*1024);
    require(partial && partial->Bytes()<=1024*1024,"partial cache exceeded its budget");
    require(!PrepareScalarPandya(large,0),"zero budget did not disable scalar cache");
    require(!PrepareCanonicalSixJs(large,1),"small budget did not disable SixJ cache");
    std::cout << "Recoupling cache: reference/copy/occupation changes, growth, partial plans and fallbacks passed\n";
    return 0;
  }
  catch (const std::exception &error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
