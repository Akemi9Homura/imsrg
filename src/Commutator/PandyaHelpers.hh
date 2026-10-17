#ifndef IMSRG_PANDYA_HELPERS_HH
#define IMSRG_PANDYA_HELPERS_HH

#include "Operator.hh"
#include "PhysicalConstants.hh"
#include <array>

namespace Commutator
{
namespace detail
{
// Channel-local matrix pointers never survive an operator/reference change.
template<class Index = size_t>
class ChannelMatrices
{
  using Blocks = std::map<std::array<Index, 2>, arma::mat>;
  struct Entry { size_t key = 0; const arma::mat *matrix = nullptr; };
  const Blocks &blocks_;
  size_t nch_;
  std::array<Entry, 256> entries_{}; // 4 KiB, including cached missing blocks
public:
  ChannelMatrices(const Blocks &blocks, size_t nch) : blocks_(blocks), nch_(nch) {}
  const arma::mat *Find(size_t bra, size_t ket)
  {
    const size_t key = bra*nch_ + ket + 1;
    auto &entry = entries_[(key ^ (key >> 8)) & (entries_.size()-1)];
    if (entry.key != key)
    {
      auto it = blocks_.find({static_cast<Index>(bra),static_cast<Index>(ket)});
      entry.matrix = it == blocks_.end() ? nullptr : &it->second;
      entry.key = key;
    }
    return entry.matrix;
  }
  const arma::mat &Get(size_t bra, size_t ket)
  {
    const auto *matrix=Find(bra,ket);
    return matrix ? *matrix : blocks_.at({static_cast<Index>(bra),static_cast<Index>(ket)});
  }
};

// A conditional expression mixing mat and a transpose expression materializes
// a mat even in the non-transposed branch. Keep a reference in that common case.
inline const arma::mat &MatrixOrTranspose(const TwoBodyME &op, size_t bra,
                                          size_t ket, int sign, arma::mat &scratch)
{
  if (bra<=ket) return op.GetMatrix(bra,ket);
  scratch=op.GetMatrix(ket,bra).t()*sign;
  return scratch;
}

// The orbit quartet is fixed throughout a J sum. Prepare its channel quantum
// numbers and global ket indices once; retain the original selection rules.
struct TwoBodyElement
{
  ModelSpace &ms;
  size_t bra, ket;
  int pb, pk, tb, tk, jmin_b, jmax_b, jmin_k, jmax_k;
  int phase_b, phase_k;
  bool flip_b, flip_k;
  double norm = 1;

  TwoBodyElement(ModelSpace &space, int a, int b, int c, int d) : ms(space)
  {
    const auto &oa=ms.GetOrbit(a), &ob=ms.GetOrbit(b), &oc=ms.GetOrbit(c), &od=ms.GetOrbit(d);
    bra = ms.GetKetIndex(std::min(a,b),std::max(a,b));
    ket = ms.GetKetIndex(std::min(c,d),std::max(c,d));
    pb=(oa.l+ob.l)%2; pk=(oc.l+od.l)%2;
    tb=(oa.tz2+ob.tz2)/2; tk=(oc.tz2+od.tz2)/2;
    jmin_b=std::abs(oa.j2-ob.j2)/2; jmax_b=(oa.j2+ob.j2)/2;
    jmin_k=std::abs(oc.j2-od.j2)/2; jmax_k=(oc.j2+od.j2)/2;
    phase_b=ms.GetKet(bra).phase_prefactor; phase_k=ms.GetKet(ket).phase_prefactor;
    flip_b=a>b; flip_k=c>d;
    if (a==b) norm *= PhysConst::SQRT2;
    if (c==d) norm *= PhysConst::SQRT2;
  }

  double Normalized(const TwoBodyME &op, ChannelMatrices<> &blocks, int Jb, int Jk) const
  {
    if (!op.IsAllocated() || (op.parity+pb+pk)%2 || std::abs(tb-tk)!=op.rank_T
        || Jb<jmin_b || Jb>jmax_b || Jk<jmin_k || Jk>jmax_k
        || std::abs(Jb-Jk)>op.rank_J || Jb+Jk<op.rank_J) return 0;
    auto cb=ms.GetTwoBodyChannelIndex(Jb,pb,tb), ck=ms.GetTwoBodyChannelIndex(Jk,pk,tk);
    const auto &bc=ms.GetTwoBodyChannel(cb), &kc=ms.GetTwoBodyChannel(ck);
    auto ib=bc.GetLocalIndex(bra), ik=kc.GetLocalIndex(ket);
    if (ib>=bc.GetNumberKets() || ik>=kc.GetNumberKets()) return 0;
    double phase=1;
    if (flip_b) phase *= phase_b * ms.phase(Jb);
    if (flip_k) phase *= phase_k * ms.phase(Jk);
    if (cb>ck)
    {
      phase *= (op.IsHermitian()?1:-1) * ms.phase(Jb-Jk);
      std::swap(cb,ck); std::swap(ib,ik);
    }
    const auto *matrix=blocks.Find(cb,ck);
    // Preserve GetMatrix's failure for an inconsistent allocated operator.
    return phase * (matrix ? (*matrix)(ib,ik) : op.GetMatrix(cb,ck)(ib,ik));
  }
  double Get(const TwoBodyME &op, ChannelMatrices<> &blocks, int Jb, int Jk) const
  { return norm * Normalized(op,blocks,Jb,Jk); }

  // Ordinary scalar X and Y share a diagonal channel and all index/phase work.
  void ScalarPair(ChannelMatrices<> &xblocks, ChannelMatrices<> &yblocks,
                  int J, double &x, double &y) const
  {
    x=0; y=0;
    if (pb!=pk || tb!=tk || J<jmin_b || J>jmax_b || J<jmin_k || J>jmax_k) return;
    const auto ch=ms.GetTwoBodyChannelIndex(J,pb,tb);
    const auto &channel=ms.GetTwoBodyChannel(ch);
    const auto ib=channel.GetLocalIndex(bra), ik=channel.GetLocalIndex(ket);
    if (ib>=channel.GetNumberKets() || ik>=channel.GetNumberKets()) return;
    double phase=1;
    if (flip_b) phase *= phase_b*ms.phase(J);
    if (flip_k) phase *= phase_k*ms.phase(J);
    x=phase*xblocks.Get(ch,ch)(ib,ik);
    y=phase*yblocks.Get(ch,ch)(ib,ik);
  }
};
}
}
#endif
