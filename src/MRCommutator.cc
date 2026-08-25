#include "MRCommutator.hh"
#include "Commutator.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <omp.h>

namespace
{
  int Hermiticity(const Operator &op)
  {
    if (op.IsHermitian())
      return +1;
    if (op.IsAntiHermitian())
      return -1;
    throw std::invalid_argument("MR lambda2 commutator requires Hermitian or anti-Hermitian inputs");
  }

  class ScalarTwoBodyCache
  {
   public:
    ScalarTwoBodyCache(const Operator &X, const Operator &Y,
                       const MRReference &reference)
        : modelspace(X.modelspace), norbits(modelspace->GetNumberOrbits()), max_J(0),
          operators{&X.TwoBody, &Y.TwoBody, &reference.Lambda2}
    {
      for (index_t a : modelspace->all_orbits)
        for (index_t b : modelspace->all_orbits)
          max_J = std::max(max_J,
                           (modelspace->GetOrbit(a).j2 + modelspace->GetOrbit(b).j2) / 2);
      const size_t size = 3 * (max_J + 1) * norbits * norbits * norbits * norbits;
      coupled.assign(size, 0.0);
      pandya.assign(size, 0.0);
      for (int which = 0; which < 3; ++which)
        for (int J = 0; J <= max_J; ++J)
          for (index_t a : modelspace->all_orbits)
            for (index_t b : modelspace->all_orbits)
              for (index_t c : modelspace->all_orbits)
                for (index_t d : modelspace->all_orbits)
                  coupled[Index(which, J, a, b, c, d)] =
                      operators[which]->GetTBME_J(J, J, a, b, c, d);

      for (int which = 0; which < 3; ++which)
        for (int J = 0; J <= max_J; ++J)
          for (index_t a : modelspace->all_orbits)
            for (index_t b : modelspace->all_orbits)
              for (index_t c : modelspace->all_orbits)
                for (index_t d : modelspace->all_orbits)
                  pandya[Index(which, J, a, b, c, d)] =
                      CalculatePandya(which, J, a, b, c, d);
    }

    double Coupled(int which, int J, index_t a, index_t b, index_t c, index_t d) const
    {
      if (J < 0 || J > max_J)
        return 0.0;
      return coupled[Index(which, J, a, b, c, d)];
    }

    double Pandya(int which, int J, index_t a, index_t b, index_t c, index_t d) const
    {
      if (J < 0 || J > max_J)
        return 0.0;
      return pandya[Index(which, J, a, b, c, d)];
    }

    ModelSpace *modelspace;
    size_t norbits;
    int max_J;

   private:
    const TwoBodyME *operators[3];
    std::vector<double> coupled;
    std::vector<double> pandya;

    size_t Index(int which, int J, index_t a, index_t b, index_t c, index_t d) const
    {
      return (((((static_cast<size_t>(which) * (max_J + 1) + J) * norbits + a) *
                  norbits +
                  b) *
                     norbits +
                 c) *
                    norbits +
                d);
    }

    double CalculatePandya(int which, int J, index_t a, index_t b,
                           index_t c, index_t d) const
    {
      const Orbit &oa = modelspace->GetOrbit(a);
      const Orbit &ob = modelspace->GetOrbit(b);
      const Orbit &oc = modelspace->GetOrbit(c);
      const Orbit &od = modelspace->GetOrbit(d);
      const int Jp_min = std::max(std::abs(oa.j2 - od.j2),
                                  std::abs(oc.j2 - ob.j2)) /
                         2;
      const int Jp_max = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
      double value = 0.0;
      for (int Jp = Jp_min; Jp <= Jp_max; ++Jp)
      {
        const double six_j = modelspace->GetSixJ(oa.j2 * 0.5, ob.j2 * 0.5, J,
                                                 oc.j2 * 0.5, od.j2 * 0.5, Jp);
        value -= (2 * Jp + 1.0) * six_j * Coupled(which, Jp, a, d, c, b);
      }
      return value;
    }
  };

  arma::vec RawTerms(index_t one, index_t two, int x_index, int y_index,
                     const ScalarTwoBodyCache &cache)
  {
    const size_t norbits = cache.norbits;
    const Orbit &o1 = cache.modelspace->GetOrbit(one);
    arma::vec result(3, arma::fill::zeros);
    const double jhat1_squared = o1.j2 + 1.0;

    for (int J = 0; J <= cache.max_J; ++J)
    {
      const double Jhat_squared = 2 * J + 1.0;
      for (index_t t = 0; t < norbits; ++t)
        for (index_t s = 0; s < norbits; ++s)
          for (index_t w = 0; w < norbits; ++w)
          {
            const double x_iv = cache.Coupled(x_index, J, one, t, s, w);
            if (x_iv != 0.0)
              for (index_t r = 0; r < norbits; ++r)
                for (index_t v = 0; v < norbits; ++v)
                  result(0) += Jhat_squared / (8.0 * jhat1_squared) * x_iv *
                               cache.Coupled(y_index, J, r, v, two, t) *
                               cache.Coupled(2, J, r, v, s, w);

            for (index_t r = 0; r < norbits; ++r)
            {
              const double x_v = cache.Pandya(x_index, J, one, t, s, r);
              if (x_v == 0.0)
                continue;
              for (index_t v = 0; v < norbits; ++v)
                result(1) += Jhat_squared / (2.0 * jhat1_squared) * x_v *
                             cache.Pandya(2, J, s, r, v, w) *
                             cache.Pandya(y_index, J, v, w, two, t);
            }
          }
    }

    for (int J1 = 0; J1 <= cache.max_J; ++J1)
      for (int J2 = 0; J2 <= cache.max_J; ++J2)
      {
        const double angular_weight = (2 * J1 + 1.0) * (2 * J2 + 1.0);
        for (index_t t = 0; t < norbits; ++t)
        {
          const Orbit &ot = cache.modelspace->GetOrbit(t);
          for (index_t s = 0; s < norbits; ++s)
          {
            if (ot.j2 != cache.modelspace->GetOrbit(s).j2)
              continue;
            const double x_vi = cache.Coupled(x_index, J1, one, t, two, s);
            if (x_vi == 0.0)
              continue;
            for (index_t w = 0; w < norbits; ++w)
              for (index_t r = 0; r < norbits; ++r)
                for (index_t v = 0; v < norbits; ++v)
                  result(2) -= angular_weight /
                               (2.0 * jhat1_squared * (ot.j2 + 1.0)) * x_vi *
                               cache.Coupled(2, J2, s, w, r, v) *
                               cache.Coupled(y_index, J2, r, v, t, w);
          }
        }
      }
    return result;
  }

  double MatrixElement(const arma::mat &matrix, const TwoBodyChannel &channel,
                       index_t a, index_t b, index_t c, index_t d)
  {
    Orbit &oa = channel.modelspace->GetOrbit(a);
    Orbit &ob = channel.modelspace->GetOrbit(b);
    Orbit &oc = channel.modelspace->GetOrbit(c);
    Orbit &od = channel.modelspace->GetOrbit(d);
    if (!channel.CheckChannel_ket(&oa, &ob) || !channel.CheckChannel_ket(&oc, &od))
      return 0.0;
    double phase = 1.0;
    if (a > b)
    {
      phase *= channel.modelspace->GetKet(b, a).Phase(channel.J);
      std::swap(a, b);
    }
    if (c > d)
    {
      phase *= channel.modelspace->GetKet(d, c).Phase(channel.J);
      std::swap(c, d);
    }
    const size_t ibra = channel.GetLocalIndex(a, b);
    const size_t iket = channel.GetLocalIndex(c, d);
    if (a == b)
      phase *= std::sqrt(2.0);
    if (c == d)
      phase *= std::sqrt(2.0);
    return phase * matrix(ibra, iket);
  }

  struct StandardProducts
  {
    std::vector<arma::mat> XLY;
    std::vector<arma::mat> YLX;
    std::vector<arma::mat> LY;
    std::vector<arma::mat> LX;
  };

  StandardProducts BuildStandardProducts(const Operator &X, const Operator &Y,
                                         const MRReference &reference)
  {
    const size_t nch = X.modelspace->GetNumberTwoBodyChannels();
    StandardProducts products;
    products.XLY.resize(nch);
    products.YLX.resize(nch);
    products.LY.resize(nch);
    products.LX.resize(nch);
    for (size_t ch = 0; ch < nch; ++ch)
    {
      const arma::mat &x = X.TwoBody.GetMatrix(ch);
      const arma::mat &y = Y.TwoBody.GetMatrix(ch);
      const arma::mat &lambda = reference.Lambda2.GetMatrix(ch);
      products.LY[ch] = lambda * y;
      products.LX[ch] = lambda * x;
      products.XLY[ch] = x * products.LY[ch];
      products.YLX[ch] = y * products.LX[ch];
    }
    return products;
  }

  std::array<double, 2> PandyaElements(
      const Operator &X, const Operator &Y,
      int J, index_t a, index_t b, index_t c, index_t d)
  {
    ModelSpace &modelspace = *X.modelspace;
    const Orbit &oa = modelspace.GetOrbit(a);
    const Orbit &ob = modelspace.GetOrbit(b);
    const Orbit &oc = modelspace.GetOrbit(c);
    const Orbit &od = modelspace.GetOrbit(d);
    const int Jp_min = std::max(std::abs(oa.j2 - od.j2),
                                std::abs(oc.j2 - ob.j2)) /
                       2;
    const int Jp_max = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
    std::array<double, 2> values{0.0, 0.0};
    for (int Jp = Jp_min; Jp <= Jp_max; ++Jp)
    {
      // All three tensors have exactly the same recoupling coefficient.  The
      // dense imsrg++ cache avoids three hash lookups per standard-coupled
      // matrix element while preserving the Pandya sum and its term order.
      const double weight = (2 * Jp + 1.0) *
                            modelspace.GetCachedSixJ(oa.j2, ob.j2, J,
                                                     oc.j2, od.j2, Jp);
      double x = 0.0;
      double y = 0.0;
      X.TwoBody.GetTBME_J_twoOps(Y.TwoBody, Jp, Jp, a, d, c, b, x, y);
      values[0] -= weight * x;
      values[1] -= weight * y;
    }
    return values;
  }

  double PandyaElement(const TwoBodyME &tensor, ModelSpace &modelspace,
                       int J, index_t a, index_t b, index_t c, index_t d)
  {
    const Orbit &oa = modelspace.GetOrbit(a);
    const Orbit &ob = modelspace.GetOrbit(b);
    const Orbit &oc = modelspace.GetOrbit(c);
    const Orbit &od = modelspace.GetOrbit(d);
    const int Jp_min = std::max(std::abs(oa.j2 - od.j2),
                                std::abs(oc.j2 - ob.j2)) /
                       2;
    const int Jp_max = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
    double value = 0.0;
    for (int Jp = Jp_min; Jp <= Jp_max; ++Jp)
      value -= (2 * Jp + 1.0) *
               modelspace.GetCachedSixJ(oa.j2, ob.j2, J,
                                        oc.j2, od.j2, Jp) *
               tensor.GetTBME_J(Jp, Jp, a, d, c, b);
    return value;
  }

  std::vector<bool> FindLambdaActiveOrbits(const MRReference &reference)
  {
    ModelSpace &modelspace = *reference.modelspace;
    std::vector<bool> active(modelspace.GetNumberOrbits(), false);
    for (size_t ch = 0; ch < modelspace.GetNumberTwoBodyChannels(); ++ch)
    {
      const TwoBodyChannel &channel = modelspace.GetTwoBodyChannel(ch);
      const arma::mat &lambda = reference.Lambda2.GetMatrix(ch);
      for (size_t ibra = 0; ibra < channel.GetNumberKets(); ++ibra)
        for (size_t iket = 0; iket < channel.GetNumberKets(); ++iket)
        {
          if (lambda(ibra, iket) == 0.0)
            continue;
          const Ket &bra = channel.GetKet(ibra);
          const Ket &ket = channel.GetKet(iket);
          active[bra.p] = true;
          active[bra.q] = true;
          active[ket.p] = true;
          active[ket.q] = true;
        }
    }
    return active;
  }

  struct OrderedCrossBlock
  {
    int J;
    int parity;
    int delta_tz;
    size_t norbits;
    std::vector<std::array<index_t, 2>> pairs;
    std::vector<int> pair_index;
    arma::uvec active_pair_indices;
    arma::mat X;
    arma::mat Y;
    arma::mat Lambda;

    int Find(index_t a, index_t b) const
    {
      return pair_index[a * norbits + b];
    }
  };

  OrderedCrossBlock BuildOrderedCrossBlock(
      const Operator &X, const Operator &Y, const MRReference &reference,
      const std::vector<bool> &lambda_active_orbits,
      int J, int parity, int delta_tz)
  {
    ModelSpace &modelspace = *X.modelspace;
    const size_t norbits = modelspace.GetNumberOrbits();
    OrderedCrossBlock block{J, parity, delta_tz, norbits, {},
                            std::vector<int>(norbits * norbits, -1),
                            {}, {}, {}, {}};
    std::vector<arma::uword> active_pair_indices;
    for (index_t a : modelspace.all_orbits)
      for (index_t b : modelspace.all_orbits)
      {
        const Orbit &oa = modelspace.GetOrbit(a);
        const Orbit &ob = modelspace.GetOrbit(b);
        if ((oa.l + ob.l) % 2 != parity ||
            oa.tz2 - ob.tz2 != 2 * delta_tz ||
            std::abs(oa.j2 - ob.j2) > 2 * J || oa.j2 + ob.j2 < 2 * J)
          continue;
        block.pair_index[a * norbits + b] = block.pairs.size();
        block.pairs.push_back({a, b});
        if (lambda_active_orbits[a] && lambda_active_orbits[b])
          active_pair_indices.push_back(block.pairs.size() - 1);
      }
    const size_t dimension = block.pairs.size();
    block.active_pair_indices = arma::conv_to<arma::uvec>::from(active_pair_indices);
    block.X.zeros(dimension, dimension);
    block.Y.zeros(dimension, dimension);
    block.Lambda.zeros(active_pair_indices.size(), active_pair_indices.size());
    for (size_t i = 0; i < dimension; ++i)
      for (size_t j = 0; j < dimension; ++j)
      {
        const auto bra = block.pairs[i];
        const auto ket = block.pairs[j];
        const std::array<double, 2> values = PandyaElements(
            X, Y, J, bra[0], bra[1], ket[0], ket[1]);
        block.X(i, j) = values[0];
        block.Y(i, j) = values[1];
      }
    for (size_t i = 0; i < active_pair_indices.size(); ++i)
      for (size_t j = 0; j < active_pair_indices.size(); ++j)
      {
        const auto bra = block.pairs[active_pair_indices[i]];
        const auto ket = block.pairs[active_pair_indices[j]];
        block.Lambda(i, j) = PandyaElement(
            reference.Lambda2, modelspace, J,
            bra[0], bra[1], ket[0], ket[1]);
      }
    return block;
  }

  MRCommutator::MR1BResult CompleteOneBodyParts(const arma::cube &raw,
                                                int output_hermiticity)
  {
    MRCommutator::MR1BResult result;
    result.IV = raw.slice(0) + output_hermiticity * raw.slice(0).t();
    result.V = raw.slice(1) + output_hermiticity * raw.slice(1).t();
    result.VI = raw.slice(2) + output_hermiticity * raw.slice(2).t();
    return result;
  }
}

namespace MRCommutator
{
  MR1BResult comm221_lambda2_reference(const Operator &X, const Operator &Y,
                                       const MRReference &reference)
  {
    if (X.modelspace != Y.modelspace || X.modelspace != reference.modelspace)
      throw std::invalid_argument("MR lambda2 commutator inputs use different ModelSpace objects");
    if (X.GetJRank() != 0 || Y.GetJRank() != 0 || X.GetTRank() != 0 ||
        Y.GetTRank() != 0 || X.GetParity() != 0 || Y.GetParity() != 0)
      throw std::invalid_argument("MR lambda2 commutator reference supports scalar inputs only");
    const int output_hermiticity = -Hermiticity(X) * Hermiticity(Y);
    ScalarTwoBodyCache cache(X, Y, reference);
    const size_t norbits = X.modelspace->GetNumberOrbits();
    arma::cube raw(norbits, norbits, 3, arma::fill::zeros);
    for (index_t one : X.modelspace->all_orbits)
      for (index_t two : X.modelspace->all_orbits)
      {
        const Orbit &o1 = X.modelspace->GetOrbit(one);
        const Orbit &o2 = X.modelspace->GetOrbit(two);
        if (o1.l != o2.l || o1.j2 != o2.j2 || o1.tz2 != o2.tz2)
          continue;
        raw.tube(one, two) = RawTerms(one, two, 0, 1, cache) -
                             RawTerms(one, two, 1, 0, cache);
      }

    return CompleteOneBodyParts(raw, output_hermiticity);
  }

  MR1BResult comm221_lambda2(const Operator &X, const Operator &Y,
                             const MRReference &reference)
  {
    const double total_start = omp_get_wtime();
    if (X.modelspace != Y.modelspace || X.modelspace != reference.modelspace)
      throw std::invalid_argument("MR lambda2 commutator inputs use different ModelSpace objects");
    if (X.GetJRank() != 0 || Y.GetJRank() != 0 || X.GetTRank() != 0 ||
        Y.GetTRank() != 0 || X.GetParity() != 0 || Y.GetParity() != 0)
      throw std::invalid_argument("MR lambda2 commutator supports scalar inputs only");
    const int output_hermiticity = -Hermiticity(X) * Hermiticity(Y);
    ModelSpace &modelspace = *X.modelspace;
    const size_t norbits = modelspace.GetNumberOrbits();
    int max_J = 0;
    for (index_t a : modelspace.all_orbits)
      for (index_t b : modelspace.all_orbits)
        max_J = std::max(max_J,
                         (modelspace.GetOrbit(a).j2 + modelspace.GetOrbit(b).j2) / 2);
    double section_start = omp_get_wtime();
    const StandardProducts products = BuildStandardProducts(X, Y, reference);
    arma::cube raw(norbits, norbits, 3, arma::fill::zeros);
    X.profiler.timer["MR comm221 lambda2 setup"] += omp_get_wtime() - section_start;

    // IV: ordinary pair-coupled X lambda Y products.
    section_start = omp_get_wtime();
    const size_t nch = modelspace.GetNumberTwoBodyChannels();
    for (size_t ch = 0; ch < nch; ++ch)
    {
      const TwoBodyChannel &channel = modelspace.GetTwoBodyChannel(ch);
      const double angular_weight = 2 * channel.J + 1.0;
      for (index_t one : modelspace.all_orbits)
        for (index_t two : modelspace.all_orbits)
        {
          const Orbit &o1 = modelspace.GetOrbit(one);
          const Orbit &o2 = modelspace.GetOrbit(two);
          if (o1.l != o2.l || o1.j2 != o2.j2 || o1.tz2 != o2.tz2)
            continue;
          for (index_t t : modelspace.all_orbits)
            raw(one, two, 0) += 0.5 * angular_weight / (o1.j2 + 1.0) *
                                (MatrixElement(products.XLY[ch], channel, one, t, two, t) -
                                 MatrixElement(products.YLX[ch], channel, one, t, two, t));
        }
    }
    X.profiler.timer["MR comm221 lambda2 IV"] += omp_get_wtime() - section_start;

    // V: full ordered particle-hole pairs keep the sign of tz_a-tz_b.
    section_start = omp_get_wtime();
    const std::vector<bool> lambda_active_orbits =
        FindLambdaActiveOrbits(reference);
    double v_build_seconds = 0.0;
    double v_blas_seconds = 0.0;
    double v_trace_seconds = 0.0;
    for (int J = 0; J <= max_J; ++J)
      for (int parity = 0; parity <= 1; ++parity)
        for (int delta_tz = -1; delta_tz <= 1; ++delta_tz)
        {
          double v_part_start = omp_get_wtime();
          const OrderedCrossBlock block = BuildOrderedCrossBlock(
              X, Y, reference, lambda_active_orbits,
              J, parity, delta_tz);
          v_build_seconds += omp_get_wtime() - v_part_start;
          if (block.pairs.empty() || block.active_pair_indices.empty())
            continue;
          v_part_start = omp_get_wtime();
          // Lambda has exact support only on pairs whose two orbits occur in
          // at least one nonzero standard-coupled lambda2 element.  Restricting
          // the intermediate index to that principal submatrix therefore
          // preserves X*Lambda*Y exactly while bounding the correlated-reference
          // work by its active space instead of the full single-particle space.
          const arma::mat xly =
              (block.X.cols(block.active_pair_indices) * block.Lambda) *
              block.Y.rows(block.active_pair_indices);
          const arma::mat ylx =
              (block.Y.cols(block.active_pair_indices) * block.Lambda) *
              block.X.rows(block.active_pair_indices);
          v_blas_seconds += omp_get_wtime() - v_part_start;
          const double angular_weight = 2 * block.J + 1.0;
          v_part_start = omp_get_wtime();
          for (index_t one : modelspace.all_orbits)
            for (index_t two : modelspace.all_orbits)
            {
              const Orbit &o1 = modelspace.GetOrbit(one);
              const Orbit &o2 = modelspace.GetOrbit(two);
              if (o1.l != o2.l || o1.j2 != o2.j2 || o1.tz2 != o2.tz2)
                continue;
              for (index_t t : modelspace.all_orbits)
              {
                const int ibra = block.Find(one, t);
                const int iket = block.Find(two, t);
                if (ibra < 0 || iket < 0)
                  continue;
                raw(one, two, 1) += 0.5 * angular_weight / (o1.j2 + 1.0) *
                                    (xly(ibra, iket) - ylx(ibra, iket));
              }
            }
          v_trace_seconds += omp_get_wtime() - v_part_start;
        }
    X.profiler.timer["MR comm221 lambda2 V build"] += v_build_seconds;
    X.profiler.timer["MR comm221 lambda2 V BLAS"] += v_blas_seconds;
    X.profiler.timer["MR comm221 lambda2 V partial trace"] += v_trace_seconds;
    X.profiler.timer["MR comm221 lambda2 V"] += omp_get_wtime() - section_start;

    // VI: lambda*Y and lambda*X remove the innermost (r,v) sum by BLAS.
    // The remaining (J2,w) trace is independent of one,two,J1, so form it
    // once for every allowed (s,t) pair instead of repeating channel lookups.
    section_start = omp_get_wtime();
    const double trace_start = section_start;
    arma::mat ly_trace(norbits, norbits, arma::fill::zeros);
    arma::mat lx_trace(norbits, norbits, arma::fill::zeros);
    for (index_t t : modelspace.all_orbits)
      for (index_t s : modelspace.all_orbits)
      {
        const Orbit &ot = modelspace.GetOrbit(t);
        const Orbit &os = modelspace.GetOrbit(s);
        if (ot.j2 != os.j2)
          continue;
        for (int J2 = 0; J2 <= max_J; ++J2)
          for (index_t w : modelspace.all_orbits)
          {
            const Orbit &ow = modelspace.GetOrbit(w);
            if ((os.l + ow.l) % 2 != (ot.l + ow.l) % 2 ||
                os.tz2 + ow.tz2 != ot.tz2 + ow.tz2)
              continue;
            const int ch = modelspace.GetTwoBodyChannelIndex(
                J2, (os.l + ow.l) % 2, (os.tz2 + ow.tz2) / 2);
            if (ch < 0 || static_cast<size_t>(ch) >= products.LY.size())
              continue;
            const TwoBodyChannel &channel = modelspace.GetTwoBodyChannel(ch);
            const double weight = (2 * J2 + 1.0) / (ot.j2 + 1.0);
            ly_trace(s, t) += weight *
                              MatrixElement(products.LY[ch], channel, s, w, t, w);
            lx_trace(s, t) += weight *
                              MatrixElement(products.LX[ch], channel, s, w, t, w);
          }
      }
    X.profiler.timer["MR comm221 lambda2 VI trace"] +=
        omp_get_wtime() - trace_start;

    for (index_t one : modelspace.all_orbits)
      for (index_t two : modelspace.all_orbits)
      {
        const Orbit &o1 = modelspace.GetOrbit(one);
        const Orbit &o2 = modelspace.GetOrbit(two);
        if (o1.l != o2.l || o1.j2 != o2.j2 || o1.tz2 != o2.tz2)
          continue;
        for (int J1 = 0; J1 <= max_J; ++J1)
          for (index_t t : modelspace.all_orbits)
          {
            const Orbit &ot = modelspace.GetOrbit(t);
            for (index_t s : modelspace.all_orbits)
            {
              if (ot.j2 != modelspace.GetOrbit(s).j2)
                continue;
              const double x = X.TwoBody.GetTBME_J(J1, J1, one, t, two, s);
              const double y = Y.TwoBody.GetTBME_J(J1, J1, one, t, two, s);
              if (x == 0.0 && y == 0.0)
                continue;
              raw(one, two, 2) -= (2 * J1 + 1.0) / (o1.j2 + 1.0) *
                                  (x * ly_trace(s, t) - y * lx_trace(s, t));
            }
          }
      }
    X.profiler.timer["MR comm221 lambda2 VI"] += omp_get_wtime() - section_start;
    MR1BResult result = CompleteOneBodyParts(raw, output_hermiticity);
    X.profiler.timer["MR comm221 lambda2 total"] += omp_get_wtime() - total_start;
    return result;
  }

  Operator Commutator(const Operator &X, const Operator &Y,
                      const MRReference &reference)
  {
    Operator result = ::Commutator::Commutator(X, Y);
    // This exact branch is the SR-degeneration gate: no duplicated SR formula
    // and no floating-point add/subtract of nominally zero MR corrections.
    if (reference.Lambda2.Norm() == 0.0)
      return result;

    const double addon_start = omp_get_wtime();
    const MR1BResult lambda2_one_body = comm221_lambda2(X, Y, reference);
    result.OneBody += lambda2_one_body.Total();
    // Hergert Eq. (49): contract lambda2 with the completed 2B commutator,
    // including 1B--2B as well as both 2B--2B topologies.
    result.ZeroBody += reference.ContractLambda2(result.TwoBody);
    X.profiler.timer["MR lambda2 commutator addons"] +=
        omp_get_wtime() - addon_start;
    return result;
  }
}
