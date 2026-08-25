#include "MRCommutator.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

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

    MR1BResult result;
    result.IV = raw.slice(0) + output_hermiticity * raw.slice(0).t();
    result.V = raw.slice(1) + output_hermiticity * raw.slice(1).t();
    result.VI = raw.slice(2) + output_hermiticity * raw.slice(2).t();
    return result;
  }
}
