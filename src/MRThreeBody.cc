#include "MRThreeBody.hh"

#include "PhysicalConstants.hh"

#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <omp.h>

namespace
{
int OscillatorQuanta(const Orbit &orbit) { return 2 * orbit.n + orbit.l; }

struct DensityElement
{
  index_t a;
  index_t b;
  double value;
};

/// lambda3 zero-body term (1/36) sum_J (2J+1) sum V^J Lambda^J over
/// unrestricted (ab)J12 c states, with Lambda rotated to the oscillator basis.
double ContractLambda3(const ThreeBodyME &V3, const MRReference &reference)
{
  ModelSpace &ms = *reference.modelspace;
  const arma::mat &C = reference.NaturalOrbitTransformation;
  const int E3max = ms.GetE3max();
  const size_t norbits = ms.GetNumberOrbits();
  std::vector<std::vector<index_t>> support(norbits);
  for (index_t old_orbit = 0; old_orbit < norbits; ++old_orbit)
    for (index_t new_orbit = 0; new_orbit < norbits; ++new_orbit)
      if (C(old_orbit, new_orbit) != 0.0)
        support[new_orbit].push_back(old_orbit);

  double total = 0.0;
  for (const Lambda3Channel &channel : reference.Lambda3)
  {
    // Unrestricted natural states: (ab) as stored and (ba) with its phase.
    std::vector<std::array<int, 4>> natural;
    std::vector<std::pair<size_t, int>> origin;
    for (size_t i = 0; i < channel.states.size(); ++i)
    {
      const auto &state = channel.states[i];
      natural.push_back(state);
      origin.push_back({i, 1});
      if (state[0] != state[1])
      {
        const int jsum = (ms.GetOrbit(state[0]).j2 + ms.GetOrbit(state[1]).j2) / 2;
        natural.push_back({state[1], state[0], state[2], state[3]});
        origin.push_back({i, ((jsum - state[3] + 1) % 2 == 0) ? 1 : -1});
      }
    }
    arma::mat lambda(natural.size(), natural.size());
    for (size_t i = 0; i < natural.size(); ++i)
      for (size_t j = 0; j < natural.size(); ++j)
        lambda(i, j) = origin[i].second * origin[j].second *
                       channel.matrix(origin[i].first, origin[j].first);

    // T[oscillator state, natural state] = C(a',a) C(b',b) C(c',c).
    std::map<std::array<int, 4>, size_t> oscillator_index;
    std::vector<std::array<int, 4>> oscillator;
    std::vector<std::tuple<size_t, size_t, double>> entries;
    for (size_t k = 0; k < natural.size(); ++k)
    {
      const auto &state = natural[k];
      for (index_t a : support[state[0]])
        for (index_t b : support[state[1]])
          for (index_t c : support[state[2]])
          {
            const std::array<int, 4> key = {static_cast<int>(a), static_cast<int>(b),
                                            static_cast<int>(c), state[3]};
            auto found = oscillator_index.find(key);
            if (found == oscillator_index.end())
            {
              found = oscillator_index.emplace(key, oscillator.size()).first;
              oscillator.push_back(key);
            }
            entries.emplace_back(found->second, k,
                                 C(a, state[0]) * C(b, state[1]) * C(c, state[2]));
          }
    }
    arma::mat T(oscillator.size(), natural.size(), arma::fill::zeros);
    for (const auto &entry : entries)
      T(std::get<0>(entry), std::get<1>(entry)) += std::get<2>(entry);
    const arma::mat lambda_oscillator = T * lambda * T.t();

    // Serial: ThreeBodyME::GetME_pn of the full storage is only used here.
    double channel_sum = 0.0;
    for (size_t i = 0; i < oscillator.size(); ++i)
    {
      const auto &bra = oscillator[i];
      const Orbit &oa = ms.GetOrbit(bra[0]);
      const Orbit &ob = ms.GetOrbit(bra[1]);
      const Orbit &oc = ms.GetOrbit(bra[2]);
      if (OscillatorQuanta(oa) + OscillatorQuanta(ob) + OscillatorQuanta(oc) > E3max)
        continue;
      if (bra[0] == bra[1] && bra[3] % 2 != 0)
        continue;
      for (size_t j = 0; j < oscillator.size(); ++j)
      {
        if (lambda_oscillator(i, j) == 0.0)
          continue;
        const auto &ket = oscillator[j];
        const Orbit &od = ms.GetOrbit(ket[0]);
        const Orbit &oe = ms.GetOrbit(ket[1]);
        const Orbit &of = ms.GetOrbit(ket[2]);
        if (OscillatorQuanta(od) + OscillatorQuanta(oe) + OscillatorQuanta(of) > E3max)
          continue;
        if (ket[0] == ket[1] && ket[3] % 2 != 0)
          continue;
        channel_sum += lambda_oscillator(i, j) *
                       V3.GetME_pn(bra[3], ket[3], channel.twoJ, bra[0], bra[1], bra[2],
                                   ket[0], ket[1], ket[2]);
      }
    }
    total += (channel.twoJ + 1.0) / 36.0 * channel_sum;
  }
  return total;
}
} // namespace

Operator MRNormalOrderThreeBody(const Operator &hbare, const MRReference &reference,
                                MRThreeBodyReport &report, double tolerance)
{
  if (hbare.modelspace != reference.modelspace)
    throw std::invalid_argument("MR-NO2B: the 3N operator and the MR reference use different ModelSpace objects");
  if (hbare.GetParticleRank() < 3 || !hbare.ThreeBody.IsAllocated())
    throw std::invalid_argument("MR-NO2B needs an operator carrying 3N matrix elements");
  if (hbare.GetJRank() != 0 || hbare.GetTRank() != 0 || hbare.GetParity() != 0)
    throw std::invalid_argument("MR-NO2B expects a scalar 3N interaction");
  const ThreeBodyME &V3 = hbare.ThreeBody;
  const std::string mode = V3.GetStorageMode();
  if (mode == "no2b" || mode == "no2bhalf")
    report.storage = "no2b";
  else if (mode == "isospin" || mode == "pn")
    report.storage = "full";
  else
    throw std::invalid_argument("MR-NO2B needs full or no2b 3N matrix elements; got storage mode '" +
                                mode + "' (mono files only hold the HF monopole part)");

  ModelSpace &ms = *reference.modelspace;
  // The pn/isospin 3N storage recouples with ModelSpace six-j symbols; fill the
  // cache before the parallel loops (as ReadWrite and the HF code paths do).
  ms.PreCalculateSixJ();
  const size_t norbits = ms.GetNumberOrbits();
  const int E3max = ms.GetE3max();
  const arma::mat &C = reference.NaturalOrbitTransformation;
  arma::vec occupations(norbits);
  for (size_t p = 0; p < norbits; ++p)
    occupations(p) = reference.occupations[p];
  const arma::mat rho = C * arma::diagmat(occupations) * C.t();

  std::vector<DensityElement> rho_elements;
  std::vector<std::pair<index_t, index_t>> same_channel;  // p <= q in one (l, j, tz) block
  for (index_t a = 0; a < norbits; ++a)
  {
    const Orbit &oa = ms.GetOrbit(a);
    for (index_t b = 0; b < norbits; ++b)
    {
      const Orbit &ob = ms.GetOrbit(b);
      if (oa.l != ob.l || oa.j2 != ob.j2 || oa.tz2 != ob.tz2)
        continue;
      if (a <= b)
        same_channel.push_back({a, b});
      if (std::abs(rho(a, b)) > 1e-14)
        rho_elements.push_back({a, b, rho(a, b)});
    }
  }

  Operator N3(ms, 0, 0, 0, 2);
  N3.SetHermitian();

  // Gamma3 = sum_ab rho_ab V3: the imsrg++ HF formula (HartreeFock::GetNormalOrderedH).
  const size_t nchannels = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ch = 0; ch < nchannels; ++ch)
  {
    TwoBodyChannel &channel = ms.GetTwoBodyChannel(ch);
    const int J = channel.J;
    arma::mat &gamma = N3.TwoBody.GetMatrix(ch);
    const size_t nkets = channel.GetNumberKets();
    for (size_t i = 0; i < nkets; ++i)
    {
      const Ket &bra = channel.GetKet(i);
      const int e_bra = OscillatorQuanta(*bra.op) + OscillatorQuanta(*bra.oq);
      for (size_t j = i; j < nkets; ++j)
      {
        const Ket &ket = channel.GetKet(j);
        const int e_ket = OscillatorQuanta(*ket.op) + OscillatorQuanta(*ket.oq);
        double value = 0.0;
        for (const DensityElement &element : rho_elements)
        {
          if (OscillatorQuanta(ms.GetOrbit(element.a)) + e_bra > E3max ||
              OscillatorQuanta(ms.GetOrbit(element.b)) + e_ket > E3max)
            continue;
          value += element.value * V3.GetME_pn_no2b(bra.p, bra.q, element.a,
                                                    ket.p, ket.q, element.b, J);
        }
        value /= 2 * J + 1;
        if (bra.p == bra.q)
          value /= PhysConst::SQRT2;
        if (ket.p == ket.q)
          value /= PhysConst::SQRT2;
        gamma(i, j) = value;
        gamma(j, i) = value;
      }
    }
  }

  // f3 = 1/4 sum V3 gamma2 with gamma2 = A(g1 g1) + lambda2:
  //   A(g1 g1) part: 1/2 sum rho_ac rho_bd V3(abp, cdq), J12 summed (HF Fock term);
  //   lambda2 part:  1/(2j_p+1) sum_J12 sum_{(ab),(cd)} lambda^J12_norm V3_no2b(abp, cdq; J12)
  //                  in the normalized pair basis, lambda2 rotated to the oscillator basis.
  Operator lambda_operator(ms, 0, 0, 0, 2);
  lambda_operator.SetHermitian();
  lambda_operator.TwoBody = reference.Lambda2;
  const Operator lambda_oscillator = lambda_operator.TransformOneAndTwoBody(C.t(), tolerance);
  struct PairDensityElement
  {
    int J;
    index_t a, b, c, d;
    double value;
  };
  std::vector<PairDensityElement> lambda_elements;
  for (size_t ch = 0; ch < nchannels; ++ch)
  {
    const TwoBodyChannel &channel = ms.GetTwoBodyChannel(ch);
    const arma::mat &matrix = lambda_oscillator.TwoBody.GetMatrix(ch);
    for (size_t i = 0; i < channel.GetNumberKets(); ++i)
      for (size_t j = 0; j < channel.GetNumberKets(); ++j)
        if (matrix(i, j) != 0.0)
        {
          const Ket &bra = channel.GetKet(i);
          const Ket &ket = channel.GetKet(j);
          double norm = 1.0;
          if (bra.p == bra.q)
            norm *= PhysConst::SQRT2;
          if (ket.p == ket.q)
            norm *= PhysConst::SQRT2;
          lambda_elements.push_back({channel.J, bra.p, bra.q, ket.p, ket.q, matrix(i, j) / norm});
        }
  }

  arma::mat f_gamma(norbits, norbits, arma::fill::zeros);
  arma::mat f_lambda(norbits, norbits, arma::fill::zeros);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t k = 0; k < same_channel.size(); ++k)
  {
    const index_t p = same_channel[k].first;
    const index_t q = same_channel[k].second;
    const Orbit &op = ms.GetOrbit(p);
    const int e_p = OscillatorQuanta(op);
    double gamma_sum = 0.0;
    for (const DensityElement &ac : rho_elements)
    {
      const Orbit &oa = ms.GetOrbit(ac.a);
      for (const DensityElement &bd : rho_elements)
      {
        const Orbit &ob = ms.GetOrbit(bd.a);
        if (OscillatorQuanta(oa) + OscillatorQuanta(ob) + e_p > E3max ||
            OscillatorQuanta(ms.GetOrbit(ac.b)) + OscillatorQuanta(ms.GetOrbit(bd.b)) + e_p > E3max)
          continue;
        // Pauli-forbidden pairs (identical orbits, odd J12) are skipped
        // explicitly: their stored value is only the rounding residual of
        // the single-precision isospin components.
        double monopole = 0.0;
        for (int J12 = std::abs(oa.j2 - ob.j2) / 2; J12 <= (oa.j2 + ob.j2) / 2; ++J12)
        {
          if ((ac.a == bd.a || ac.b == bd.b) && J12 % 2 != 0)
            continue;
          monopole += V3.GetME_pn_no2b(ac.a, bd.a, p, ac.b, bd.b, q, J12);
        }
        gamma_sum += ac.value * bd.value * monopole;
      }
    }
    double lambda_sum = 0.0;
    for (const PairDensityElement &element : lambda_elements)
    {
      if (OscillatorQuanta(ms.GetOrbit(element.a)) + OscillatorQuanta(ms.GetOrbit(element.b)) + e_p > E3max ||
          OscillatorQuanta(ms.GetOrbit(element.c)) + OscillatorQuanta(ms.GetOrbit(element.d)) + e_p > E3max)
        continue;
      lambda_sum += element.value * V3.GetME_pn_no2b(element.a, element.b, p,
                                                     element.c, element.d, q, element.J);
    }
    f_gamma(p, q) = f_gamma(q, p) = 0.5 * gamma_sum / (op.j2 + 1.0);
    f_lambda(p, q) = f_lambda(q, p) = lambda_sum / (op.j2 + 1.0);
  }
  N3.OneBody = f_gamma + f_lambda;

  // E3 = 1/36 sum V3 [A(ggg) + A(g lambda2) + lambda3].  A(ggg) follows from the
  // Fock term, (1/3) Tr(rho f_gamma); A(g lambda2) = ContractLambda2(Gamma3)
  // after the rotation to natural orbits.
  double e3_gamma = 0.0;
  for (const auto &pq : same_channel)
  {
    const double weight = (pq.first == pq.second) ? 1.0 : 2.0;
    e3_gamma += weight * (ms.GetOrbit(pq.first).j2 + 1.0) * rho(pq.first, pq.second) *
                f_gamma(pq.first, pq.second) / 3.0;
  }
  report.e3_gamma1_cubed = e3_gamma;
  report.lambda3_in_reference = reference.has_lambda3;
  report.lambda3_term_included = false;
  report.e3_lambda3 = 0.0;
  if (!reference.has_lambda3)
    report.lambda3_reason = "the reference file has no lambda3";
  else if (report.storage != "full")
    report.lambda3_reason = "the 3N file is no2b, which lacks the matrix elements of the lambda3 term";
  else
  {
    report.e3_lambda3 = ContractLambda3(V3, reference);
    report.lambda3_term_included = true;
    report.lambda3_reason.clear();
  }
  N3.ZeroBody = report.e3_gamma1_cubed + report.e3_lambda3;

  Operator N3_natural = N3.TransformOneAndTwoBody(C, tolerance);
  report.e3_gamma1_lambda2 = reference.ContractLambda2(N3_natural.TwoBody);
  N3_natural.ZeroBody += report.e3_gamma1_lambda2;
  return N3_natural;
}
