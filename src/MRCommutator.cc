#include "MRCommutator.hh"
#include "Commutator.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
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

  bool MatrixElementCoordinates(const TwoBodyChannel &channel,
                                index_t a, index_t b, index_t c, index_t d,
                                size_t &ibra, size_t &iket, double &phase)
  {
    Orbit &oa = channel.modelspace->GetOrbit(a);
    Orbit &ob = channel.modelspace->GetOrbit(b);
    Orbit &oc = channel.modelspace->GetOrbit(c);
    Orbit &od = channel.modelspace->GetOrbit(d);
    if (!channel.CheckChannel_ket(&oa, &ob) || !channel.CheckChannel_ket(&oc, &od))
      return false;
    phase = 1.0;
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
    ibra = channel.GetLocalIndex(a, b);
    iket = channel.GetLocalIndex(c, d);
    if (a == b)
      phase *= std::sqrt(2.0);
    if (c == d)
      phase *= std::sqrt(2.0);
    return true;
  }

  void MatrixElements(const arma::mat &first, const arma::mat &second,
                      const TwoBodyChannel &channel,
                      index_t a, index_t b, index_t c, index_t d,
                      double &first_value, double &second_value)
  {
    first_value = 0.0;
    second_value = 0.0;
    size_t ibra = 0;
    size_t iket = 0;
    double phase = 1.0;
    if (!MatrixElementCoordinates(channel, a, b, c, d,
                                  ibra, iket, phase))
      return;
    first_value = phase * first(ibra, iket);
    second_value = phase * second(ibra, iket);
  }

  struct StandardProducts
  {
    std::vector<arma::mat> XLY;
    std::vector<arma::mat> YLX;
    std::vector<arma::mat> LY_active_rows;
    std::vector<arma::mat> LX_active_rows;
    std::vector<arma::mat> X_active_columns;
    std::vector<arma::mat> Y_active_columns;
    std::vector<std::vector<int>> pair_to_active_row;
    std::vector<char> full_active;
    size_t dense_ly_lx_bytes = 0;
    size_t dense_xly_ylx_bytes = 0;

    size_t DataSize() const
    {
      size_t bytes = 0;
      for (size_t ch = 0; ch < XLY.size(); ++ch)
      {
        bytes += sizeof(double) *
                 (XLY[ch].n_elem + YLX[ch].n_elem +
                  LY_active_rows[ch].n_elem + LX_active_rows[ch].n_elem);
        bytes += sizeof(double) *
                 (X_active_columns[ch].n_elem + Y_active_columns[ch].n_elem);
        bytes += sizeof(int) * pair_to_active_row[ch].size();
      }
      return bytes;
    }

    size_t ActiveLYLXDataSize() const
    {
      size_t bytes = 0;
      for (size_t ch = 0; ch < LY_active_rows.size(); ++ch)
        bytes += sizeof(double) *
                 (LY_active_rows[ch].n_elem + LX_active_rows[ch].n_elem);
      return bytes;
    }

    size_t StoredIVOutputDataSize() const
    {
      size_t bytes = 0;
      for (size_t ch = 0; ch < XLY.size(); ++ch)
        bytes += sizeof(double) * (XLY[ch].n_elem + YLX[ch].n_elem);
      return bytes;
    }
  };

  void ActiveRowMatrixElements(const arma::mat &first_active_rows,
                               const arma::mat &second_active_rows,
                               const std::vector<int> &pair_to_active_row,
                               const TwoBodyChannel &channel,
                               index_t a, index_t b, index_t c, index_t d,
                               double &first_value, double &second_value)
  {
    first_value = 0.0;
    second_value = 0.0;
    size_t ibra = 0;
    size_t iket = 0;
    double phase = 1.0;
    if (!MatrixElementCoordinates(channel, a, b, c, d,
                                  ibra, iket, phase))
      return;
    const int active_row = pair_to_active_row[ibra];
    if (active_row < 0)
      return;
    first_value = phase * first_active_rows(active_row, iket);
    second_value = phase * second_active_rows(active_row, iket);
  }

  void SkinnyProductMatrixElements(
      const arma::mat &first_left, const arma::mat &first_right,
      const arma::mat &second_left, const arma::mat &second_right,
      const TwoBodyChannel &channel,
      index_t a, index_t b, index_t c, index_t d,
      double &first_value, double &second_value)
  {
    first_value = 0.0;
    second_value = 0.0;
    size_t ibra = 0;
    size_t iket = 0;
    double phase = 1.0;
    if (!MatrixElementCoordinates(channel, a, b, c, d,
                                  ibra, iket, phase))
      return;
    first_value = phase *
                  arma::dot(first_left.row(ibra), first_right.col(iket));
    second_value = phase *
                   arma::dot(second_left.row(ibra), second_right.col(iket));
  }

  StandardProducts BuildStandardProducts(const Operator &X, const Operator &Y,
                                         const MRReference &reference)
  {
    const size_t nch = X.modelspace->GetNumberTwoBodyChannels();
    StandardProducts products;
    products.XLY.resize(nch);
    products.YLX.resize(nch);
    products.LY_active_rows.resize(nch);
    products.LX_active_rows.resize(nch);
    products.X_active_columns.resize(nch);
    products.Y_active_columns.resize(nch);
    products.pair_to_active_row.resize(nch);
    products.full_active.assign(nch, 0);
    std::vector<size_t> dense_ly_lx_bytes(nch, 0);
    std::vector<size_t> dense_xly_ylx_bytes(nch, 0);
    // Channels are independent and each iteration writes only its own slots,
    // as in the imsrg++ per-channel commutator loops.
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch = 0; ch < nch; ++ch)
    {
      const arma::mat &x = X.TwoBody.GetMatrix(ch);
      const arma::mat &y = Y.TwoBody.GetMatrix(ch);
      const arma::mat &lambda = reference.Lambda2.GetMatrix(ch);
      dense_ly_lx_bytes[ch] =
          2 * sizeof(double) * lambda.n_rows * lambda.n_cols;
      dense_xly_ylx_bytes[ch] =
          sizeof(double) *
          (x.n_rows * y.n_cols + y.n_rows * x.n_cols);
      std::vector<bool> pair_is_active(lambda.n_rows, false);
      size_t active_count = 0;
      for (arma::uword i = 0; i < lambda.n_rows; ++i)
        for (arma::uword j = 0; j < lambda.n_cols; ++j)
        {
          if (lambda(i, j) == 0.0)
            continue;
          if (!pair_is_active[i])
          {
            pair_is_active[i] = true;
            ++active_count;
          }
          if (!pair_is_active[j])
          {
            pair_is_active[j] = true;
            ++active_count;
          }
        }

      if (active_count == 0)
      {
        products.LY_active_rows[ch].zeros(0, y.n_cols);
        products.LX_active_rows[ch].zeros(0, x.n_cols);
        products.pair_to_active_row[ch].assign(lambda.n_rows, -1);
        continue;
      }

      std::vector<arma::uword> active_storage;
      active_storage.reserve(active_count);
      for (arma::uword i = 0; i < pair_is_active.size(); ++i)
        if (pair_is_active[i])
          active_storage.push_back(i);
      const arma::uvec active =
          arma::conv_to<arma::uvec>::from(active_storage);
      products.pair_to_active_row[ch].assign(lambda.n_rows, -1);
      for (size_t i = 0; i < active_storage.size(); ++i)
        products.pair_to_active_row[ch][active_storage[i]] =
            static_cast<int>(i);

      if (active_count == lambda.n_rows)
      {
        products.full_active[ch] = 1;
        products.LY_active_rows[ch] = lambda * y;
        products.LX_active_rows[ch] = lambda * x;
        products.XLY[ch] = x * products.LY_active_rows[ch];
        products.YLX[ch] = y * products.LX_active_rows[ch];
        continue;
      }

      const arma::mat lambda_active = lambda.submat(active, active);
      products.LY_active_rows[ch] = lambda_active * y.rows(active);
      products.LX_active_rows[ch] = lambda_active * x.rows(active);

      // Lambda is exactly zero outside the active principal submatrix.  Keep
      // the exact active rows needed by term VI and the matching active columns
      // needed by the IV partial trace.  No numerical cutoff or MR-IMSRG term
      // is introduced or removed here; the full-active branch above preserves
      // the original complete matrix products.
      products.X_active_columns[ch] = x.cols(active);
      products.Y_active_columns[ch] = y.cols(active);
    }
    for (size_t ch = 0; ch < nch; ++ch)
    {
      products.dense_ly_lx_bytes += dense_ly_lx_bytes[ch];
      products.dense_xly_ylx_bytes += dense_xly_ylx_bytes[ch];
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
    bool full_active;
    arma::mat X;
    arma::mat Y;
    arma::mat X_active_columns;
    arma::mat X_active_rows;
    arma::mat Y_active_columns;
    arma::mat Y_active_rows;
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
                            {}, false, {}, {}, {}, {}, {}, {}, {}};
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
    if (active_pair_indices.empty())
      return block;
    block.full_active = active_pair_indices.size() == dimension;
    block.Lambda.zeros(active_pair_indices.size(), active_pair_indices.size());
    if (block.full_active)
    {
      block.X.zeros(dimension, dimension);
      block.Y.zeros(dimension, dimension);
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
    }
    else
    {
      const size_t active_dimension = active_pair_indices.size();
      block.X_active_columns.zeros(dimension, active_dimension);
      block.Y_active_columns.zeros(dimension, active_dimension);
      block.X_active_rows.zeros(active_dimension, dimension);
      block.Y_active_rows.zeros(active_dimension, dimension);
      for (size_t i = 0; i < dimension; ++i)
        for (size_t j = 0; j < active_dimension; ++j)
        {
          const auto bra = block.pairs[i];
          const auto ket = block.pairs[active_pair_indices[j]];
          const std::array<double, 2> values = PandyaElements(
              X, Y, J, bra[0], bra[1], ket[0], ket[1]);
          block.X_active_columns(i, j) = values[0];
          block.Y_active_columns(i, j) = values[1];
        }
      for (size_t i = 0; i < active_dimension; ++i)
        for (size_t j = 0; j < dimension; ++j)
        {
          const auto bra = block.pairs[active_pair_indices[i]];
          const auto ket = block.pairs[j];
          const std::array<double, 2> values = PandyaElements(
              X, Y, J, bra[0], bra[1], ket[0], ket[1]);
          block.X_active_rows(i, j) = values[0];
          block.Y_active_rows(i, j) = values[1];
        }
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

  // ---------------------------------------------------------------------------
  // Tensor lambda2 one-body terms: X scalar, Y an isospin-conserving tensor of
  // rank L > 0.  Hergert Eq. (50) with x = X, y = Y gives Z = F(X,Y) - F(Y,X),
  //   F(x,y)_ij = 1/4 x_iabc y_deja l_debc + x_iabc y_bejd l_aecd
  //             - 1/2 x_iajb y_cdae l_cdbe + 1/2 x_iajb y_bcde l_acde,
  // which in pair products reads
  //   Z = PT(X l Y - Y l X) + C(X, PT(Y l - l Y)) - D(Y, PT(X l - l X))
  //       + R(X,Y) - R(Y,X),
  // with PT_ij = sum_c W_{ci,cj}, C_ij = sum_ab N_ab X_biaj,
  // D_ij = sum_ab S_ab Y_biaj and the ring R(x,y)_ij = sum x_iabc y_bejd l_aecd.
  // R uses cross-coupled matrices Wbar_{(ab),(cd)} <-> W_{adcb}, with b and d
  // coupled as holes (j_b -m_b) times (-1)^(j_b - m_b):
  //   R(x,y)_ij = -sum_t sum_{JJ'} Jhat J'hat {J J' L; j_j j_i j_t}
  //               (-1)^(j_i+j_t+J'+L) [xbar lbar ybar]^{JJ'}_{(it),(jt)},
  //   scalar  Xbar^J_{(ab),(cd)} = sum_J3 (2J3+1) {j_a j_b J; j_c j_d J3} X^J3_{(ad),(cb)},
  //   tensor  Xbar^{JJ'}_{(ab),(cd)} = sum_{J3 J4} Jhat J'hat J3hat J4hat
  //           (-1)^(j_b+j_d+J'+J4) {j_a j_d J3; j_b j_c J4; J J' L} X^{J3J4}_{(ad),(cb)}.
  // Hermitian conjugation maps PT(Y l X) and R(Y,X) onto PT(X l Y) and R(X,Y):
  // their sum enters as P_ij + h_Z (-1)^(j_i-j_j) P_ji.  Reduced matrix elements
  // follow the imsrg++ (Edmonds) Wigner-Eckart convention and the
  // TwoBodyME::GetTBME_J normalization; every coefficient was checked against
  // brute-force m-scheme sums for L = 1, 2, 3 and both parities.
  //
  // Speed: matrix elements are read through PairTable/BlockAccess (channels and
  // block pointers resolved once, instead of the std::map lookups inside
  // GetTBME_J); 6j/9j coefficients are tabulated per tuple of j values and kept
  // across calls; every product with lambda2 runs over active orbits only;
  // parallel loops write disjoint rows or per-channel buffers that are summed
  // in a fixed order, so the result does not depend on the thread count.

  // Channel, local ket index and GetTBME_J factor (exchange phase, sqrt2 for
  // a == b) of every ordered orbit pair |ab;J>.
  class PairTable
  {
   public:
    struct Entry
    {
      int channel = -1;
      int index = 0;
      double factor = 0.0;
    };

    explicit PairTable(ModelSpace &space)
        : modelspace(&space), norbits(space.GetNumberOrbits()), max_J(0)
    {
      for (index_t a = 0; a < norbits; ++a)
        max_J = std::max(max_J, space.GetOrbit(a).j2);
      entries.resize((max_J + 1) * norbits * norbits);
      for (int J = 0; J <= max_J; ++J)
        for (index_t a = 0; a < norbits; ++a)
          for (index_t b = 0; b < norbits; ++b)
          {
            const Orbit &oa = space.GetOrbit(a);
            const Orbit &ob = space.GetOrbit(b);
            if (2 * J < std::abs(oa.j2 - ob.j2) || 2 * J > oa.j2 + ob.j2)
              continue;
            const int channel = space.GetTwoBodyChannelIndex(J, (oa.l + ob.l) % 2, (oa.tz2 + ob.tz2) / 2);
            const TwoBodyChannel &tbc = space.GetTwoBodyChannel(channel);
            const size_t index = tbc.GetLocalIndex(std::min(a, b), std::max(a, b));
            if (index >= tbc.GetNumberKets())
              continue;
            Entry &entry = entries[Index(J, a, b)];
            entry.channel = channel;
            entry.index = index;
            entry.factor = (a > b ? space.phase((oa.j2 + ob.j2) / 2 + J + 1) : 1) *
                           (a == b ? PhysConst::SQRT2 : 1.0);
          }
    }

    const Entry &operator()(int J, index_t a, index_t b) const { return entries[Index(J, a, b)]; }

    ModelSpace *modelspace;
    size_t norbits;
    int max_J;

   private:
    size_t Index(int J, index_t a, index_t b) const { return (J * norbits + a) * norbits + b; }
    std::vector<Entry> entries;
  };

  // Normalized two-body blocks for any channel order, read as <ab;J1|W|cd;J2>
  // in the GetTBME_J convention.
  class BlockAccess
  {
   public:
    explicit BlockAccess(const PairTable &table)
        : pairs(&table), nchannels(table.modelspace->GetNumberTwoBodyChannels()),
          blocks(nchannels * nchannels) {}

    // Stored blocks of an operator plus their Hermitian conjugates.
    void SetOperator(const Operator &op) { SetTwoBody(op.TwoBody, Hermiticity(op)); }

    void SetTwoBody(const TwoBodyME &two_body, int h)
    {
      ModelSpace &modelspace = *pairs->modelspace;
      for (const auto &entry : two_body.MatEl)
      {
        const size_t bra = entry.first[0], ket = entry.first[1];
        blocks[bra * nchannels + ket] = {&entry.second, false, 1.0};
        if (bra != ket)
          blocks[ket * nchannels + bra] = {
              &entry.second, true,
              double(h * modelspace.phase(modelspace.GetTwoBodyChannel(bra).J -
                                          modelspace.GetTwoBodyChannel(ket).J))};
      }
    }

    void SetBlock(size_t bra, size_t ket, const arma::mat *matrix)
    {
      blocks[bra * nchannels + ket] = {matrix, false, 1.0};
    }

    double operator()(int J1, int J2, index_t a, index_t b, index_t c, index_t d) const
    {
      const PairTable::Entry &bra = (*pairs)(J1, a, b);
      if (bra.channel < 0)
        return 0.0;
      const PairTable::Entry &ket = (*pairs)(J2, c, d);
      if (ket.channel < 0)
        return 0.0;
      const Block &block = blocks[bra.channel * nchannels + ket.channel];
      if (block.matrix == nullptr)
        return 0.0;
      const double value = block.transposed ? (*block.matrix)(ket.index, bra.index)
                                            : (*block.matrix)(bra.index, ket.index);
      return block.sign * bra.factor * ket.factor * value;
    }

   private:
    struct Block
    {
      const arma::mat *matrix = nullptr;
      bool transposed = false;
      double sign = 1.0;
    };
    const PairTable *pairs;
    size_t nchannels;
    std::vector<Block> blocks;
  };

  // 6j/9j coefficients over tuples of single-particle j values (index
  // (j2 - 1) / 2), built once per (largest j, L) and kept for later calls.
  class AngularTables
  {
   public:
    struct SixEntry
    {
      int J, J3;
      double coefficient;
    };
    struct NineEntry
    {
      int J, Jp, J3, J4;
      double coefficient;
    };

    static AngularTables &Get(ModelSpace &modelspace, int L)
    {
      static std::mutex mutex;
      static std::map<std::pair<int, int>, std::unique_ptr<AngularTables>> cache;
      int j2max = 1;
      for (index_t a : modelspace.all_orbits)
        j2max = std::max(j2max, modelspace.GetOrbit(a).j2);
      std::lock_guard<std::mutex> lock(mutex);
      auto &slot = cache[{j2max, L}];
      if (!slot)
        slot.reset(new AngularTables(modelspace, j2max, L));
      return *slot;
    }

    // Jhat1 Jhat2 {J1 J2 L; j_b j_a j_c} (-1)^(j_b+j_c+J1+L): partial trace and D.
    double Trace(int J1, int J2, int ja, int jb, int jc) const
    {
      return trace[(((J1 * (max_J + 1) + J2) * nj + ja) * nj + jb) * nj + jc];
    }
    // (2J+1) {j_a j_b L; j_c j_d J}: contraction C.
    double Contract(int ja, int jb, int jc, int jd, int J) const
    {
      return contract[Quad(ja, jb, jc, jd) * (max_J + 1) + J];
    }
    // Jhat J'hat {J J' L; j_j j_i j_t} (-1)^(j_i+j_t+J'+L): ring trace.
    double RingTrace(int J, int Jp, int ji, int jj, int jt) const
    {
      return ring_trace[(((J * (max_J + 1) + Jp) * nj + ji) * nj + jj) * nj + jt];
    }
    // scalar cross coupling for (j_a, j_b, j_c, j_d)
    const std::vector<SixEntry> &Six(int ja, int jb, int jc, int jd) const
    {
      return six[Quad(ja, jb, jc, jd)];
    }
    // tensor cross coupling; Prepare() must have built the tuple
    const std::vector<NineEntry> &Nine(int ja, int jb, int jc, int jd) const
    {
      return nine[Quad(ja, jb, jc, jd)];
    }
    // Build the 9j lists for the requested tuples (call outside parallel regions).
    void PrepareNine(const std::vector<std::array<int, 4>> &tuples)
    {
      std::vector<size_t> missing;
      for (const auto &t : tuples)
      {
        const size_t q = Quad(t[0], t[1], t[2], t[3]);
        if (!nine_built[q])
        {
          nine_built[q] = 1;
          missing.push_back(q);
        }
      }
#pragma omp parallel for schedule(dynamic, 1)
      for (size_t m = 0; m < missing.size(); ++m)
      {
        size_t q = missing[m];
        const int jd = q % nj, jc = (q / nj) % nj, jb = (q / nj / nj) % nj, ja = q / nj / nj / nj;
        const double a = ja + 0.5, b = jb + 0.5, c = jc + 0.5, d = jd + 0.5;
        std::vector<NineEntry> &entries = nine[q];
        for (int J = std::abs(ja - jb); J <= ja + jb + 1; ++J)
          for (int Jp = std::max(std::abs(jc - jd), std::abs(J - L)); Jp <= std::min(jc + jd + 1, J + L); ++Jp)
            for (int J3 = std::abs(ja - jd); J3 <= ja + jd + 1; ++J3)
              for (int J4 = std::max(std::abs(jc - jb), std::abs(J3 - L)); J4 <= std::min(jc + jb + 1, J3 + L); ++J4)
              {
                const double ninej = AngMom::NineJ(a, d, J3, b, c, J4, J, Jp, L);
                if (std::abs(ninej) < 1e-14)
                  continue;
                entries.push_back({J, Jp, J3, J4,
                                   std::sqrt((2 * J + 1.0) * (2 * Jp + 1.0) * (2 * J3 + 1.0) * (2 * J4 + 1.0)) *
                                       ((jb + jd + 1 + Jp + J4) % 2 == 0 ? 1.0 : -1.0) * ninej});
              }
      }
    }

    int nj, max_J, L;

   private:
    AngularTables(ModelSpace &modelspace, int j2max, int L_)
        : nj((j2max + 1) / 2), max_J(j2max), L(L_)
    {
      modelspace.PreCalculateSixJ();
      const int nJ = max_J + 1;
      trace.assign(nJ * nJ * nj * nj * nj, 0.0);
      ring_trace.assign(nJ * nJ * nj * nj * nj, 0.0);
      contract.assign(nj * nj * nj * nj * nJ, 0.0);
      six.resize(nj * nj * nj * nj);
      nine.resize(nj * nj * nj * nj);
      nine_built.assign(nj * nj * nj * nj, 0);
      for (int J1 = 0; J1 < nJ; ++J1)
        for (int J2 = std::max(0, J1 - L); J2 <= std::min(nJ - 1, J1 + L); ++J2)
          for (int ja = 0; ja < nj; ++ja)
            for (int jb = 0; jb < nj; ++jb)
              for (int jc = 0; jc < nj; ++jc)
              {
                const size_t index = (((J1 * nJ + J2) * nj + ja) * nj + jb) * nj + jc;
                const double hats = std::sqrt((2 * J1 + 1.0) * (2 * J2 + 1.0));
                trace[index] = hats * modelspace.GetSixJ(J1, J2, L, jb + 0.5, ja + 0.5, jc + 0.5) *
                               modelspace.phase(jb + jc + 1 + J1 + L);
                // ring: (J, J', j_i = ja, j_j = jb, j_t = jc)
                ring_trace[index] = hats * modelspace.GetSixJ(J1, J2, L, jb + 0.5, ja + 0.5, jc + 0.5) *
                                    modelspace.phase(ja + jc + 1 + J2 + L);
              }
      for (int ja = 0; ja < nj; ++ja)
        for (int jb = 0; jb < nj; ++jb)
          for (int jc = 0; jc < nj; ++jc)
            for (int jd = 0; jd < nj; ++jd)
            {
              const size_t q = Quad(ja, jb, jc, jd);
              for (int J = 0; J < nJ; ++J)
                contract[q * nJ + J] = (2 * J + 1.0) *
                    modelspace.GetSixJ(ja + 0.5, jb + 0.5, L, jc + 0.5, jd + 0.5, J);
              for (int J = std::abs(ja - jb); J <= ja + jb + 1; ++J)
                for (int J3 = std::max(std::abs(ja - jd), std::abs(jc - jb)); J3 <= std::min(ja + jd, jc + jb) + 1; ++J3)
                {
                  const double sixj = modelspace.GetSixJ(ja + 0.5, jb + 0.5, J, jc + 0.5, jd + 0.5, J3);
                  if (std::abs(sixj) > 1e-14)
                    six[q].push_back({J, J3, (2 * J3 + 1.0) * sixj});
                }
            }
    }

    size_t Quad(int ja, int jb, int jc, int jd) const { return ((ja * nj + jb) * nj + jc) * nj + jd; }

    std::vector<double> trace, ring_trace, contract;
    std::vector<std::vector<SixEntry>> six;
    std::vector<std::vector<NineEntry>> nine;
    std::vector<char> nine_built;
  };

  int JIndex(const Orbit &orbit) { return (orbit.j2 - 1) / 2; }

  // Kets of every two-body channel that contain a given orbit c, written as
  // |c other; J> = factor * stored ket (factor: exchange phase, sqrt2 for c == other).
  struct OrbitKet
  {
    arma::uword index;
    index_t other;
    double factor;
  };
  std::vector<std::vector<std::vector<OrbitKet>>> KetsByOrbit(ModelSpace &modelspace)
  {
    const size_t nchannels = modelspace.GetNumberTwoBodyChannels();
    std::vector<std::vector<std::vector<OrbitKet>>> result(
        nchannels, std::vector<std::vector<OrbitKet>>(modelspace.GetNumberOrbits()));
    for (size_t ch = 0; ch < nchannels; ++ch)
    {
      TwoBodyChannel &channel = modelspace.GetTwoBodyChannel(ch);
      for (size_t k = 0; k < channel.GetNumberKets(); ++k)
      {
        const Ket &ket = channel.GetKet(k);
        if (ket.p == ket.q)
        {
          result[ch][ket.p].push_back({k, static_cast<index_t>(ket.q), PhysConst::SQRT2});
          continue;
        }
        result[ch][ket.p].push_back({k, static_cast<index_t>(ket.q), 1.0});
        const int jsum = (modelspace.GetOrbit(ket.p).j2 + modelspace.GetOrbit(ket.q).j2) / 2;
        result[ch][ket.q].push_back({k, static_cast<index_t>(ket.p),
                                     double(modelspace.phase(jsum + channel.J + 1))});
      }
    }
    return result;
  }

  arma::uvec ActiveKets(const TwoBodyChannel &channel, const std::vector<bool> &active)
  {
    std::vector<arma::uword> kets;
    for (size_t k = 0; k < channel.GetNumberKets(); ++k)
    {
      const Ket &ket = channel.GetKet(k);
      if (active[ket.p] && active[ket.q])
        kets.push_back(k);
    }
    return arma::uvec(kets);
  }

  // Everything one tensor lambda2 evaluation shares.
  struct TensorContext
  {
    TensorContext(const Operator &X_, const Operator &Y_, const MRReference &reference_)
        : X(X_), Y(Y_), reference(reference_), modelspace(*X_.modelspace), L(Y_.GetJRank()),
          norbits(modelspace.GetNumberOrbits()),
          orbits(modelspace.all_orbits.begin(), modelspace.all_orbits.end()),
          active(FindLambdaActiveOrbits(reference_)), pairs(modelspace),
          x_access(pairs), y_access(pairs), lambda_access(pairs),
          tables(AngularTables::Get(modelspace, L)), kets_by_orbit(KetsByOrbit(modelspace)),
          allowed(norbits, std::vector<char>(norbits, 0))
    {
      x_access.SetOperator(X);
      y_access.SetOperator(Y);
      lambda_access.SetTwoBody(reference.Lambda2, +1);
      for (index_t i : orbits)
      {
        const Orbit &oi = modelspace.GetOrbit(i);
        for (index_t j : Y.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
          allowed[i][j] = 1;
      }
      for (size_t ch = 0; ch < modelspace.GetNumberTwoBodyChannels(); ++ch)
        active_kets.push_back(ActiveKets(modelspace.GetTwoBodyChannel(ch), active));
    }

    const Operator &X, &Y;
    const MRReference &reference;
    ModelSpace &modelspace;
    int L;
    size_t norbits;
    std::vector<index_t> orbits;
    std::vector<bool> active;
    PairTable pairs;
    BlockAccess x_access, y_access, lambda_access;
    AngularTables &tables;
    std::vector<std::vector<std::vector<OrbitKet>>> kets_by_orbit;
    std::vector<std::vector<char>> allowed;
    std::vector<arma::uvec> active_kets;
  };

  std::vector<std::pair<size_t, size_t>> TensorChannelPairs(const Operator &Y)
  {
    ModelSpace &modelspace = *Y.modelspace;
    const size_t nchannels = modelspace.GetNumberTwoBodyChannels();
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t ch_bra = 0; ch_bra < nchannels; ++ch_bra)
    {
      const TwoBodyChannel &bra = modelspace.GetTwoBodyChannel(ch_bra);
      if (bra.GetNumberKets() == 0)
        continue;
      for (size_t ch_ket = 0; ch_ket < nchannels; ++ch_ket)
      {
        const TwoBodyChannel &ket = modelspace.GetTwoBodyChannel(ch_ket);
        if (ket.GetNumberKets() == 0 || bra.Tz != ket.Tz ||
            (bra.parity + ket.parity + Y.GetParity()) % 2 != 0 ||
            std::abs(bra.J - ket.J) > Y.GetJRank() || bra.J + ket.J < Y.GetJRank())
          continue;
        pairs.emplace_back(ch_bra, ch_ket);
      }
    }
    return pairs;
  }

  // Stored normalized block of a scalar or tensor operator for either order.
  arma::mat StoredBlock(const Operator &op, size_t ch_bra, size_t ch_ket)
  {
    ModelSpace &modelspace = *op.modelspace;
    const TwoBodyChannel &bra = modelspace.GetTwoBodyChannel(ch_bra);
    const TwoBodyChannel &ket = modelspace.GetTwoBodyChannel(ch_ket);
    const bool ordered = ch_bra <= ch_ket;
    const auto it = op.TwoBody.MatEl.find(
        {ordered ? ch_bra : ch_ket, ordered ? ch_ket : ch_bra});
    if (it == op.TwoBody.MatEl.end())
      return arma::mat(bra.GetNumberKets(), ket.GetNumberKets(), arma::fill::zeros);
    if (ordered)
      return it->second;
    return (Hermiticity(op) * modelspace.phase(bra.J - ket.J)) * it->second.t();
  }

  // Partial trace sum_c W_{ci,cj} over blocks given for (bra, ket) channel
  // pairs; rank L > 0 gives reduced elements (Trace coefficient), L = 0 the
  // ordinary ones.  Each bra channel accumulates into its own buffer.
  arma::mat PartialTrace(const TensorContext &ctx,
                         const std::vector<std::pair<size_t, size_t>> &channel_pairs,
                         const std::vector<arma::mat> &blocks, int L)
  {
    ModelSpace &modelspace = ctx.modelspace;
    const size_t nchannels = modelspace.GetNumberTwoBodyChannels();
    std::vector<std::vector<size_t>> by_bra(nchannels);
    for (size_t p = 0; p < channel_pairs.size(); ++p)
      by_bra[channel_pairs[p].first].push_back(p);
    std::vector<arma::mat> buffers(nchannels);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch_bra = 0; ch_bra < nchannels; ++ch_bra)
    {
      if (by_bra[ch_bra].empty())
        continue;
      arma::mat &buffer = buffers[ch_bra];
      buffer.zeros(ctx.norbits, ctx.norbits);
      const int J1 = modelspace.GetTwoBodyChannel(ch_bra).J;
      for (size_t p : by_bra[ch_bra])
      {
        const size_t ch_ket = channel_pairs[p].second;
        const int J2 = modelspace.GetTwoBodyChannel(ch_ket).J;
        const arma::mat &W = blocks[p];
        for (index_t c : ctx.orbits)
        {
          const auto &bra_kets = ctx.kets_by_orbit[ch_bra][c];
          const auto &ket_kets = ctx.kets_by_orbit[ch_ket][c];
          if (bra_kets.empty() || ket_kets.empty())
            continue;
          const int jc = JIndex(modelspace.GetOrbit(c));
          for (const OrbitKet &bra : bra_kets)
          {
            const Orbit &oi = modelspace.GetOrbit(bra.other);
            for (const OrbitKet &ket : ket_kets)
            {
              if (!ctx.allowed[bra.other][ket.other] && L > 0)
                continue;
              const Orbit &oj = modelspace.GetOrbit(ket.other);
              double coefficient;
              if (L > 0)
                coefficient = ctx.tables.Trace(J1, J2, JIndex(oi), JIndex(oj), jc);
              else if (oi.j2 == oj.j2 && oi.l == oj.l && oi.tz2 == oj.tz2)
                coefficient = (2 * J1 + 1.0) / (oi.j2 + 1.0);
              else
                continue;
              buffer(bra.other, ket.other) += coefficient * bra.factor * ket.factor * W(bra.index, ket.index);
            }
          }
        }
      }
    }
    arma::mat result(ctx.norbits, ctx.norbits, arma::fill::zeros);
    for (const arma::mat &buffer : buffers)
      if (!buffer.empty())
        result += buffer;
    return result;
  }

  // reduced sum_ab N_ab X_biaj: X scalar two-body, N rank-L one-body.
  arma::mat ContractScalarTwoTensorOne(const TensorContext &ctx, const arma::mat &N)
  {
    ModelSpace &modelspace = ctx.modelspace;
    struct Element { index_t a, b; double value; };
    std::vector<Element> nonzero;
    for (index_t a : ctx.orbits)
      for (index_t b : ctx.orbits)
        if (N(a, b) != 0.0)
          nonzero.push_back({a, b, N(a, b)});
    arma::mat result(ctx.norbits, ctx.norbits, arma::fill::zeros);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ii = 0; ii < ctx.orbits.size(); ++ii)
    {
      const index_t i = ctx.orbits[ii];
      const Orbit &oi = modelspace.GetOrbit(i);
      for (index_t j : ctx.orbits)
      {
        if (!ctx.allowed[i][j])
          continue;
        const Orbit &oj = modelspace.GetOrbit(j);
        double sum = 0.0;
        for (const Element &e : nonzero)
        {
          const Orbit &oa = modelspace.GetOrbit(e.a);
          const Orbit &ob = modelspace.GetOrbit(e.b);
          for (int J = std::max(std::abs(oi.j2 - ob.j2), std::abs(oa.j2 - oj.j2)) / 2;
               J <= std::min(oi.j2 + ob.j2, oa.j2 + oj.j2) / 2; ++J)
          {
            const double x = ctx.x_access(J, J, i, e.b, e.a, j);
            if (x != 0.0)
              sum -= e.value * ctx.tables.Contract(JIndex(oi), JIndex(oj), JIndex(oa), JIndex(ob), J) * x;
          }
        }
        result(i, j) = sum;
      }
    }
    return result;
  }

  // reduced sum_ab S_ab Y_biaj: Y rank-L two-body, S scalar one-body.
  arma::mat ContractTensorTwoScalarOne(const TensorContext &ctx, const arma::mat &S)
  {
    ModelSpace &modelspace = ctx.modelspace;
    const int L = ctx.L;
    struct Element { index_t a, b; double value; };
    std::vector<Element> nonzero;
    for (index_t a : ctx.orbits)
      for (index_t b : ctx.orbits)
        if (S(a, b) != 0.0)
          nonzero.push_back({a, b, S(a, b)});
    arma::mat result(ctx.norbits, ctx.norbits, arma::fill::zeros);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ii = 0; ii < ctx.orbits.size(); ++ii)
    {
      const index_t i = ctx.orbits[ii];
      const Orbit &oi = modelspace.GetOrbit(i);
      for (index_t j : ctx.orbits)
      {
        if (!ctx.allowed[i][j])
          continue;
        const Orbit &oj = modelspace.GetOrbit(j);
        double sum = 0.0;
        for (const Element &e : nonzero)
        {
          const Orbit &oa = modelspace.GetOrbit(e.a);
          const Orbit &ob = modelspace.GetOrbit(e.b);
          for (int J1 = std::abs(ob.j2 - oi.j2) / 2; J1 <= (ob.j2 + oi.j2) / 2; ++J1)
            for (int J2 = std::max(std::abs(oa.j2 - oj.j2) / 2, std::abs(J1 - L));
                 J2 <= std::min((oa.j2 + oj.j2) / 2, J1 + L); ++J2)
            {
              const double y = ctx.y_access(J1, J2, e.b, i, e.a, j);
              if (y != 0.0)
                sum += ctx.tables.Trace(J1, J2, JIndex(oi), JIndex(oj), JIndex(oa)) * e.value * y;
            }
        }
        result(i, j) = sum;
      }
    }
    return result;
  }

  // Cross-coupled channels (J, parity, Tz difference) of all ordered orbit
  // pairs, the active-active subset, and rows grouped by the second orbit.
  struct CrossCoupledChannels
  {
    struct Channel
    {
      int J = 0;
      std::vector<std::pair<index_t, index_t>> all, active;
      std::vector<std::vector<arma::uword>> rows_by_second;  // per orbit t
    };

    CrossCoupledChannels(const TensorContext &ctx) : max_J(ctx.pairs.max_J), norbits(ctx.norbits)
    {
      ModelSpace &modelspace = ctx.modelspace;
      channels.resize((max_J + 1) * 6);
      position_all.assign((max_J + 1) * norbits * norbits, -1);
      position_active.assign((max_J + 1) * norbits * norbits, -1);
      for (int J = 0; J <= max_J; ++J)
        for (index_t a : ctx.orbits)
          for (index_t b : ctx.orbits)
          {
            const Orbit &oa = modelspace.GetOrbit(a);
            const Orbit &ob = modelspace.GetOrbit(b);
            if (2 * J < std::abs(oa.j2 - ob.j2) || 2 * J > oa.j2 + ob.j2)
              continue;
            Channel &channel = channels[Id(J, oa, ob)];
            channel.J = J;
            if (channel.rows_by_second.empty())
              channel.rows_by_second.resize(norbits);
            position_all[(J * norbits + a) * norbits + b] = channel.all.size();
            channel.rows_by_second[b].push_back(channel.all.size());
            channel.all.emplace_back(a, b);
            if (ctx.active[a] && ctx.active[b])
            {
              position_active[(J * norbits + a) * norbits + b] = channel.active.size();
              channel.active.emplace_back(a, b);
            }
          }
    }

    static int Id(int J, int parity, int dtz2)  // dtz2 = tz2_a - tz2_b in {-2, 0, 2}
    {
      return (J * 2 + parity) * 3 + (dtz2 / 2 + 1);
    }
    static int Id(int J, const Orbit &oa, const Orbit &ob)
    {
      return Id(J, (oa.l + ob.l) % 2, oa.tz2 - ob.tz2);
    }
    int All(int J, index_t a, index_t b) const { return position_all[(J * norbits + a) * norbits + b]; }
    int Active(int J, index_t a, index_t b) const { return position_active[(J * norbits + a) * norbits + b]; }

    int max_J;
    size_t norbits;
    std::vector<Channel> channels;
    std::vector<int> position_all, position_active;
  };

  // Ring R(X,Y)_ij = sum X_iabc Y_bejd l_aecd (X scalar, Y rank L) with
  // cross-coupled matrices; the lambda2 orbits a, e, c, d are active.
  arma::mat RingTerm(TensorContext &ctx)
  {
    ModelSpace &modelspace = ctx.modelspace;
    const int L = ctx.L;
    const int Y_parity = ctx.Y.GetParity();
    const CrossCoupledChannels cc(ctx);
    const int nchannels = cc.channels.size();
    const int max_J = cc.max_J;

    // lbar and Xbar (rows: all pairs, columns: active pairs) per scalar channel.
    std::vector<arma::mat> lambda_bar(nchannels), x_bar(nchannels);
    for (int id = 0; id < nchannels; ++id)
    {
      const auto &channel = cc.channels[id];
      lambda_bar[id].zeros(channel.active.size(), channel.active.size());
      x_bar[id].zeros(channel.all.size(), channel.active.size());
    }
    std::vector<std::pair<index_t, index_t>> active_pairs, all_pairs;
    for (index_t a : ctx.orbits)
      for (index_t b : ctx.orbits)
      {
        all_pairs.emplace_back(a, b);
        if (ctx.active[a] && ctx.active[b])
          active_pairs.emplace_back(a, b);
      }
    // Xbar^J_{(it),(sr)} = sum_J3 six * X^J3_{(ir),(st)}; also lbar for active rows.
    // Each task owns the rows of one ordered pair (i,t).
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t p = 0; p < all_pairs.size(); ++p)
    {
      const index_t i = all_pairs[p].first, t = all_pairs[p].second;
      const Orbit &oi = modelspace.GetOrbit(i), &ot = modelspace.GetOrbit(t);
      const bool row_active = ctx.active[i] && ctx.active[t];
      std::vector<double> x_values(max_J + 1), l_values(max_J + 1);
      for (const auto &sr : active_pairs)
      {
        const index_t s = sr.first, r = sr.second;
        const Orbit &os = modelspace.GetOrbit(s), &orb = modelspace.GetOrbit(r);
        if ((oi.l + ot.l + os.l + orb.l) % 2 != 0 || oi.tz2 - ot.tz2 != os.tz2 - orb.tz2)
          continue;
        const auto &six = ctx.tables.Six(JIndex(oi), JIndex(ot), JIndex(os), JIndex(orb));
        const int J3min = std::max(std::abs(oi.j2 - orb.j2), std::abs(os.j2 - ot.j2)) / 2;
        const int J3max = std::min(oi.j2 + orb.j2, os.j2 + ot.j2) / 2;
        bool any = false;
        for (int J3 = J3min; J3 <= J3max; ++J3)
        {
          x_values[J3] = ctx.x_access(J3, J3, i, r, s, t);
          l_values[J3] = row_active ? ctx.lambda_access(J3, J3, i, r, s, t) : 0.0;
          any = any || x_values[J3] != 0.0 || l_values[J3] != 0.0;
        }
        if (!any)
          continue;
        for (const auto &entry : six)
        {
          const int id = CrossCoupledChannels::Id(entry.J, oi, ot);
          const int row = cc.All(entry.J, i, t), column = cc.Active(entry.J, s, r);
          if (row < 0 || column < 0)
            continue;
          x_bar[id](row, column) += entry.coefficient * x_values[entry.J3];
          if (row_active)
            lambda_bar[id](cc.Active(entry.J, i, t), column) += entry.coefficient * l_values[entry.J3];
        }
      }
    }

    // Ybar^{JJ'}_{(vw),(jt)} = sum nine * Y^{J3J4}_{(vt),(jw)} for active (v,w):
    // one matrix per (row channel, J').  Each task owns the rows of one (v,w).
    std::vector<std::array<int, 4>> tuples;
    for (const auto &vw : active_pairs)
      for (index_t j : ctx.orbits)
        for (index_t t : ctx.orbits)
          tuples.push_back({JIndex(modelspace.GetOrbit(vw.first)), JIndex(modelspace.GetOrbit(vw.second)),
                            JIndex(modelspace.GetOrbit(j)), JIndex(modelspace.GetOrbit(t))});
    std::sort(tuples.begin(), tuples.end());
    tuples.erase(std::unique(tuples.begin(), tuples.end()), tuples.end());
    ctx.tables.PrepareNine(tuples);
    auto y_bar_index = [&](int id, int Jp) { return id * (2 * L + 1) + (Jp - cc.channels[id].J + L); };
    std::vector<arma::mat> y_bar(nchannels * (2 * L + 1));
    for (int id = 0; id < nchannels; ++id)
    {
      const auto &channel = cc.channels[id];
      if (channel.active.empty())
        continue;
      const int parity = (id / 3) % 2, dtz = id % 3;
      for (int Jp = std::abs(channel.J - L); Jp <= std::min(channel.J + L, max_J); ++Jp)
      {
        const int column_id = (Jp * 2 + (parity + Y_parity) % 2) * 3 + dtz;
        y_bar[y_bar_index(id, Jp)].zeros(channel.active.size(), cc.channels[column_id].all.size());
      }
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t p = 0; p < active_pairs.size(); ++p)
    {
      const index_t v = active_pairs[p].first, w = active_pairs[p].second;
      const Orbit &ov = modelspace.GetOrbit(v), &ow = modelspace.GetOrbit(w);
      std::vector<double> y_values((max_J + 1) * (max_J + 1));
      for (const auto &jt : all_pairs)
      {
        const index_t j = jt.first, t = jt.second;
        const Orbit &oj = modelspace.GetOrbit(j), &ot = modelspace.GetOrbit(t);
        if ((ov.l + ow.l + oj.l + ot.l + Y_parity) % 2 != 0 || ov.tz2 - ow.tz2 != oj.tz2 - ot.tz2)
          continue;
        bool any = false;
        for (int J3 = std::abs(ov.j2 - ot.j2) / 2; J3 <= (ov.j2 + ot.j2) / 2; ++J3)
          for (int J4 = std::max(std::abs(oj.j2 - ow.j2) / 2, std::abs(J3 - L));
               J4 <= std::min((oj.j2 + ow.j2) / 2, J3 + L); ++J4)
          {
            const double y = ctx.y_access(J3, J4, v, t, j, w);
            y_values[J3 * (max_J + 1) + J4] = y;
            any = any || y != 0.0;
          }
        if (!any)
          continue;
        for (const auto &entry : ctx.tables.Nine(JIndex(ov), JIndex(ow), JIndex(oj), JIndex(ot)))
        {
          const double y = y_values[entry.J3 * (max_J + 1) + entry.J4];
          if (y == 0.0)
            continue;
          const int id = CrossCoupledChannels::Id(entry.J, ov, ow);
          y_bar[y_bar_index(id, entry.Jp)](cc.Active(entry.J, v, w), cc.All(entry.Jp, j, t)) +=
              entry.coefficient * y;
        }
      }
    }

    // H = Xbar lbar per channel, then -sum_t RingTrace * [H Ybar]_{(it),(jt)};
    // each row channel accumulates into its own buffer.
    std::vector<arma::mat> buffers(nchannels);
#pragma omp parallel for schedule(dynamic, 1)
    for (int id = 0; id < nchannels; ++id)
    {
      const auto &channel = cc.channels[id];
      if (channel.active.empty() || channel.all.empty())
        continue;
      const arma::mat H = x_bar[id] * lambda_bar[id];
      arma::mat &buffer = buffers[id];
      buffer.zeros(ctx.norbits, ctx.norbits);
      const int parity = (id / 3) % 2, dtz = id % 3;
      for (int Jp = std::abs(channel.J - L); Jp <= std::min(channel.J + L, max_J); ++Jp)
      {
        const auto &column_channel = cc.channels[(Jp * 2 + (parity + Y_parity) % 2) * 3 + dtz];
        const arma::mat &Ybar = y_bar[y_bar_index(id, Jp)];
        if (column_channel.all.empty())
          continue;
        for (index_t t : ctx.orbits)
        {
          const auto &rows = channel.rows_by_second[t];
          const auto &columns = column_channel.rows_by_second[t];
          if (rows.empty() || columns.empty())
            continue;
          const arma::mat G = H.rows(arma::uvec(rows)) * Ybar.cols(arma::uvec(columns));
          const int jt = JIndex(modelspace.GetOrbit(t));
          for (size_t r = 0; r < rows.size(); ++r)
          {
            const index_t i = channel.all[rows[r]].first;
            const int ji = JIndex(modelspace.GetOrbit(i));
            for (size_t c = 0; c < columns.size(); ++c)
            {
              const index_t j = column_channel.all[columns[c]].first;
              if (ctx.allowed[i][j])
                buffer(i, j) -= ctx.tables.RingTrace(channel.J, Jp, ji, JIndex(modelspace.GetOrbit(j)), jt) * G(r, c);
            }
          }
        }
      }
    }
    arma::mat result(ctx.norbits, ctx.norbits, arma::fill::zeros);
    for (const arma::mat &buffer : buffers)
      if (!buffer.empty())
        result += buffer;
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
    const size_t standard_product_bytes = products.DataSize();
    const size_t avoided_dense_ly_lx_bytes =
        products.dense_ly_lx_bytes - products.ActiveLYLXDataSize();
    const size_t avoided_dense_iv_output_bytes =
        products.dense_xly_ylx_bytes - products.StoredIVOutputDataSize();
    X.profiler.counter["MR standard products peak KiB"] = std::max(
        X.profiler.counter["MR standard products peak KiB"],
        static_cast<int>((standard_product_bytes + 1023) / 1024));
    X.profiler.counter["MR dense LY LX avoided KiB"] = std::max(
        X.profiler.counter["MR dense LY LX avoided KiB"],
        static_cast<int>((avoided_dense_ly_lx_bytes + 1023) / 1024));
    X.profiler.counter["MR IV dense output avoided KiB"] = std::max(
        X.profiler.counter["MR IV dense output avoided KiB"],
        static_cast<int>((avoided_dense_iv_output_bytes + 1023) / 1024));
    X.profiler.timer["MR comm221 lambda2 setup"] += omp_get_wtime() - section_start;

    // The loops below follow the imsrg++ commutator pattern: independent
    // outer iterations under "omp parallel for schedule(dynamic)", each
    // writing only its own output slots.  IV and VI parallelize over the
    // first one-body index, so every raw(one,two,.) element is accumulated by
    // one thread in the original term order and the result does not depend on
    // the thread count.
    const std::vector<index_t> orbit_list(modelspace.all_orbits.begin(),
                                          modelspace.all_orbits.end());
    const int norbit_list = static_cast<int>(orbit_list.size());

    // IV: ordinary pair-coupled X lambda Y products.
    section_start = omp_get_wtime();
    const size_t nch = modelspace.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (int io = 0; io < norbit_list; ++io)
    {
      const index_t one = orbit_list[io];
      const Orbit &o1 = modelspace.GetOrbit(one);
      for (size_t ch = 0; ch < nch; ++ch)
      {
        if (!products.full_active[ch] &&
            products.X_active_columns[ch].n_cols == 0)
          continue;
        const TwoBodyChannel &channel = modelspace.GetTwoBodyChannel(ch);
        const double angular_weight = 2 * channel.J + 1.0;
        for (index_t two : X.GetOneBodyChannel(o1.l, o1.j2, o1.tz2))
        {
          for (index_t t : orbit_list)
          {
            double xly = 0.0;
            double ylx = 0.0;
            if (products.full_active[ch])
              MatrixElements(products.XLY[ch], products.YLX[ch], channel,
                             one, t, two, t, xly, ylx);
            else
              SkinnyProductMatrixElements(
                  products.X_active_columns[ch], products.LY_active_rows[ch],
                  products.Y_active_columns[ch], products.LX_active_rows[ch],
                  channel, one, t, two, t, xly, ylx);
            raw(one, two, 0) += 0.5 * angular_weight / (o1.j2 + 1.0) *
                                (xly - ylx);
          }
        }
      }
    }
    X.profiler.timer["MR comm221 lambda2 IV"] += omp_get_wtime() - section_start;

    // V: full ordered particle-hole pairs keep the sign of tz_a-tz_b.
    // Blocks (J, parity, delta_tz) are independent, like the cross-coupled
    // channels of comm222_phss.  Each block's partial trace goes to its own
    // buffer; the buffers are then added in fixed block order, so the result
    // is independent of the thread count.  The V sub-timers report summed
    // thread time.
    section_start = omp_get_wtime();
    const std::vector<bool> lambda_active_orbits =
        FindLambdaActiveOrbits(reference);
    std::vector<std::array<int, 3>> block_keys;
    for (int J = 0; J <= max_J; ++J)
      for (int parity = 0; parity <= 1; ++parity)
        for (int delta_tz = -1; delta_tz <= 1; ++delta_tz)
          block_keys.push_back({J, parity, delta_tz});
    const int nblocks = static_cast<int>(block_keys.size());
    std::vector<arma::mat> v_block_trace(nblocks);
    double v_build_seconds = 0.0;
    double v_blas_seconds = 0.0;
    double v_trace_seconds = 0.0;
    int v_avoided_kib = 0;
#pragma omp parallel for schedule(dynamic, 1) \
    reduction(+ : v_build_seconds, v_blas_seconds, v_trace_seconds) \
    reduction(max : v_avoided_kib)
    for (int ib = 0; ib < nblocks; ++ib)
    {
      double v_part_start = omp_get_wtime();
      const OrderedCrossBlock block = BuildOrderedCrossBlock(
          X, Y, reference, lambda_active_orbits,
          block_keys[ib][0], block_keys[ib][1], block_keys[ib][2]);
      v_build_seconds += omp_get_wtime() - v_part_start;
      if (block.pairs.empty() || block.active_pair_indices.empty())
        continue;
      v_part_start = omp_get_wtime();
      // Lambda has exact support only on pairs whose two orbits occur in
      // at least one nonzero standard-coupled lambda2 element.  Restricting
      // the intermediate index to that principal submatrix therefore
      // preserves X*Lambda*Y exactly while bounding the correlated-reference
      // work by its active space instead of the full single-particle space.
      arma::mat xly;
      arma::mat ylx;
      arma::mat xlambda;
      arma::mat ylambda;
      if (block.full_active)
      {
        xly = (block.X * block.Lambda) * block.Y;
        ylx = (block.Y * block.Lambda) * block.X;
      }
      else
      {
        // The one-body partial trace below reads only entries whose two
        // ordered pairs share the spectator orbit t.  Retain the skinny
        // left factors and evaluate exactly those dots instead of forming
        // two mostly unused dimension-by-dimension output matrices.
        xlambda = block.X_active_columns * block.Lambda;
        ylambda = block.Y_active_columns * block.Lambda;
        const size_t avoided_bytes =
            2 * sizeof(double) * block.pairs.size() * block.pairs.size();
        v_avoided_kib = std::max(
            v_avoided_kib, static_cast<int>((avoided_bytes + 1023) / 1024));
      }
      v_blas_seconds += omp_get_wtime() - v_part_start;
      const double angular_weight = 2 * block.J + 1.0;
      v_part_start = omp_get_wtime();
      arma::mat trace(norbits, norbits, arma::fill::zeros);
      for (index_t one : orbit_list)
      {
        const Orbit &o1 = modelspace.GetOrbit(one);
        for (index_t two : X.GetOneBodyChannel(o1.l, o1.j2, o1.tz2))
        {
          for (index_t t : orbit_list)
          {
            const int ibra = block.Find(one, t);
            const int iket = block.Find(two, t);
            if (ibra < 0 || iket < 0)
              continue;
            const double xly_value = block.full_active
                                         ? xly(ibra, iket)
                                         : arma::dot(
                                               xlambda.row(ibra),
                                               block.Y_active_rows.col(iket));
            const double ylx_value = block.full_active
                                         ? ylx(ibra, iket)
                                         : arma::dot(
                                               ylambda.row(ibra),
                                               block.X_active_rows.col(iket));
            trace(one, two) += 0.5 * angular_weight / (o1.j2 + 1.0) *
                               (xly_value - ylx_value);
          }
        }
      }
      v_block_trace[ib] = std::move(trace);
      v_trace_seconds += omp_get_wtime() - v_part_start;
    }
    for (int ib = 0; ib < nblocks; ++ib)
      if (!v_block_trace[ib].is_empty())
        raw.slice(1) += v_block_trace[ib];
    X.profiler.counter["MR V dense output avoided KiB"] = std::max(
        X.profiler.counter["MR V dense output avoided KiB"], v_avoided_kib);
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
#pragma omp parallel for schedule(dynamic, 1)
    for (int it = 0; it < norbit_list; ++it)
    {
      const index_t t = orbit_list[it];
      const Orbit &ot = modelspace.GetOrbit(t);
      for (index_t s : orbit_list)
      {
        const Orbit &os = modelspace.GetOrbit(s);
        if (ot.j2 != os.j2)
          continue;
        for (int J2 = 0; J2 <= max_J; ++J2)
          for (index_t w : orbit_list)
          {
            const Orbit &ow = modelspace.GetOrbit(w);
            if ((os.l + ow.l) % 2 != (ot.l + ow.l) % 2 ||
                os.tz2 + ow.tz2 != ot.tz2 + ow.tz2)
              continue;
            const int ch = modelspace.GetTwoBodyChannelIndex(
                J2, (os.l + ow.l) % 2, (os.tz2 + ow.tz2) / 2);
            if (ch < 0 || static_cast<size_t>(ch) >= products.XLY.size())
              continue;
            const TwoBodyChannel &channel = modelspace.GetTwoBodyChannel(ch);
            const double weight = (2 * J2 + 1.0) / (ot.j2 + 1.0);
            double ly = 0.0;
            double lx = 0.0;
            ActiveRowMatrixElements(
                products.LY_active_rows[ch], products.LX_active_rows[ch],
                products.pair_to_active_row[ch], channel,
                s, w, t, w, ly, lx);
            ly_trace(s, t) += weight * ly;
            lx_trace(s, t) += weight * lx;
          }
      }
    }
    X.profiler.timer["MR comm221 lambda2 VI trace"] +=
        omp_get_wtime() - trace_start;

#pragma omp parallel for schedule(dynamic, 1)
    for (int io = 0; io < norbit_list; ++io)
    {
      const index_t one = orbit_list[io];
      const Orbit &o1 = modelspace.GetOrbit(one);
      for (index_t two : X.GetOneBodyChannel(o1.l, o1.j2, o1.tz2))
      {
        for (int J1 = 0; J1 <= max_J; ++J1)
          for (index_t t : orbit_list)
          {
            const Orbit &ot = modelspace.GetOrbit(t);
            if (2 * J1 < std::abs(o1.j2 - ot.j2) ||
                2 * J1 > o1.j2 + ot.j2)
              continue;
            for (index_t s : orbit_list)
            {
              if (ot.j2 != modelspace.GetOrbit(s).j2)
                continue;
              const double ly = ly_trace(s, t);
              const double lx = lx_trace(s, t);
              if (ly == 0.0 && lx == 0.0)
                continue;
              double x = 0.0;
              double y = 0.0;
              X.TwoBody.GetTBME_J_twoOps(
                  Y.TwoBody, J1, J1, one, t, two, s, x, y);
              if (x == 0.0 && y == 0.0)
                continue;
              raw(one, two, 2) -= (2 * J1 + 1.0) / (o1.j2 + 1.0) *
                                  (x * ly - y * lx);
            }
          }
      }
    }
    X.profiler.timer["MR comm221 lambda2 VI"] += omp_get_wtime() - section_start;
    MR1BResult result = CompleteOneBodyParts(raw, output_hermiticity);
    X.profiler.timer["MR comm221 lambda2 total"] += omp_get_wtime() - total_start;
    return result;
  }

  arma::mat comm221_lambda2_tensor(const Operator &X, const Operator &Y,
                                   const MRReference &reference)
  {
    const auto is_scalar = [](const Operator &op)
    { return op.GetJRank() == 0 && op.GetTRank() == 0 && op.GetParity() == 0; };
    if (!is_scalar(X) && is_scalar(Y))
      return -comm221_lambda2_tensor(Y, X, reference);
    if (X.modelspace != Y.modelspace || X.modelspace != reference.modelspace)
      throw std::invalid_argument("MR lambda2 commutator inputs use different ModelSpace objects");
    if (!is_scalar(X) || Y.GetJRank() == 0 || Y.GetTRank() != 0)
      throw std::invalid_argument(
          "MR lambda2 tensor terms need a scalar and an isospin-conserving tensor of rank > 0");
    const double total_start = omp_get_wtime();
    TensorContext ctx(X, Y, reference);
    ModelSpace &modelspace = ctx.modelspace;
    const int L = ctx.L;

    // Products with lambda2 over active kets only: P = X l Y and W2 = Y l - l Y
    // for every tensor channel order, W3 = X l - l X per scalar channel.
    double section_start = omp_get_wtime();
    const auto channel_pairs = TensorChannelPairs(Y);
    std::vector<arma::mat> P(channel_pairs.size()), W2(channel_pairs.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t p = 0; p < channel_pairs.size(); ++p)
    {
      const size_t ch_bra = channel_pairs[p].first, ch_ket = channel_pairs[p].second;
      const arma::uvec &bra_active = ctx.active_kets[ch_bra];
      const arma::uvec &ket_active = ctx.active_kets[ch_ket];
      const arma::mat y = StoredBlock(Y, ch_bra, ch_ket);
      P[p].zeros(y.n_rows, y.n_cols);
      W2[p].zeros(y.n_rows, y.n_cols);
      if (!bra_active.empty())
      {
        const arma::mat ly = reference.Lambda2.GetMatrix(ch_bra, ch_bra)(bra_active, bra_active) *
                             y.rows(bra_active);
        P[p] = X.TwoBody.GetMatrix(ch_bra, ch_bra).cols(bra_active) * ly;
        W2[p].rows(bra_active) -= ly;
      }
      if (!ket_active.empty())
        W2[p].cols(ket_active) +=
            y.cols(ket_active) * reference.Lambda2.GetMatrix(ch_ket, ch_ket)(ket_active, ket_active);
    }
    std::vector<std::pair<size_t, size_t>> scalar_pairs;
    std::vector<arma::mat> W3;
    for (size_t ch = 0; ch < modelspace.GetNumberTwoBodyChannels(); ++ch)
    {
      const arma::uvec &active = ctx.active_kets[ch];
      if (active.empty())
        continue;
      const arma::mat &x = X.TwoBody.GetMatrix(ch, ch);
      const arma::mat &l = reference.Lambda2.GetMatrix(ch, ch);
      arma::mat w(x.n_rows, x.n_cols, arma::fill::zeros);
      w.cols(active) += x.cols(active) * l(active, active);
      w.rows(active) -= l(active, active) * x.rows(active);
      scalar_pairs.emplace_back(ch, ch);
      W3.push_back(std::move(w));
    }
    X.profiler.timer["MR comm221 lambda2 tensor products"] += omp_get_wtime() - section_start;

    section_start = omp_get_wtime();
    arma::mat result = ContractScalarTwoTensorOne(ctx, PartialTrace(ctx, channel_pairs, W2, L));
    result -= ContractTensorTwoScalarOne(ctx, PartialTrace(ctx, scalar_pairs, W3, 0));
    arma::mat first = PartialTrace(ctx, channel_pairs, P, L);
    X.profiler.timer["MR comm221 lambda2 tensor traces"] += omp_get_wtime() - section_start;

    section_start = omp_get_wtime();
    first += RingTerm(ctx);
    X.profiler.timer["MR comm221 lambda2 tensor ring"] += omp_get_wtime() - section_start;

    // PT(Y l X) and R(Y,X) are the Hermitian conjugates of PT(X l Y) and R(X,Y).
    const int hZ = -Hermiticity(X) * Hermiticity(Y);
    for (index_t i : ctx.orbits)
      for (index_t j : ctx.orbits)
        if (ctx.allowed[i][j])
          result(i, j) += first(i, j) +
                          hZ * modelspace.phase((modelspace.GetOrbit(i).j2 - modelspace.GetOrbit(j).j2) / 2) *
                              first(j, i);
    X.profiler.timer["MR comm221 lambda2 tensor total"] += omp_get_wtime() - total_start;
    return result;
  }

  Operator Commutator(const Operator &X, const Operator &Y,
                      const MRReference &reference)
  {
    // The lambda2 addons implement MR-IMSRG(2) only.  With a three-body input
    // or the global IMSRG(3) switch, the SR commutator would return 3B parts
    // without their MR terms, so refuse instead of mixing truncations.
    if (X.GetParticleRank() > 2 || Y.GetParticleRank() > 2 ||
        ::Commutator::use_imsrg3)
      throw std::invalid_argument(
          "MR commutator supports MR-IMSRG(2) only: 3B operators and IMSRG(3) are not implemented");
    Operator result = ::Commutator::Commutator(X, Y);
    // This exact branch is the SR-degeneration gate: no duplicated SR formula
    // and no floating-point add/subtract of nominally zero MR corrections.
    if (reference.Lambda2.Norm() == 0.0)
      return result;

    const double addon_start = omp_get_wtime();
    if (X.GetJRank() != 0 || Y.GetJRank() != 0 || X.GetTRank() != 0 ||
        Y.GetTRank() != 0 || X.GetParity() != 0 || Y.GetParity() != 0)
    {
      // Tensor output: lambda2 enters only the one-body part (the zero-body
      // contraction exists for scalars only).
      result.OneBody += comm221_lambda2_tensor(X, Y, reference);
      X.profiler.timer["MR lambda2 commutator addons"] +=
          omp_get_wtime() - addon_start;
      return result;
    }
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
