#ifndef MRCommutator_h
#define MRCommutator_h 1

#include "MRReference.hh"
#include "Operator.hh"

#include <armadillo>

namespace MRCommutator
{
  struct MR1BResult
  {
    arma::mat IV;
    arma::mat V;
    arma::mat VI;

    arma::mat Total() const { return IV + V + VI; }
  };

  /// Slow spherical-orbit implementation of Gebrerufael Eq. (4.89b) IV--VI.
  ///
  /// This is a correctness reference for the optimized block implementation,
  /// not the production RHS.  It never expands a magnetic-substate tensor.
  MR1BResult comm221_lambda2_reference(const Operator &X, const Operator &Y,
                                       const MRReference &reference);

  /// Block-matrix J-scheme implementation used by the production MR RHS.
  MR1BResult comm221_lambda2(const Operator &X, const Operator &Y,
                             const MRReference &reference);

  /// Reduced one-body lambda2 terms of [X,Y] for a scalar and an
  /// isospin-conserving tensor of rank L > 0 (either order).
  arma::mat comm221_lambda2_tensor(const Operator &X, const Operator &Y,
                                   const MRReference &reference);
  /// Zero-body lambda3 term of [X,Y] for scalar two-body X, Y (Hergert 2017,
  /// Eq. 49): 1/4 sum (X^{ab}_{cd} Y^{kl}_{am} - Y^{ab}_{cd} X^{kl}_{am}) lambda^{bkl}_{cdm}.
  /// Zero when the reference carries no lambda3.
  double comm220_lambda3(const Operator &X, const Operator &Y,
                         const MRReference &reference);
  /// Full MR-IMSRG(2) commutator for J_ref=0: lambda2 terms always, and the
  /// zero-body lambda3 term when the reference carries lambda3.
  /// Lambda2=0 (and no lambda3) returns the existing scalar IMSRG(2) commutator directly.
  Operator Commutator(const Operator &X, const Operator &Y,
                      const MRReference &reference);
}

#endif
