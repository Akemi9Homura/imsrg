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

  /// Full MR-IMSRG(2) commutator for J_ref=0 and lambda3=0.
  /// Lambda2=0 returns the existing scalar IMSRG(2) commutator directly.
  Operator Commutator(const Operator &X, const Operator &Y,
                      const MRReference &reference);
}

#endif
