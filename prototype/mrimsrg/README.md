# Rapid MR-IMSRG prototype

This directory contains the deliberately small m-scheme prototype described in
`docs/MR-IMSRG-快速结果计划.md`.  The C++ bridge does not implement another
interaction reader or NCSM solver: it builds against the validated
`/home/mengziyan/shell-model-obs` `Hamiltonian` and `simpleFCI` libraries.

Build the input bridge:

```bash
cmake -S prototype/mrimsrg -B prototype/mrimsrg/build \
  -DSHELL_MODEL_OBS_ROOT=/home/mengziyan/shell-model-obs \
  -DCMAKE_BUILD_TYPE=Release
cmake --build prototype/mrimsrg/build --target mrimsrg_prepare mrimsrg_validate -j2
```

Prepare and check the first He4 reference:

```bash
prototype/mrimsrg/build/mrimsrg_prepare \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --output prototype/mrimsrg/data/He4_Nrefmax0 --Z 2 --N 2 --nrefmax 0
python3 prototype/mrimsrg/densities.py prototype/mrimsrg/data/He4_Nrefmax0
```

Generated Hamiltonians, wavefunctions and densities belong under `data/` and
are ignored by Git.  The bridge refuses to overwrite an existing bundle.

Run the dependency-free unit tests with:

```bash
python3 -m unittest discover -s prototype/mrimsrg/tests -v
```

`commutator.py` implements the natural-orbital MR-IMSRG(2) commutator with
`lambda3=0`.  Its equations were regenerated independently with QCombo 0.2.0:
for a random four-orbit correlated reference, direct summation of the QCombo
0B/1B/2B expressions agreed with the numerical implementation to
`1.8e-15`.  The committed tests additionally compare the correlated 1B--2B
zero-body term and the single-reference single/double-excitation blocks with
explicit Fock-space matrix commutators. Direct calls reject a non-diagonal
`gamma1` instead of silently treating its diagonal as natural occupations.
For the required `He4/O16, Nrefmax=2` flows, `flow.py` diagonalizes only the
connected off-diagonal blocks of `gamma1` and evaluates the published equations
and `Delta e` mask in that temporary spherical natural-orbital basis, as stated
in Vobig Sec. 6.5.3.  The natural orbitals inherit the `e` labels of their
output slots. Final and checkpoint Hamiltonians are transformed back to the
unchanged input HO orbit ordering; the exact
orthogonal matrix is saved as `formula_basis_vectors.npy`. This is an internal
basis-covariant evaluation, not a natural-orbital NCSM optimization or a
replacement of the HO quantum-number labels.

`generator.py` implements one generator only: the modified White generator
used in the IM-NCSM literature, restricted by the relaxed `Delta e != 0`
one- and two-body masks. It follows the production `White-NCSM` truncation
defined in Vobig Sec. 6.5.4: all irreducible-density terms are omitted from
the generator numerator and the displayed `O(lambda2)` terms are omitted
from the Epstein--Nesbet denominators. The MR-IMSRG(2) commutator itself still
retains `lambda2`. A sign-preserving `1e-6 MeV` denominator cutoff is used;
in the Slater limit the denominators reduce term by term to
`src/Generator.cc`. Independently, the strict directional
`D1/D2=<Psi|H:A:|Psi>` diagnostic retains the published terms linear in
`lambda2` with `lambda3=0`. Its anti-Hermitian combination is monitored as
the strict diagnostic and agrees with `<Psi|[H,:A:]|Psi>` from the QCombo and
explicit-commutator checks.

`flow.py` directly integrates `dH/ds=[eta,H]` with an adaptive DOP853
stepper.  Like the existing C++ `IMSRGSolver`, it updates `eta` from every
trial Hamiltonian and stops when the norm of the actual masked White-NCSM
generator in Vobig Eqs. (6.5.28)--(6.5.29) falls to the requested fraction of
its initial value. The ODE state stores only independent
Hermitian one-body and antisymmetric/Hermitian pair-space two-body elements;
this mirrors the established operator storage and prevents redundant tensor
components from developing symmetry-violating numerical modes. The default
`max_step=10` is only an upper bound: DOP853 continues to reduce trial steps
to satisfy `rtol/atol`.

Every trajectory point stores three separately labeled quantities:
`generator_residual_ratio` is the norm of the denominator-weighted,
anti-Hermitian generator and is the formal acceptance condition;
`generator_numerator_residual_ratio` is the unweighted lambda-free numerator;
and `residual_ratio` is the strict lambda2-dependent `D-D^dagger` diagnostic.
The latter two are diagnostics and need not vanish at the approximate
generator's fixed point. This matches the `Eta.Norm()` convergence criterion
of the existing C++ `IMSRGSolver` and avoids claiming that the stricter
quantity vanished when it did not.

Run and materialize one calculation with:

```bash
PYTHONPATH=prototype/mrimsrg python3 prototype/mrimsrg/run_mrimsrg.py \
  prototype/mrimsrg/data/He4_Nrefmax0_final \
  prototype/mrimsrg/data/He4_flow_rtol1e-6 \
  --checkpoint-s 1.0
```

The output directory contains `gamma1/gamma2/lambda2`, initial and final
MR-normal-ordered tensors, and the final ordinary vacuum-normal-ordered
`E0+t+V` tensors.  `metadata.json` records the fixed interaction identity,
reference metadata, `lambda3=0`, generator/mask, ODE settings and every
accepted-step residual.  Existing output directories are never overwritten.
The output also materializes the ordinary vacuum Hamiltonian at `s=0` and at
the requested intermediate `--checkpoint-s` under `checkpoints/`; the root
payload is the final point.  Thus a production run retains the three points
required by the acceptance plan without keeping every large ODE state.
The command returns status 2 after saving if `smax` is reached before the
White-NCSM generator target, allowing a Slurm job to distinguish a diagnostic from an
accepted result.

Read a materialized ordinary Hamiltonian back into the existing NCSM solver:

```bash
prototype/mrimsrg/build/mrimsrg_validate \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --flow-output prototype/mrimsrg/data/He4_flow_rtol1e-6 \
  --Z 2 --N 2 --nmax 8 --states 3
```

The validator uses the original interaction only to initialize the identical
orbit/channel ordering, replaces its m-scheme matrix elements with the saved
vacuum `E0+t+V`, and then calls the existing `simpleFCI`.  At `s=0` this path
reproduces the complete 3060-dimensional He4 result as
`-20.3388325043 MeV, J=0`. By default the validator requests the three
lowest states (or the complete space when its dimension is smaller) and
prints each energy, excitation energy, and `2J`; `--states` changes that
small readback list without introducing another solver.

The bridge requires NumPy and SciPy on the Python side.  The first checked
integration run used `shell-model-obs` revision `1687f16` and reproduced the
full He4 `emax=2` ground-state benchmark as `-20.3388325043 MeV`; contraction
of the independently constructed `gamma1/gamma2` gave the same value.
The dense m-scheme readback interface is revision `4b10afa` of that dependency.

On point7, generate and inspect exactly one Slurm job at a time with:

```bash
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 --submit
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 --nrefmax 2
```

The dedicated generator is for this Python prototype; it deliberately does
not alter the existing production `gen_job.py` used by `imsrg++`.  It fixes the
documented point7 paths, one node/task, 64 cores, QOS/log conventions, sources
`sourceme.sh`, and prints the complete generated script before an optional
submission. Its single-point defaults are the accepted emax2 production window
`smax=25000`, `checkpoint_s=40`, `rtol=1e-6`, and `atol=1e-8`; the flow stops
early when its generator target is reached.  The longer upper bound leaves
room for the successive `2e-4`-scale occupation modes measured in the He4
correlated-reference case; it is only a bound, not a request to integrate past
convergence.
The accepted-step safety cap is 8000; neither upper bound changes the ODE when
the generator criterion is reached earlier. Use explicit
smaller values only for a labeled diagnostic. Reinvoke it for one different
nucleus or tolerance at a time;
there is no scan or batch-generation mode.
Use an explicit single-run label such as `--label smax5` when preserving a
short diagnostic and generating a longer flow with otherwise identical
tolerances; an existing materialized output is still never overwritten.

Reference metadata records the original development-machine interaction path.
When the exact fixed file is relocated, pass `--interaction <path>`; the runner
does not trust the path alone and still requires the frozen SHA-256 before it
starts a flow.  The point7 generator uses the task-specific copy under
`/tns/mengziyan/mr-imsrg-inputs/` because the shared force tree has larger
N2LOopt spaces but not this exact emax2 file.
The production path is compatible with point7's system Python 3.9 and the
generated batch script initializes Environment Modules explicitly before it
sources the machine-aware repository environment.
Compute nodes do not carry the login node's `/usr/local` NumPy/SciPy install,
so this standalone prototype calls the shared
`/opt/library/miniconda-3.12.9/bin/python3` interpreter explicitly and logs its
versions.  It does not activate or load miniconda and does not import the
system-Python pyIMSRG module.
