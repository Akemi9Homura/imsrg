# Repository Notes

## `fmt2=no2bpack`

This repository can read the packed binary NO2B interaction format produced by the `normal-order` code.

Use it by setting:

```bash
fmt2=no2bpack 2bme=<normal-order-output.bin> 3bme=none
```

The reader is `ReadWrite::Read_no2bpack()`. It reads the binary layout written by `normal-order`'s `Write_minipack()` implementation:

- oscillator frequency and `emax`
- orbit table in `(n, l, 2j, 2tz)`
- zero-body term
- upper-triangular one-body matrix elements
- packed J-coupled two-body matrix elements
- optional center-of-mass TBME payloads, which are consumed and ignored

The reader remaps file orbit indices into the active `ModelSpace` using `(n, l, 2j, 2tz)`, so it does not depend on matching raw orbit numbering between executables.

Because `no2bpack` files already contain the normal-ordered Hamiltonian pieces from `normal-order`, the main IMSRG driver does not add another `Trel_Op` for this format.
