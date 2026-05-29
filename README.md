# Conservation Spectral SDK — Common Lisp

**Symbolic Theorem Proving for Spectral Graph Conservation**

The Common Lisp implementation of the Conservation Spectral framework — and the **only implementation that can do algebra on conservation**. Numerical computation, symbolic manipulation, and automated theorem proving all in one system. Seven source files covering graphs, matrices, symbolic algebra, a prover, debugging, and testing.

## The Aha Moment

Lisp treats code as data — S-expressions are both valid programs and valid algebraic expressions. `(+ (* x y) 3)` is simultaneously a Lisp form and a polynomial. This means we can write a **theorem prover** that manipulates conservation expressions symbolically, simplifies them, substitutes values, differentiates, and verifies properties like "conservation is bounded by 1 for all attributes on this graph." No other language in the SDK can prove things — they can only compute numbers. Lisp can compute *truths*. The `prover.lisp` module can verify conservation boundedness, monotonicity, edge sensitivity, and even universal quantification over all possible attributes.

## How to Use

```bash
# Install SBCL
sudo apt install sbcl

# Load and run
sbcl --load conservation.lisp --eval '(in-package :conservation-spectral)' --eval '(full-analysis)' --quit
```

### Run the test suite

```bash
sbcl --load test.lisp --quit
```

### Symbolic analysis

```lisp
;; Build a symbolic Laplacian and analyze conservation algebraically
(symbolic-conservation '((a b 1.0) (b c 1.0) (a c 1.0)))
;; Returns symbolic expressions, not numbers
```

### Theorem proving

```lisp
;; Prove conservation is bounded for a specific graph
(prove-conservation-bounded '(0 1 2) '((0 1) (1 2) (0 2)))
;; => #S(CONSERVATION-PROOF-RESULT :VERIFIED T ...)
```

## Architecture

| File | Purpose |
|------|---------|
| `conservation.lisp` | Core: graph construction, Laplacian, eigendecomposition, conservation ratios, anomaly detection |
| `matrix.lisp` | Dense matrix library: construction, arithmetic, eigendecomposition via QR iteration |
| `symbolic.lisp` | Symbolic algebra engine: S-expressions as algebraic expressions, simplification, differentiation, LaTeX output |
| `prover.lisp` | Theorem prover: verify conservation boundedness, monotonicity, edge sensitivity, universal quantification |
| `debug.lisp` / `debug2.lisp` | Debugging utilities and step-by-step trace |
| `test.lisp` | Test suite |

## Key Exports

- **Numerical**: `build-laplacian`, `eigendecompose`, `power-iteration`, `conservation-ratio`, `spectral-gap`, `cheeger-constant`, `full-analysis`
- **Symbolic**: `symbolic-laplacian`, `symbolic-conservation`, `sym-simplify`, `sym->latex`
- **Prover**: `prove-conservation-bounded`, `prove-monotonicity`, `prove-edge-sensitivity`, `prove-conservation-for-all-attributes`
- **Graph**: `make-graph`, `graph-add-edge`, `graph-add-vertex`

## Connection to the Conservation Spectral Framework

The conservation ratio α(G,a) = (a^T L a) / (λ_max ‖a‖²) is normally computed as a number. Lisp lets us keep it as an algebraic expression — the numerator is a quadratic form in the attribute variables, the denominator involves eigenvalue expressions. We can simplify, differentiate with respect to edge weights, and prove properties that hold for all possible attribute vectors. This is the "brain" of the SDK.

## Related Repos

- [conservation-spectral-forth](https://github.com/SuperInstance/conservation-spectral-forth) — Stack = sheaf stalk (the other "structural" language)
- [conservation-spectral-ada](https://github.com/SuperInstance/conservation-spectral-ada) — Range types and pre/post conditions
- [conservation-spectral-pascal](https://github.com/SuperInstance/conservation-spectral-pascal) — SET type for graph operations
- [conservation-spectral-apl](https://github.com/SuperInstance/conservation-spectral-apl) — One-line vector thinking
- [conservation-spectral-v2](https://github.com/SuperInstance/conservation-spectral-v2) — Reference Python implementation

## License

MIT

Part of the [SuperInstance OpenConstruct](https://github.com/SuperInstance/OpenConstruct) ecosystem.
