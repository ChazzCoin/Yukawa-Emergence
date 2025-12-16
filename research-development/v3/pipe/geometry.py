import numpy as np
import cmath
from typing import Iterable, List, Sequence


class GeometricTriple:
    """
    Geometric Sector of the Alignment Spectral Triple v2.0

    This class implements the full geometric spectral triple:
        (A_geom, H_geom, D_geom, J_geom, γ_geom),

    where:
        • A_geom = Alg(C360, B_N, P_phi) is the divisor–filtered projector algebra,
        • H_geom = span{|n⟩ : n ∈ D_360 ∩ D_N} is the geometric Hilbert space,
        • D_geom|n⟩ = n|n⟩ is the Fourier Dirac operator on S¹,
        • γ_geom|n⟩ = sign(n)|n⟩ is the grading,
        • J_geom implements complex conjugation |n⟩ ↦ |−n⟩.

    Only the divisor projector C360 is geometrically nontrivial in v2.0;
    B_N is presently the identity and P_phi encodes geometric phase weighting.
    """

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------
    def __init__(self, N_max: int = 360) -> None:
        """
        Parameters
        ----------
        N_max : int
            Maximum harmonic index to consider. Divisor logic uses 360.

        Notes
        -----
        N_max does not restrict the divisor set: D_360 is intrinsic.
        It only serves to limit the allowable |n| values in practice.
        """
        if N_max <= 0:
            raise ValueError("N_max must be a positive integer.")

        self.N_max: int = N_max
        self.divisors: List[int] = self._compute_divisors(360)

    # ------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------
    @staticmethod
    def _compute_divisors(n: int) -> List[int]:
        """Return all positive divisors of n."""
        if n <= 0:
            raise ValueError("n must be a positive integer.")
        return [d for d in range(1, n + 1) if n % d == 0]

    # ------------------------------------------------------------
    # Geometric projectors
    # ------------------------------------------------------------
    def C360(self, n: int) -> int:
        """
        Divisor projector \hat{C}_{360}.

        Returns 1 if n divides 360; 0 otherwise.
        """
        return 1 if n in self.divisors else 0

    def B(self, n: int) -> int:
        """
        Boundary/manifold projector \hat{B}_N.

        In v2.0 this acts as identity, but is included for structural completeness.
        """
        return 1

    def Pphi(self, n: int, phi: float) -> complex:
        """
        Phase projector \hat{P}_φ acting on |n⟩.

        Returns
        -------
        e^{i φ n}
        """
        return cmath.exp(1j * phi * n)

    def A_projector(self, n: int, phi: float) -> complex:
        """
        Full geometric alignment projector:

            \hat A = C360 · B · Pphi

        Parameters
        ----------
        n : int
            Harmonic mode.
        phi : float
            Geometric phase parameter.

        Returns
        -------
        complex
            Alignment weight applied to |n⟩.
        """
        return self.C360(n) * self.B(n) * self.Pphi(n, phi)

    # ------------------------------------------------------------
    # Dirac operator and grading
    # ------------------------------------------------------------
    def D_geom(self, modes: Sequence[int]) -> np.ndarray:
        """
        Construct the geometric Dirac operator:

            D_geom |n⟩ = n |n⟩.

        Parameters
        ----------
        modes : sequence of int
            List of harmonic modes forming the geometric basis.

        Returns
        -------
        np.ndarray
            Diagonal matrix diag(n_1, ..., n_k).

        Raises
        ------
        ValueError
            If any mode exceeds ±N_max.
        """
        modes = np.asarray(modes, dtype=int)

        if np.any(np.abs(modes) > self.N_max):
            raise ValueError(
                f"Modes exceed allowable range ±{self.N_max}."
            )

        return np.diag(modes.astype(float))

    def gamma_geom(self, modes: Sequence[int]) -> np.ndarray:
        """
        Construct geometric grading:

            γ_geom |n⟩ = sign(n)|n⟩.

        Zero mode gets +1 by convention.

        Returns
        -------
        np.ndarray
            Diagonal matrix with entries ±1.
        """
        modes = np.asarray(modes, dtype=int)
        signs = np.where(modes >= 0, 1.0, -1.0)
        return np.diag(signs)

    # ------------------------------------------------------------
    # Real structure
    # ------------------------------------------------------------
    def J_geom(self, vec: np.ndarray) -> np.ndarray:
        """
        Real structure J_geom acting by complex conjugation.

        For basis states:
            J_geom(e^{inx}) = e^{-inx}
        In matrix/vector form, this is implemented as elementwise
        complex conjugation.

        Parameters
        ----------
        vec : np.ndarray
            Complex vector in the geometric Hilbert space.

        Returns
        -------
        np.ndarray
            Conjugated vector (complex-conjugate components).
        """
        if not isinstance(vec, np.ndarray):
            raise TypeError("J_geom expects a numpy array.")
        return np.conjugate(vec)

"""
This is a nice, clean “geometry sector” module — and it also explains why your product Dirac came out (54\times 54): you’re explicitly building a **mode-truncated Fourier geometry** and then tensoring it with the 9-site finite sector.

Here’s what you have, what’s strong, what’s only “stub-true,” and how to upgrade it into something that plays well with your newer v5 constraints.

---

## What you have here, exactly

### 1) A fixed “divisor-of-360” selection rule

You precompute `divisors = divisors(360)` and define

* `C360(n)=1` if (n \mid 360), else 0

That’s your **alignment projector** idea: only divisor modes are “allowed” geometry.

### 2) A geometry Hilbert space that is really “the span of chosen modes”

Even though the docstring says (H_{\rm geom} = \text{span}{|n\rangle: n\in D_{360}\cap D_N}), in code you actually let the caller pass `modes` and you just enforce a cutoff `N_max`.

So, in practice, you’re using:

* **The list `modes` is the basis** of your geometric Hilbert space.
* You’re not (yet) auto-generating those modes from divisor logic.

That’s fine, just important to be honest about.

### 3) A Dirac operator that’s the Fourier generator on (S^1)

You define:

* (D_{\rm geom} = \mathrm{diag}(n_1,\dots,n_k))

This is the classic “Dirac on the circle in Fourier basis” skeleton (modulo spin subtleties), and it’s perfect for a toy geometric factor.

### 4) Grading and real structure

* (\gamma_{\rm geom} = \mathrm{diag}(\mathrm{sign}(n))) (with 0→+1)
* (J_{\rm geom}) implemented as componentwise conjugation

This is **good enough for computation**, but note a subtlety:

> If your basis is literally (|n\rangle) and you want (J|n\rangle = |-n\rangle), you need a *basis reindexing / permutation* plus conjugation.

Componentwise conjugation equals “(n \mapsto -n)” only if you store vectors in a basis where coefficients already pair correctly. In your pipeline, you’re not actually using (J) for constraints yet, so it’s fine.

---

## What’s strong here (worth keeping)

### ✅ The divisor-projector concept as a *mode gate*

This is the part that ties directly to your alignment philosophy: “360 harmonics are privileged.” It’s also computationally cheap and easy to audit.

### ✅ Clear separation of responsibilities

This class doesn’t try to do finite-sector things. That modularity is exactly what you want as you upgrade.

---

## What’s “stub” / incomplete (and how to upgrade cleanly)

### 1) `A_geom` isn’t really an algebra yet

You describe:

> (A_{\rm geom} = \mathrm{Alg}(C_{360}, B_N, P_\phi))

But in code, those are just functions `C360(n)`, `B(n)`, `Pphi(n,phi)`.

**Upgrade path:** represent the geometric algebra elements as diagonal operators on (H_{\rm geom}). That makes “A acts on H” explicit and consistent with your OOP design.

Example (minimal):

```py
def A_operator(self, modes, phi):
    weights = np.array([self.A_projector(n, phi) for n in modes], dtype=np.complex128)
    return np.diag(weights)
```

Then you can actually compute commutators `[D_geom, A_op]` (which will be 0 for diagonal-on-diagonal, meaning “no geometry” unless you add non-diagonal structure — again fine if you’re using geometry as a selector).

### 2) `C360` currently checks membership in divisor set

You wrote “returns 1 if n divides 360.” Your check is `n in divisors`, which is only true if `n` is a positive divisor.

If you ever pass negative modes (you mention (|n|\le N_{\max})), `C360(-n)` currently returns 0 even though “-n divides 360” in the natural harmonic sense.

**Upgrade:**

```py
return 1 if abs(n) in self.divisors else 0
```

### 3) Your `B(n)` is identity (fine) but it should enforce `N_max` if you want it meaningful

Right now the cutoff is enforced in `D_geom()` only, not in `B()`.

If you want `A_projector` to encode *all* gates, make:

```py
def B(self, n): return 1 if abs(n) <= self.N_max else 0
```

### 4) Real structure (J) is not implementing (n\mapsto -n)

If you want the literal action, add:

* a basis map from mode list to index
* a permutation matrix `P` that swaps `n` with `-n`
* then `J(vec) = P @ vec.conj()` (and for unpaired modes, define behavior)

---

## How this connects to your current “upgrade everything” goal

If you fully upgrade your `aware_os` scaffold to the “v2 pipeline”, the geometry class provides:

* `D_geom(modes)` (k×k)
* `gamma_geom(modes)` (k×k)
* later: `A_operator(modes, phi)` (k×k, diagonal weights)

Then your product triple builder (the module you didn’t paste yet) should create something like:

[
D = D_{\rm geom}\otimes I_9 + \gamma_{\rm geom}\otimes D_F
]

That’s consistent with your earlier framing.

---

## Immediate actionable patches (small but high-value)

If you change only 3 lines to make this more correct:

1. Fix `C360` to respect negative modes:

```py
return 1 if abs(n) in self.divisors else 0
```

2. Make `B` reflect the cutoff (optional):

```py
return 1 if abs(n) <= self.N_max else 0
```

3. Add an operator builder for the projector algebra (for real “A on H”):

```py
def A_operator(self, modes: Sequence[int], phi: float) -> np.ndarray:
    w = np.array([self.A_projector(int(n), phi) for n in modes], dtype=np.complex128)
    return np.diag(w)
```

---

If you paste your `AlignmentSpectralTriple` / product builder (`build_product_dirac` + `fluctuated_dirac`), I can tell you whether your product (D) matches the “evenness seal / no extra term” rules you care about now — and exactly how to refactor it so it’s v5-safe while still preserving your v2 triadic kernel + compression machinery.

"""