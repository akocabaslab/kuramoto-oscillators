# Kuramoto Oscillators

Simulation code for coupled phase-oscillator models of active biofilm carpets. Each paper's code is listed below. The models differ, so use the code that belongs to the paper you are reproducing.

## 1. Synchronization and symmetry breaking of bacterial pili (eLife, 2026)

`kuramoto.m` is a MATLAB simulation of a 100 × 100 lattice of phase oscillators with no-flux boundaries,

```math
\frac{d\theta_i}{dt} = \omega - b \sin\theta_i + K \sum_{j \in N(i)} \sin\left(\theta_j - \theta_i - \alpha\right).
```

The parameters are set at the top of the script. To run it, open MATLAB in the repository folder and run:

```matlab
kuramoto
```

If you use this code, please cite:

Altın B, Günay ET, Yaman YI, Ünlü A, Gediz Y, Gedik N, Karadaş B, Başaran M, Kocabaş C, Özdemir Ş, Kocabaş A. *Controlling the synchronization and symmetry breaking of coupled bacterial pili on active biofilm carpets*. eLife. 2026. [doi:10.7554/eLife.107609](https://doi.org/10.7554/eLife.107609.2)

## 2. Self-organized travelling waves and mechanical computation (submitted)

The [`travelling-waves`](travelling-waves) folder contains Python code for the Kuramoto–Sakaguchi lattice of Eq. (1) of the manuscript *Self-Organized Travelling Waves and Emergent Mechanical Computation in Bacterial Active Matter*,

```math
\frac{d\phi_i}{dt} = \omega_0 - K \sum_{j \in N(i)} \sin\left(\phi_i - \phi_j - \alpha\right).
```

Installation, a demo with its expected output, and usage are described in [`travelling-waves/README.md`](travelling-waves/README.md).

The two papers write the phase lag with opposite signs: Eq. (1) is the coupling of `kuramoto.m` with α replaced by −α and b = 0.

## License

This repository is distributed under the MIT License. See `LICENSE` for details.
