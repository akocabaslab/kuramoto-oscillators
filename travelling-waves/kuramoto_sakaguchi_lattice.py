"""
Minimal simulation of the two-dimensional Kuramoto-Sakaguchi lattice, Eq. (1):

    dphi_i/dt = omega0 - K * sum_{j in N(i)} sin(phi_i - phi_j - alpha)

N(i) are the four nearest neighbours of site i on an L x L lattice with
periodic boundaries. All oscillators share the natural frequency omega0;
K is the coupling strength and alpha the phase lag. Phases start uniformly
at random and are integrated with the explicit Euler method. The run stops
when the lattice is frequency-locked (spread of instantaneous frequencies
below `tol`) or when t_max is reached.

At the end the script prints the locked frequency Omega and the number of
topological defects, and saves snapshots of the phase field. As a consistency
check it also prints omega0 + K sin(alpha) * Cbar, where Cbar is the lattice
average of C_i = sum_j cos(phi_i - phi_j). On the periodic lattice this equals
the mean instantaneous frequency at every time, locked or not, so it checks the
implementation, not locking; locking is judged by the frequency spread.

Requirements: Python 3, NumPy, Matplotlib (tested with Python 3.13,
NumPy 2.5 and Matplotlib 3.11). A default run takes roughly 10 s.

Examples
    python kuramoto_sakaguchi_lattice.py                  # alpha = 0.15*pi, seed 42
    python kuramoto_sakaguchi_lattice.py --alpha 0        # reciprocal limit
    python kuramoto_sakaguchi_lattice.py --alpha 0.10 --L 64 --seed 1
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt

NEIGHBOURS = ((1, 0), (-1, 0), (1, 1), (-1, 1))   # (shift, axis) for np.roll


def rhs(phi, omega0, K, alpha):
    """Right-hand side of Eq. (1) for every site at once."""
    coupling = np.zeros_like(phi)
    for shift, axis in NEIGHBOURS:
        phi_j = np.roll(phi, shift, axis=axis)       # periodic nearest neighbour
        coupling += np.sin(phi - phi_j - alpha)
    return omega0 - K * coupling


def local_coherence(phi):
    """C_i = sum over the four neighbours of cos(phi_i - phi_j)."""
    return sum(np.cos(phi - np.roll(phi, s, axis=a)) for s, a in NEIGHBOURS)


def winding_numbers(phi):
    """Topological charge (+1, -1 or 0) of every plaquette."""
    def wrap(x):
        return (x + np.pi) % (2 * np.pi) - np.pi
    a = phi
    b = np.roll(phi, -1, axis=1)
    c = np.roll(b, -1, axis=0)
    d = np.roll(phi, -1, axis=0)
    total = wrap(b - a) + wrap(c - b) + wrap(d - c) + wrap(a - d)
    return np.rint(total / (2 * np.pi)).astype(int)


def simulate(L=32, K=0.5, omega0=1.0, alpha=0.15 * np.pi, dt=0.01,
             t_max=5000.0, tol=1e-5, seed=42, snapshot_times=(0.0, 10.0, 100.0)):
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=(L, L))
    snapshot_steps = {int(round(t / dt)) for t in snapshot_times}
    snapshots = {}
    n_max = int(round(t_max / dt))
    for n in range(n_max + 1):
        if n in snapshot_steps:
            snapshots[n * dt] = phi.copy()            # actual time of the snapshot
        dphi = rhs(phi, omega0, K, alpha)
        spread = dphi.max() - dphi.min()
        if spread < tol or n == n_max:               # phi, dphi and n*dt stay consistent
            break
        phi = phi + dt * dphi                         # explicit Euler step
    return phi, n * dt, dphi, spread, snapshots


def main():
    p = argparse.ArgumentParser(description="2D Kuramoto-Sakaguchi lattice, Eq. (1)")
    p.add_argument("--alpha", type=float, default=0.15, help="phase lag in units of pi")
    p.add_argument("--L", type=int, default=32, help="lattice size (L x L)")
    p.add_argument("--K", type=float, default=0.5, help="coupling strength")
    p.add_argument("--omega0", type=float, default=1.0, help="natural frequency")
    p.add_argument("--dt", type=float, default=0.01, help="Euler time step")
    p.add_argument("--tmax", type=float, default=5000.0, help="maximum time")
    p.add_argument("--tol", type=float, default=1e-5, help="locking tolerance")
    p.add_argument("--seed", type=int, default=42, help="random seed")
    p.add_argument("--out", default="kuramoto_lattice.png", help="output figure")
    p.add_argument("--no-show", action="store_true", help="do not open a window")
    args = p.parse_args()
    if args.L < 2:
        p.error("--L must be at least 2")
    if args.dt <= 0 or args.tmax < 0 or args.tol <= 0:
        p.error("--dt and --tol must be positive and --tmax non-negative")

    alpha = args.alpha * np.pi
    phi, t_end, dphi, spread, snapshots = simulate(
        L=args.L, K=args.K, omega0=args.omega0, alpha=alpha, dt=args.dt,
        t_max=args.tmax, tol=args.tol, seed=args.seed)

    locked = spread < args.tol
    omega = dphi.mean()
    cbar = local_coherence(phi).mean()
    w = winding_numbers(phi)
    print(f"alpha = {args.alpha:g} pi, L = {args.L}, K = {args.K}, seed = {args.seed}")
    if locked:
        print(f"locked at t = {t_end:.2f} (frequency spread {spread:.2e})")
        print(f"Omega = {omega:.4f}")
    else:
        print(f"NOT locked by t = {t_end:.2f} (frequency spread {spread:.2e})")
        print(f"mean instantaneous frequency = {omega:.4f}")
    print(f"consistency check: omega0 + K sin(alpha) Cbar = "
          f"{args.omega0 + args.K * np.sin(alpha) * cbar:.4f}")
    print(f"defects: {(w == 1).sum()} with charge +1, {(w == -1).sum()} with charge -1")

    snapshots[t_end] = phi
    times = sorted(snapshots)
    fig, axes = plt.subplots(1, len(times), figsize=(3.2 * len(times), 3.4),
                             constrained_layout=True)
    for ax, t in zip(np.atleast_1d(axes), times):
        field = np.angle(np.exp(1j * snapshots[t]))   # wrap to (-pi, pi]
        im = ax.imshow(field, cmap="twilight", vmin=-np.pi, vmax=np.pi,
                       origin="lower", interpolation="nearest")
        if t == t_end:
            for charge, colour in ((1, "black"), (-1, "white")):
                rows, cols = np.nonzero(winding_numbers(snapshots[t]) == charge)
                ax.scatter(cols + 0.5, rows + 0.5, s=40, facecolors="none",
                           edgecolors=colour, linewidths=1.5)
        ax.set_title(f"t = {t:g}")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=axes, shrink=0.8, label="phase (rad)")
    fig.suptitle(f"Kuramoto-Sakaguchi lattice, alpha = {args.alpha:g} pi "
                 f"(circles: black +1, white -1 defects)")
    fig.savefig(args.out, dpi=150)
    print(f"saved {args.out}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
