"""Plotting layer -- separate from computation (headless-safe).

Every function takes a result object and returns a matplotlib Figure, so
calculations can run without creating figures.
"""
from __future__ import annotations

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:      # matplotlib optional
    plt = None


def _need_mpl():
    if plt is None:
        raise ImportError("matplotlib is required for plotting")


def plot_ampa(result, cmap="viridis"):
    """Concentration, open-state and desensitization maps for an AmpaResult."""
    _need_mpl()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ext = [0, result.time_ms[-1], result.distance_um[0], result.distance_um[-1]]

    with np.errstate(divide="ignore"):
        logc = np.log(np.where(result.concentration_uM > 0, result.concentration_uM, np.nan))
    im0 = axes[0].imshow(logc, aspect="auto", origin="lower", cmap=cmap,
                         extent=[0, result.glu_time_ms[-1],
                                 result.distance_um[0], result.distance_um[-1]])
    axes[0].set_title("log Glu concentration (uM)")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(result.desensitization, aspect="auto", origin="lower",
                         cmap=cmap, extent=ext)
    axes[1].set_title(f"{result.model_type} desensitization")
    fig.colorbar(im1, ax=axes[1])

    im2 = axes[2].imshow(result.open_state, aspect="auto", origin="lower",
                         cmap=cmap, extent=ext, vmin=0, vmax=result.open_state.max() or 1)
    axes[2].set_title(f"{result.model_type} open state")
    fig.colorbar(im2, ax=axes[2])

    for ax in axes:
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Distance (um)")
    fig.tight_layout()
    return fig


def plot_nmda(result, cmap="jet"):
    _need_mpl()
    fig, ax = plt.subplots(figsize=(7, 5))
    ext = [0, result.time_ms[-1], result.distance_um[0], result.distance_um[-1]]
    im = ax.imshow(result.open_state, aspect="auto", origin="lower", cmap=cmap, extent=ext)
    ax.set_title("NMDA open-state probability")
    ax.set_xlabel("Time (ms)"); ax.set_ylabel("Distance (um)")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    return fig


def plot_spillover(result):
    _need_mpl()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(result.time_ms, result.fraction_spillover, label="all")
    ax.plot(result.time_ms, result.fraction_free_spillover, label="free")
    ax.set_xlabel("Time (ms)"); ax.set_ylabel("Spillover fraction")
    ax.legend(); ax.grid(True)
    fig.tight_layout()
    return fig


def plot_particles(snapshot, lim=2000):
    """3D cloud from a snapshot array [X;Y;Z;State] (State 1=free)."""
    _need_mpl()
    X, Y, Z, S = snapshot
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    free = S == 1
    ax.scatter(X[~free], Y[~free], Z[~free], s=6, alpha=0.15, c="steelblue")
    ax.scatter(X[free], Y[free], Z[free], s=18, alpha=0.9, c="crimson")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_xlabel("X (nm)"); ax.set_ylabel("Y (nm)"); ax.set_zlabel("Z (nm)")
    fig.tight_layout()
    return fig
