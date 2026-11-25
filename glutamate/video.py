"""On-demand export of the diffusion Monte Carlo as a web-ready video.

Renders the 3D field -- glutamate particles diffusing among the (semi-
transparent) obstacle spheres -- with an optional slow camera rotation, and
encodes an H.264 MP4 whose pixel size and file size are both controllable, so
the clip drops straight into the website tiles (e.g. 540x540, <= ~1 MB,
muted/looping, web-compatible yuv420p + faststart).

Typical use
-----------
    from glutamate import default_config
    from glutamate.video import export_diffusion_video
    cfg = default_config()
    size_bytes = export_diffusion_video(cfg, "StochasticSynapse.mp4",
                                        n_frames=80, size_px=540,
                                        fps=20, target_mb=1.0)

`target_mb` sets the approximate output size: the encoder bitrate is derived
from it (bitrate = 0.92 * target_MB * 8e6 / duration), so a smaller target
gives a smaller file. Pass target_mb=None for a quality-driven (larger) clip.
Use a ".gif" extension for a GIF fallback (larger, but needs no H.264).
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

from .diffusion import run_diffusion


def _collect_frames(cfg, seed, n_frames, progress_cb=None):
    """Run the diffusion once, capturing about ``n_frames`` field snapshots."""
    dt = cfg.computation.GTimeStep
    n_steps = int(round(cfg.basic.TotalTime / dt))
    every = max(1, n_steps // max(1, n_frames))
    frames, hold = [], {}

    def cb(fr):
        if "balls" in fr:
            hold["balls"] = fr["balls"]
        frames.append(fr)

    run_diffusion(cfg, rng=np.random.default_rng(seed),
                  frame_cb=cb, graphics_every=every,
                  progress_cb=(lambda p: progress_cb(0.7 * p)) if progress_cb else None)
    return frames, hold.get("balls")


_FG = "black"  # match the black axis labels of the other site figures
# MATLAB-style sans-serif so the overlays match the ASTRO/BRAINCELL/ARACHNE
# tiles. Arial is used on Windows; Liberation Sans is the Arial look-alike on
# Linux; DejaVu Sans is the final matplotlib fallback.
_FONT = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]

# legend marker colours -- same hue as the obstacles, a bit more opaque so a
# single marker is visible (a lone obstacle at its real alpha is invisible).
_ASTRO_RGBA = (1.0, 0.0, 0.0, 0.45)          # red = astroglia (adhesive)
_NEURON_RGBA = (0.254, 0.411, 0.882, 0.45)   # royalblue = neuron (non-adhesive)


def _synapse_surfaces(ax, R, gap, color=(0.5, 0.5, 0.5), alpha=0.32):
    """Two grey, semi-transparent hemispheres (pre- and post-synaptic) framing
    the cleft, oriented along z (cleft in the xy-plane). White wire-edges make
    the domes read against the crowded obstacle cloud."""
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi / 2, 20)
    U, V = np.meshgrid(u, v)
    xs = R * np.sin(V) * np.cos(U)
    ys = R * np.sin(V) * np.sin(U)
    zc = R * np.cos(V)               # R at pole -> 0 at the flat (cleft) face
    h = gap / 2.0
    for zsign in (1.0, -1.0):        # +1 presynaptic (up), -1 postsynaptic (down)
        a = 0.62 if zsign > 0 else alpha   # top (presynaptic) dome less transparent
        ax.plot_surface(xs, ys, zsign * (h + zc), color=color, alpha=a,
                        linewidth=0.4, edgecolors=(1, 1, 1, 0.35),
                        antialiased=True, shade=True)
    # bright white cut surfaces (the two membrane faces bounding the cleft);
    # the gap is exaggerated (not to scale) so the thin cleft reads clearly.
    rad = np.linspace(0, R, 6)
    ang = np.linspace(0, 2 * np.pi, 44)
    Rr, Ang = np.meshgrid(rad, ang)
    dx = Rr * np.cos(Ang); dy = Rr * np.sin(Ang)
    for zz in (h, -h):
        ax.plot_surface(dx, dy, np.full_like(dx, zz), color="white", alpha=0.9,
                        linewidth=0, antialiased=True, shade=False)


def _add_legend(ax, fontsize):
    """Static bottom legend: astroglia (red) and neuron (blue) marker spheres."""
    y = 0.055
    ax.text2D(0.20, y, "●", transform=ax.transAxes, color=_ASTRO_RGBA,
              fontsize=fontsize * 2.6, ha="center", va="center")
    ax.text2D(0.27, y, "astroglia", transform=ax.transAxes, color=_FG,
              fontfamily=_FONT, fontsize=fontsize * 0.9, ha="left", va="center")
    ax.text2D(0.60, y, "●", transform=ax.transAxes, color=_NEURON_RGBA,
              fontsize=fontsize * 2.6, ha="center", va="center")
    ax.text2D(0.67, y, "neuron", transform=ax.transAxes, color=_FG,
              fontfamily=_FONT, fontsize=fontsize * 0.9, ha="left", va="center")


def _render(fig, canvas, cfg, fr, balls, elev, azim, ball_mode="both",
            scale_um=1.0, show_time=True, fontsize=13.0,
            show_synapse=True, synapse_radius=300.0, synapse_gap=80.0,
            show_legend=True):
    """Draw one 3D field frame and return it as an HxWx3 uint8 array."""
    fig.clear()
    ax = fig.add_subplot(111, projection="3d")
    lim = cfg.cube.AxisSize
    try:
        ax.set_proj_type("ortho")  # true scale, no perspective foreshortening
    except Exception:
        pass
    if balls is not None and balls["xyz"].shape[0] > 0:
        xyz, r, adh = balls["xyz"], balls["radius"], balls["adhesive"]
        # marker area ~ physical diameter relative to the +/-lim cube
        ppn = (fig.get_size_inches()[0] * 72.0 * 0.55) / (2.0 * lim)
        s = (2.0 * r * ppn) ** 2
        if ball_mode in ("both", "adhesive"):
            ax.scatter(xyz[adh, 0], xyz[adh, 1], xyz[adh, 2], s=s[adh],
                       c="red", alpha=0.10, edgecolors="none")
        if ball_mode in ("both", "non-adhesive"):
            ax.scatter(xyz[~adh, 0], xyz[~adh, 1], xyz[~adh, 2], s=s[~adh],
                       c="royalblue", alpha=0.07, edgecolors="none")
    if show_synapse:
        _synapse_surfaces(ax, synapse_radius, synapse_gap)
    x, y, z = fr["x"], fr["y"], fr["z"]
    act, bnd = fr["active"], fr["bound"]
    free = act & ~bnd
    ax.scatter(x[free], y[free], z[free], s=6, c="deepskyblue",
               alpha=0.9, edgecolors="none")
    ax.scatter(x[act & bnd], y[act & bnd], z[act & bnd], s=12, c="crimson",
               alpha=0.9, edgecolors="none")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass

    # --- spatial scale bar: a vertical (z) world ruler. Because the camera
    # pans about the vertical axis at fixed elevation, a z-aligned bar keeps a
    # constant on-screen length, so it reads as an honest micron scale. ---
    if scale_um and scale_um > 0:
        L = scale_um * 1000.0  # nm
        x0, y0 = lim * 0.86, lim * 0.86
        z0 = -lim * 0.9
        ax.plot([x0, x0], [y0, y0], [z0, z0 + L], c=_FG, lw=2, solid_capstyle="butt")
        ax.text(x0, y0, z0 + L / 2.0, "  %g µm" % scale_um, color=_FG,
                fontfamily=_FONT, fontsize=fontsize * 0.85, ha="left", va="center")

    # --- time counter (screen overlay, top-left) ---
    if show_time:
        ax.text2D(0.04, 0.95, "t = %.2f ms" % float(fr.get("t", 0.0)),
                  transform=ax.transAxes, color=_FG, fontfamily=_FONT,
                  fontsize=fontsize, ha="left", va="top")

    if show_legend:
        _add_legend(ax, fontsize)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[..., :3].copy()


def encode_frames(imgs, path, fps=20, target_mb=1.0):
    """Encode a list of RGB frames to ``path`` (.mp4/.webm/.gif). Returns bytes."""
    try:
        import imageio
    except ImportError:
        raise ImportError(
            "Video export needs the optional imageio package. Install it with:\n"
            "    pip install imageio imageio-ffmpeg\n"
            "(or: pip install -e \".[video]\"). Everything else works without it."
        ) from None
    ext = os.path.splitext(path)[1].lower()
    dur = max(1e-3, len(imgs) / float(fps))
    if ext == ".gif":
        imageio.mimwrite(path, imgs, format="GIF", fps=fps, loop=0)
        return os.path.getsize(path)
    kw = dict(fps=fps, codec="libx264", quality=None, macro_block_size=1,
              pixelformat="yuv420p", ffmpeg_params=["-movflags", "+faststart"])
    if target_mb:
        kw["bitrate"] = int(0.92 * target_mb * 8e6 / dur)  # ~size headroom
    writer = imageio.get_writer(path, **kw)
    for im in imgs:
        writer.append_data(im)
    writer.close()
    return os.path.getsize(path)


def export_diffusion_video(cfg, path, n_frames=80, size_px=540, fps=20,
                           target_mb=1.0, seed=None, rotate=True, elev=18.0,
                           azim0=45.0, azim_span=90.0, ball_mode="both",
                           scale_um=1.0, show_time=True, font_scale=1.5,
                           show_synapse=True, synapse_radius=300.0,
                           synapse_gap=80.0, show_legend=True, progress_cb=None):
    """Run the diffusion and write a web-ready video. Returns the file size (bytes).

    Parameters
    ----------
    cfg        : Config for the simulation.
    path       : output file (.mp4 recommended; .webm or .gif also work).
    n_frames   : number of frames captured across the run (site clips ~80).
    size_px    : square output resolution in pixels (site tiles use 540).
    fps        : frames per second of the clip (loop length = n_frames / fps).
    target_mb  : approximate target file size in MB (None -> quality-driven).
    rotate     : slowly pan the camera azimuth from azim0 by azim_span degrees.
    ball_mode  : 'both' | 'adhesive' | 'non-adhesive' obstacles to show.
    progress_cb: optional callback(percent 0..100).
    """
    frames, balls = _collect_frames(cfg, seed, n_frames, progress_cb)
    if not frames:
        raise RuntimeError("No frames were produced by the diffusion run.")
    dpi = 100.0
    inch = size_px / dpi
    fig = Figure(figsize=(inch, inch), dpi=dpi)
    canvas = FigureCanvasAgg(fig)
    fontsize = max(8.0, size_px / 40.0) * font_scale
    imgs = []
    N = len(frames)
    for k, fr in enumerate(frames):
        azim = azim0 + (azim_span * k / max(1, N - 1) if rotate else 0.0)
        imgs.append(_render(fig, canvas, cfg, fr, balls, elev, azim, ball_mode,
                            scale_um=scale_um, show_time=show_time, fontsize=fontsize,
                            show_synapse=show_synapse, synapse_radius=synapse_radius,
                            synapse_gap=synapse_gap, show_legend=show_legend))
        if progress_cb:
            progress_cb(70.0 + 30.0 * (k + 1) / N)
    return encode_frames(imgs, path, fps=fps, target_mb=target_mb)
