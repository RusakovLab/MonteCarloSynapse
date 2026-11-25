"""Monte Carlo glutamate diffusion (port of DiffusionGlutamateBalls.m).

Vectorized over particles. Statistically faithful, not bit-identical to MATLAB
(different RNG). Coordinates in nm, time in ms. An optional frame_cb enables a
live 3D field view (obstacles + free/bound molecules) during the run.

Fix vs the MATLAB original: a particle is treated as "inside the cleft" only
when it is within the cleft radius AND within the cleft half-width in z. The
MATLAB code tested the xy-radius only, which mis-classified particles that sit
above/below the synapse (radial_xy <= cleft but |z| large); the cleft z-reflection
then threw them far and the solid-volume check reverted the step, freezing them.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np

from .config import Config


@dataclass
class DiffusionResult:
    average_free: np.ndarray
    average_bound: np.ndarray
    time_ms: np.ndarray
    number_spillover: np.ndarray
    fraction_spillover: np.ndarray
    number_free_spillover: np.ndarray
    fraction_free_spillover: np.ndarray
    first_passage_time: np.ndarray
    first_passage_latency: np.ndarray
    release_time: np.ndarray
    occupied_fraction: float
    free_fraction: float
    occupancy_se: float
    mean_diffusion: float
    std_diffusion: float
    balls_xyz: np.ndarray
    balls_radius: np.ndarray
    balls_adhesive: np.ndarray
    metadata: Dict = field(default_factory=dict)
    snapshots: Dict[str, np.ndarray] = field(default_factory=dict)


def _assign_releases(np_particles, release_times, release_fractions, rng):
    raw = np_particles * np.asarray(release_fractions, float)
    counts = np.floor(raw).astype(int)
    remaining = np_particles - counts.sum()
    order = np.argsort(-(raw - counts))
    counts[order[:remaining]] += 1
    perm = rng.permutation(np_particles)
    release_time = np.full(np_particles, np.nan)
    start = 0
    for k, t in enumerate(release_times):
        idx = perm[start:start + counts[k]]
        release_time[idx] = t
        start += counts[k]
    return release_time, counts


def _place_balls(cfg, rng, first_ball_index, astro_center=None, astro_radius=None):
    B, cube, cleft = cfg.balls, cfg.cube, cfg.cleft
    n = cfg.basic.Numbersphere
    xd = np.zeros(n); yd = np.zeros(n); zd = np.zeros(n)
    radius = np.zeros(n); adhesive = np.zeros(n, bool)
    L = cube.Size
    for j in range(first_ball_index, n):
        placed = False
        for _ in range(B.MaximumPlacementAttempts):
            cx = L * (rng.random() - 0.5)
            cy = L * (rng.random() - 0.5)
            cz = L * (rng.random() - 0.5)
            rr = B.MinBallRadius + rng.random() * (B.MaxBallRadius - B.MinBallRadius)
            if cfg.basic.ReleaseSource.lower() == "synapse":
                ok = np.sqrt(cx * cx + cy * cy + cz * cz) - rr >= cleft.SynapticEnvelopeRadius
            else:
                ok = (np.sqrt((cx - astro_center[0]) ** 2 + (cy - astro_center[1]) ** 2
                              + (cz - astro_center[2]) ** 2) - rr
                      >= astro_radius + cleft.DistanceAstrocyte)
            if not ok:
                continue
            if B.PlacementMode.lower() == "overlapping":
                placed = True
            else:
                d = np.sqrt((cx - xd[:j]) ** 2 + (cy - yd[:j]) ** 2 + (cz - zd[:j]) ** 2)
                placed = bool(np.all(d - rr - radius[:j] >= B.MinimumBallSeparation))
            if placed:
                xd[j], yd[j], zd[j], radius[j] = cx, cy, cz, rr
                adhesive[j] = rng.random() < B.ProbabilityofAstrocytes
                break
        if not placed:
            raise RuntimeError("Unable to place ball %d" % j)
    return xd, yd, zd, radius, adhesive


def _occupancy(cfg, rng, xd, yd, zd, radius, first_ball_index, astro_center=None, astro_radius=None):
    B, cube, cleft = cfg.balls, cfg.cube, cfg.cleft
    n = B.OccupancySampleCount
    L = cube.Size
    bx = xd[first_ball_index:]; by = yd[first_ball_index:]; bz = zd[first_ball_index:]
    br2 = radius[first_ball_index:] ** 2
    px = np.empty(0); py = np.empty(0); pz = np.empty(0)
    while px.size < n:
        m = int((n - px.size) * 1.5) + 16
        cx = L * (rng.random(m) - 0.5)
        cy = L * (rng.random(m) - 0.5)
        cz = L * (rng.random(m) - 0.5)
        if cfg.basic.ReleaseSource.lower() == "synapse":
            keep = np.sqrt(cx * cx + cy * cy + cz * cz) > cleft.SynapticEnvelopeRadius
        else:
            keep = (np.sqrt((cx - astro_center[0]) ** 2 + (cy - astro_center[1]) ** 2
                            + (cz - astro_center[2]) ** 2) > astro_radius + cleft.DistanceAstrocyte)
        px = np.concatenate([px, cx[keep]])
        py = np.concatenate([py, cy[keep]])
        pz = np.concatenate([pz, cz[keep]])
    px = px[:n]; py = py[:n]; pz = pz[:n]
    inside = np.zeros(n, bool)
    chunk = max(1, int(2_000_000 / max(1, bx.size)))
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        d2 = ((px[s:e, None] - bx[None, :]) ** 2 + (py[s:e, None] - by[None, :]) ** 2
              + (pz[s:e, None] - bz[None, :]) ** 2)
        inside[s:e] = np.any(d2 < br2[None, :], axis=1)
    betta = inside.mean()
    se = np.sqrt(betta * (1 - betta) / n)
    return float(betta), float(se)


def _containment(px, py, pz, bx, by, bz, br2, adhesive_b):
    m = px.size
    inside = np.zeros(m, bool); inside_adh = np.zeros(m, bool)
    if m == 0 or bx.size == 0:
        return inside, inside_adh
    chunk = max(1, int(2_000_000 / bx.size))
    for s in range(0, m, chunk):
        e = min(m, s + chunk)
        d2 = ((px[s:e, None] - bx[None, :]) ** 2 + (py[s:e, None] - by[None, :]) ** 2
              + (pz[s:e, None] - bz[None, :]) ** 2)
        within = d2 < br2[None, :]
        inside[s:e] = within.any(axis=1)
        inside_adh[s:e] = (within & adhesive_b[None, :]).any(axis=1)
    return inside, inside_adh


def run_diffusion(cfg, rng=None, progress=False, progress_cb=None,
                  frame_cb=None, graphics_every=None):
    if rng is None:
        rng = np.random.default_rng()
    B, C, comp, cleft, cube, adh = (cfg.basic, cfg.balls, cfg.computation,
                                    cfg.cleft, cfg.cube, cfg.adhesion)
    npar = B.NumberOfParticles
    dt = comp.GTimeStep
    n_steps = int(round(B.TotalTime / dt))
    gevery = int(graphics_every) if graphics_every else max(1, n_steps // 50)
    step_rad = comp.StepRad
    dist_bins, time_bins = comp.DistanceBins, comp.TimeBins
    synapse = B.ReleaseSource.lower() == "synapse"
    half_w = cleft.CleftWidth / 2.0

    p_bound_removal = (1 - np.exp(-dt / B.TimeDecay)) if B.UseExponentialDecay else 0.0
    p_adh_loss = (1 - np.exp(-dt / C.AdhesionDecayTime)) if C.UseAdhesionDecay else 0.0

    dead_steps = np.clip(np.round(4000 + 2000 * rng.standard_normal(npar)), 1000, 7000).astype(int)
    prob_uptaken = 1 - B.UnboundProb
    deadhesion_thresh = np.where(rng.random(npar) < prob_uptaken, 0, dead_steps)

    release_time, release_counts = _assign_releases(npar, B.ReleaseTimes, B.ReleaseFractions, rng)
    released = release_time <= 0

    astro_center = astro_radius = None
    x = np.zeros(npar); y = np.zeros(npar); z = np.zeros(npar)
    if synapse:
        first_ball_index = 0
    elif B.ReleaseSource.lower() == "astrocyte":
        first_ball_index = 1
        astro_radius = B.AstrocyteRadius
        astro_center = np.array([0.0, B.AstrocyteCenterDistance, 0.0])
        theta = 2 * np.pi * rng.random(npar)
        phi_max = np.arcsin(B.ReleasePatchRadius / astro_radius)
        phi = phi_max * (2 * rng.random(npar) - 1)
        xp = astro_radius * np.sin(phi) * np.cos(theta)
        yp = astro_radius * np.sin(phi) * np.sin(theta)
        zp = astro_radius * np.cos(phi)
        d0 = -astro_center / np.linalg.norm(astro_center)
        init = np.array([0.0, 0.0, 1.0])
        axis = np.cross(init, d0); s = np.linalg.norm(axis); c = init @ d0
        if s < 1e-12:
            Rm = np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
        else:
            axis = axis / s
            K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
            Rm = np.eye(3) + s * K + (1 - c) * (K @ K)
        pts = (Rm @ np.vstack([xp, yp, zp])).T + astro_center
        x = pts[:, 0] + 10.0 * (rng.random(npar) - 0.5)
        y = pts[:, 1] + 10.0 * (rng.random(npar) - 0.5)
        z = pts[:, 2] + 10.0 * (rng.random(npar) - 0.5)
    else:
        raise ValueError("ReleaseSource must be 'synapse' or 'astrocyte'")

    xd, yd, zd, radius, adhesive = _place_balls(cfg, rng, first_ball_index, astro_center, astro_radius)
    if B.ReleaseSource.lower() == "astrocyte":
        xd[0], yd[0], zd[0], radius[0], adhesive[0] = (0.0, B.AstrocyteCenterDistance, 0.0, astro_radius, True)
    initial_adhesive = adhesive.copy()
    bx = xd[first_ball_index:]; by = yd[first_ball_index:]; bz = zd[first_ball_index:]
    br2 = radius[first_ball_index:] ** 2

    betta, occ_se = _occupancy(cfg, rng, xd, yd, zd, radius, first_ball_index, astro_center, astro_radius)

    stak1 = np.ones(npar, int)
    active = released.copy()
    time_adhesive = np.zeros(npar, int)

    time_ms = np.arange(n_steps + 1) * dt
    num_spill = np.zeros(n_steps + 1); frac_spill = np.zeros(n_steps + 1)
    num_free_spill = np.zeros(n_steps + 1); frac_free_spill = np.zeros(n_steps + 1)
    fp_time = np.full(npar, np.nan); fp_lat = np.full(npar, np.nan)
    difx = np.zeros(n_steps + 1)

    radial_edges = np.linspace(cleft.RadiusofCleft, cube.AxisSize, dist_bins + 1)
    avg_free = np.zeros((dist_bins, time_bins + 1)); avg_bound = np.zeros((dist_bins, time_bins + 1))
    avg_free[:, 0] = radial_edges[1:] / 1000.0
    avg_bound[:, 0] = radial_edges[1:] / 1000.0
    sample_every = int(round(B.TotalTime / time_bins / dt))
    sample_col = 0

    snapshots = {}
    snap_steps = {int(round(t / dt)): ("%g" % t) for t in comp.SnapshotTimes}
    loc = comp.StepLocalMovement * step_rad * 2

    for i in range(1, n_steps + 1):
        t_start = (i - 1) * dt
        t_now = i * dt

        newly = (~released) & (release_time <= t_start + 1e-15)
        released[newly] = True
        active[newly] = True

        rem = active & (stak1 > 1) & (rng.random(npar) < p_bound_removal)
        active[rem] = False

        if p_adh_loss > 0 and adhesive.any():
            lost = adhesive & (rng.random(adhesive.size) < p_adh_loss)
            adhesive[lost] = False
        adhesive_b = adhesive[first_ball_index:]

        ox = x.copy(); oy = y.copy(); oz = z.copy()
        if synapse:
            in_cleft = active & (np.hypot(x, y) <= cleft.RadiusofCleft) & (np.abs(z) <= half_w)
        else:
            in_cleft = np.zeros(npar, bool)
        free = active & (stak1 == 1)
        bound = active & (stak1 > 1)
        free_ext = free & ~in_cleft
        bound_ext = bound & ~in_cleft
        ext = free_ext | bound_ext

        nx = ox.copy(); ny = oy.copy(); nz = oz.copy()
        nx[free_ext] = ox[free_ext] + step_rad * 2 * (rng.random(free_ext.sum()) - 0.5)
        ny[free_ext] = oy[free_ext] + step_rad * 2 * (rng.random(free_ext.sum()) - 0.5)
        nz[free_ext] = oz[free_ext] + step_rad * 2 * (rng.random(free_ext.sum()) - 0.5)

        if in_cleft.any():
            m = in_cleft.sum()
            nx[in_cleft] = ox[in_cleft] + step_rad * 2 * (rng.random(m) - 0.5)
            ny[in_cleft] = oy[in_cleft] + step_rad * 2 * (rng.random(m) - 0.5)
            zc = oz[in_cleft] + 0.1 * step_rad * 2 * (rng.random(m) - 0.5)
            zc = np.where(zc > half_w, 2 * half_w - zc, zc)
            zc = np.where(zc < -half_w, -2 * half_w - zc, zc)
            nz[in_cleft] = zc

        inside_any = np.zeros(npar, bool); inside_adh = np.zeros(npar, bool)
        if ext.any():
            ia, iad = _containment(nx[ext], ny[ext], nz[ext], bx, by, bz, br2, adhesive_b)
            inside_any[ext] = ia; inside_adh[ext] = iad

        cand = free_ext & inside_adh
        if cand.any():
            time_adhesive[cand] += 1
            p = adh.MaxProbAdhesive * (1 - np.exp(-(dt * time_adhesive) / adh.TimeINsideAdhesiveZone))
            do_bind = cand & (rng.random(npar) < p)
            stak1[do_bind] = 2
            time_adhesive[do_bind] = 0

        lm = ext & inside_any
        if lm.any():
            k = lm.sum()
            nx[lm] = ox[lm] + loc * (rng.random(k) - 0.5)
            ny[lm] = oy[lm] + loc * (rng.random(k) - 0.5)
            nz[lm] = oz[lm] + loc * (rng.random(k) - 0.5)

        if bound_ext.any():
            stak1[bound_ext] += 1
            reset = bound_ext & (stak1 == deadhesion_thresh)
            stak1[reset] = 1

        if synapse:
            cd = np.sqrt(nx * nx + ny * ny + nz * nz)
            reject = active & (cd <= cleft.SynapticEnvelopeRadius) & (np.abs(nz) > half_w)
            nx[reject] = ox[reject]; ny[reject] = oy[reject]; nz[reject] = oz[reject]

        x, y, z = nx, ny, nz

        age = t_now - release_time
        inbox = (np.abs(x) < cube.Size / 2) | (np.abs(y) < cube.Size / 2) | (np.abs(z) < cube.Size / 2)
        valid = active & (age > 0) & inbox
        cnt = int(valid.sum())
        difx[i] = (np.sum(x[valid] ** 2 / (2 * age[valid])) / cnt) / 1e6 if cnt else 0.0

        rad = np.sqrt(x * x + y * y + z * z)
        spill = active & (rad >= cleft.SpilloverRadius)
        free_spill = spill & (stak1 == 1)
        num_spill[i] = spill.sum(); frac_spill[i] = num_spill[i] / npar
        num_free_spill[i] = free_spill.sum(); frac_free_spill[i] = num_free_spill[i] / npar
        new_cross = np.isnan(fp_time) & spill
        fp_time[new_cross] = t_now
        fp_lat[new_cross] = t_now - release_time[new_cross]

        if i % sample_every == 0 and sample_col < time_bins:
            sample_col += 1
            rr = np.where(rad < cleft.RadiusofCleft, 0.0, rad)
            hb, _ = np.histogram(np.where(active & (stak1 > 1), rr, 0.0), bins=radial_edges)
            hf, _ = np.histogram(np.where(active & (stak1 <= 1), rr, 0.0), bins=radial_edges)
            avg_bound[:, sample_col] = hb
            avg_free[:, sample_col] = hf

        if i in snap_steps:
            state = np.where((stak1 > 1) | (~active), 0, stak1).astype(float)
            snapshots[snap_steps[i]] = np.vstack([x, y, z, state])

        if frame_cb is not None and (i == 1 or i % gevery == 0 or i == n_steps):
            frame = dict(t=t_now, x=x.copy(), y=y.copy(), z=z.copy(),
                         active=active.copy(), bound=(stak1 > 1).copy())
            if i == 1:
                frame["balls"] = dict(xyz=np.column_stack([xd, yd, zd]),
                                      radius=radius.copy(), adhesive=initial_adhesive.copy())
            frame_cb(frame)

        if i % max(1, n_steps // 20) == 0:
            if progress:
                print("  diffusion %.0f%%" % (100 * i / n_steps), flush=True)
            if progress_cb is not None:
                progress_cb(100 * i / n_steps)

    begin = int(B.BeginTime * B.TotalTime / dt)
    end = int(B.EndTime * B.TotalTime / dt)
    window = difx[begin:end]
    mean_d = float(np.mean(window)) if window.size else float("nan")
    std_d = float(np.std(window, ddof=1)) if window.size > 1 else float("nan")

    radial_inner_um = radial_edges[:-1] / 1000.0
    radial_outer_um = radial_edges[1:] / 1000.0
    metadata = dict(
        FormatVersion=1, TotalTime_ms=B.TotalTime, TimeBins=time_bins, DistanceBins=dist_bins,
        TimePoints_ms=np.arange(1, time_bins + 1) * (B.TotalTime / time_bins),
        RadialEdges_nm=radial_edges, InnerRadius_um=radial_inner_um, OuterRadius_um=radial_outer_um,
        ShellVolume_um3=(4 / 3) * np.pi * (radial_outer_um ** 3 - radial_inner_um ** 3),
        NumberOfParticles=npar, ReleaseSource=B.ReleaseSource,
        ReleaseTimes_ms=np.asarray(B.ReleaseTimes), ReleaseFractions=np.asarray(B.ReleaseFractions),
        DiffusionCoefficient_nm2_per_ms=comp.DifCoefINIT, BrownianTimeStep_ms=dt,
        CleftRadius_nm=cleft.RadiusofCleft, SynapticEnvelopeRadius_nm=cleft.SynapticEnvelopeRadius,
        SpilloverRadius_nm=cleft.SpilloverRadius, RadialAnalysisLimit_nm=cube.AxisSize,
        OccupiedFraction=betta,
    )

    return DiffusionResult(
        average_free=avg_free, average_bound=avg_bound, time_ms=time_ms,
        number_spillover=num_spill, fraction_spillover=frac_spill,
        number_free_spillover=num_free_spill, fraction_free_spillover=frac_free_spill,
        first_passage_time=fp_time, first_passage_latency=fp_lat,
        release_time=release_time, occupied_fraction=betta, free_fraction=1 - betta,
        occupancy_se=occ_se, mean_diffusion=mean_d, std_diffusion=std_d,
        balls_xyz=np.column_stack([xd, yd, zd]), balls_radius=radius,
        balls_adhesive=initial_adhesive, metadata=metadata, snapshots=snapshots)
