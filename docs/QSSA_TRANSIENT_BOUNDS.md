# QSSA Transient vs Dynamic PDE Bounds (Issue #472 / GAP-11)

This note records the measured receptor-occupancy and cumulative kill-probability errors introduced by GutIBM's Quasi-Steady-State Approximation (QSSA) bacteriocin solver relative to the exact transient 3D advection-diffusion-decay PDE over a single biological timestep (`bio_dt = 60 s`) and a full lysis burst window (`burst_release_tau = 300 s`, five `60 s` steps).

All numbers below are produced and continuously verified by the `qssa_transient_bounds` CTest target (`tests/test_qssa_transient_bounds.cpp`).

## 1. Mathematical Reference & Self-Consistency

For a point lysis burst at $\mathbf{x}_0 = (x_0, y_0, z_0)$ with exponential release rate
$$Q(t) = Q_0 \exp(-t / \tau_{\text{rel}}), \qquad Q_0 = \frac{N_{\text{burst}} / N_A}{\tau_{\text{rel}}}, \qquad \tau_{\text{rel}} = 300\text{ s},$$
in a sealed Neumann mucus channel $z \in [0, H]$ with uniform wall-parallel advection $\mathbf{U} = (u_x, 0, 0)$, effective diffusivity $D_{\text{eff}} = D_{\text{free}} / R$, and first-order protease decay $k_{\text{dec}} = \ln(2) / t_{1/2}$, the exact transient PDE concentration is the Duhamel convolution of the bounded Method-of-Images 3D heat kernel:
$$C_{\text{PDE}}(\mathbf{x}, t) = \int_0^t Q(t - s)\, e^{-k_{\text{dec}} s} \sum_{n=-\infty}^{\infty} \frac{\exp\!\left(-\frac{(x - x_0 - u_x s)^2 + (y - y_0)^2}{4 D_{\text{eff}} s}\right)}{(4\pi D_{\text{eff}} s)^{3/2}} \left[e^{-\frac{(z - z_0 - 2nH)^2}{4 D_{\text{eff}} s}} + e^{-\frac{(z + z_0 - 2nH)^2}{4 D_{\text{eff}} s}}\right] ds.$$

- **Steady-state verification**: Evaluating $C_{\text{PDE}}$ for a constant source ($\tau_{\text{rel}} \to \infty$) at $t = 2 \times 10^4\text{ s}$ reproduces `GreensFunction::concentration_bounded` to relative error **`8.42e-14`**, confirming exact operator agreement. Any discrepancy at finite $t$ is strictly the QSSA transient-arrival and discrete-step approximation error.
- **Receptor occupancy & kill probability**:
  - $\theta(C) = C / (C + K_d^{\text{app}})$, evaluated both at uncompeted BtuB affinity ($K_d^{\text{app}} = 5.0 \times 10^{-7}\text{ mol/m}^3$, high-sensitivity upper bound) and at the shipped corrinoid-competed affinity ($K_d^{\text{app}} = 5.0 \times 10^{-4}\text{ mol/m}^3$ at $[B_{12}] = 10^{-3}\text{ mol/m}^3$).
  - Single-hit kill rate $k_{\text{kill}} = 1.0 \times 10^{-3}\text{ s}^{-1}$ (`ReceptorConfig::kill_rate_colicin`).
  - Dynamic PDE kill probability integrates $k_{\text{kill}}\,\theta(C_{\text{PDE}}(t))$ at 5 s trapezoidal substeps; GutIBM QSSA evaluates $C_{\text{QSSA}}(t_k) = C_{\text{GF}} e^{-t_k / \tau_{\text{rel}}}$ at $t_k = k \cdot 60\text{ s}$ ($k = 0, \dots, 4$).

## 2. Test Geometries & Plasmids

- **Geometries**:
  - `60um_compact`: $60 \times 60 \times 50\text{ }\mu\text{m}$ ($z_0 = 25\text{ }\mu\text{m}$), $u_x = 0.833\text{ }\mu\text{m/s}$.
  - `200um_mucus`: $200 \times 200 \times 100\text{ }\mu\text{m}$ ($z_0 = 50\text{ }\mu\text{m}$), $u_x = 0.833\text{ }\mu\text{m/s}$.
- **Plasmids** (from `PlasmidLibrary`):
  - **ColE1** (retarded pore-former): $N_{\text{burst}} = 10^5$, $R = 50.2$, $D_{\text{eff}} = 7.96 \times 10^{-13}\text{ m}^2/\text{s}$, $t_{1/2} = 1800\text{ s}$.
  - **ColE2** (weakly retarded nuclease complex): $N_{\text{burst}} = 5 \times 10^4$, $R = 2.04$, $D_{\text{eff}} = 1.96 \times 10^{-11}\text{ m}^2/\text{s}$, $t_{1/2} = 1800\text{ s}$.

## 3. Measured Bounds

### 3.1 Compact `60 um` Geometry ($60 \times 60 \times 50\text{ }\mu\text{m}$)

| Plasmid | Direction | Distance (um) | Max abs(d_theta) (60 s / 300 s, uncomp.) | P_kill,QSSA (60 s) | P_kill,PDE (60 s) | abs(dP_kill) (60 s) | P_kill,QSSA (300 s) | P_kill,PDE (300 s) | abs(dP_kill) (300 s) | abs(dP_kill) (300 s, shipped B12) |
|---------|-----------|-------------------:|---------------------------------------------:|------------------------------:|-----------------------------:|----------------------------------:|-------------------------------:|------------------------------:|-----------------------------------:|------------------------------------------------------:|
| ColE1 | downstream | 5.0 | 0.0412 | 0.0558 | 0.0531 | 0.0027 | 0.2445 | 0.2413 | 0.0032 | 0.0005 |
| ColE1 | downstream | 10.0 | 0.6783 | 0.0535 | 0.0472 | 0.0063 | 0.2314 | 0.2250 | 0.0064 | 0.0003 |
| ColE1 | downstream | 20.0 | 0.8457 | 0.0495 | 0.0350 | 0.0145 | 0.2092 | 0.1964 | 0.0128 | 0.0002 |
| ColE1 | downstream | 25.0 | 0.8139 | 0.0477 | 0.0292 | 0.0185 | 0.1996 | 0.1840 | 0.0156 | 0.0001 |
| ColE1 | transverse | 10.0 | 0.0539 | 0.0033 | 0.0025 | 0.0008 | 0.0117 | 0.0104 | 0.0013 | — |
| ColE1 | upstream | 10.0 | 0.0003 | 0.0000 | 0.0000 | 0.0000 | 0.0001 | 0.0001 | 0.0000 | — |
| ColE1 | transverse | 20.0 | 0.0002 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | — |
| ColE1 | upstream | 20.0 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | — |
| ColE2 | downstream | 5.0 | 0.0713 | 0.0211 | 0.0176 | 0.0035 | 0.0787 | 0.0721 | 0.0066 | < 0.0001 |
| ColE2 | downstream | 10.0 | 0.1062 | 0.0137 | 0.0100 | 0.0037 | 0.0498 | 0.0443 | 0.0055 | < 0.0001 |
| ColE2 | downstream | 20.0 | 0.1249 | 0.0089 | 0.0049 | 0.0040 | 0.0320 | 0.0272 | 0.0048 | < 0.0001 |
| ColE2 | downstream | 25.0 | 0.1223 | 0.0079 | 0.0038 | 0.0041 | 0.0283 | 0.0234 | 0.0049 | < 0.0001 |
| ColE2 | transverse | 10.0 | 0.0901 | 0.0113 | 0.0082 | 0.0031 | 0.0408 | 0.0362 | 0.0046 | — |
| ColE2 | upstream | 10.0 | 0.0753 | 0.0093 | 0.0066 | 0.0027 | 0.0332 | 0.0293 | 0.0039 | — |
| ColE2 | transverse | 20.0 | 0.0824 | 0.0059 | 0.0031 | 0.0028 | 0.0208 | 0.0175 | 0.0033 | — |
| ColE2 | upstream | 20.0 | 0.0531 | 0.0037 | 0.0020 | 0.0017 | 0.0132 | 0.0110 | 0.0022 | — |

### 3.2 Colonic Mucus `200 um` Geometry ($200 \times 200 \times 100\text{ }\mu\text{m}$)

| Plasmid | Direction | Distance (um) | Max abs(d_theta) (60 s / 300 s, uncomp.) | P_kill,QSSA (60 s) | P_kill,PDE (60 s) | abs(dP_kill) (60 s) | P_kill,QSSA (300 s) | P_kill,PDE (300 s) | abs(dP_kill) (300 s) | abs(dP_kill) (300 s, shipped B12) |
|---------|-----------|-------------------:|---------------------------------------------:|------------------------------:|-----------------------------:|----------------------------------:|-------------------------------:|------------------------------:|-----------------------------------:|------------------------------------------------------:|
| ColE1 | downstream | 5.0 | 0.0412 | 0.0558 | 0.0531 | 0.0027 | 0.2445 | 0.2413 | 0.0032 | 0.0005 |
| ColE1 | downstream | 10.0 | 0.6783 | 0.0535 | 0.0472 | 0.0063 | 0.2314 | 0.2250 | 0.0064 | 0.0003 |
| ColE1 | downstream | 20.0 | 0.8457 | 0.0495 | 0.0350 | 0.0145 | 0.2092 | 0.1964 | 0.0128 | 0.0002 |
| ColE1 | downstream | 30.0 | 0.7843 | 0.0460 | 0.0237 | 0.0223 | 0.1908 | 0.1728 | 0.0180 | 0.0001 |
| ColE1 | downstream | 50.0 | 0.6837 | 0.0402 | 0.0060 | 0.0342 | 0.1622 | 0.1362 | 0.0260 | 0.0001 |
| ColE1 | transverse | 10.0 | 0.0539 | 0.0033 | 0.0025 | 0.0008 | 0.0117 | 0.0104 | 0.0013 | — |
| ColE1 | upstream | 10.0 | 0.0003 | 0.0000 | 0.0000 | 0.0000 | 0.0001 | 0.0001 | 0.0000 | — |
| ColE1 | transverse | 20.0 | 0.0002 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | — |
| ColE1 | upstream | 20.0 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | — |
| ColE2 | downstream | 5.0 | 0.0567 | 0.0202 | 0.0173 | 0.0029 | 0.0752 | 0.0692 | 0.0060 | < 0.0001 |
| ColE2 | downstream | 10.0 | 0.0833 | 0.0123 | 0.0096 | 0.0027 | 0.0447 | 0.0401 | 0.0046 | < 0.0001 |
| ColE2 | downstream | 20.0 | 0.0932 | 0.0071 | 0.0044 | 0.0027 | 0.0251 | 0.0217 | 0.0034 | < 0.0001 |
| ColE2 | downstream | 30.0 | 0.0808 | 0.0051 | 0.0024 | 0.0027 | 0.0179 | 0.0148 | 0.0031 | < 0.0001 |
| ColE2 | downstream | 50.0 | 0.0569 | 0.0034 | 0.0009 | 0.0025 | 0.0120 | 0.0091 | 0.0029 | < 0.0001 |
| ColE2 | transverse | 10.0 | 0.0703 | 0.0101 | 0.0078 | 0.0023 | 0.0365 | 0.0326 | 0.0039 | — |
| ColE2 | upstream | 10.0 | 0.0585 | 0.0083 | 0.0063 | 0.0020 | 0.0296 | 0.0263 | 0.0033 | — |
| ColE2 | transverse | 20.0 | 0.0607 | 0.0046 | 0.0028 | 0.0018 | 0.0161 | 0.0138 | 0.0023 | — |
| ColE2 | upstream | 20.0 | 0.0388 | 0.0029 | 0.0018 | 0.0011 | 0.0102 | 0.0087 | 0.0015 | — |

## 4. Key Takeaways & Envelope

1. **Direction of bias (causality)**: Instantaneous QSSA always overestimates first-step (`60 s`) and burst-integrated (`300 s`) kill probability relative to the causal transient PDE because QSSA sets the full downstream plume at $t = 0$ instead of waiting for advective-diffusive arrival ($\tau_{\text{adv}} = r / u_x$).
2. **Instantaneous occupancy vs integrated kill probability**:
   - For strongly retarded **ColE1** ($D_{\text{eff}} = 7.96 \times 10^{-13}\text{ m}^2/\text{s}$) under uncompeted $K_d$, peak instantaneous receptor-occupancy error $|\Delta\theta|$ reaches **0.846** at $r = 20\text{ }\mu\text{m}$ downstream during the early arrival lag ($t < 24\text{ s}$). Because receptors rapidly saturate once the plume arrives ($\tau_{\text{adv}} \approx 24\text{ s}$ at $20\text{ }\mu\text{m}$, $60\text{ s}$ at $50\text{ }\mu\text{m}$), the **cumulative 300 s kill-probability error** stays bounded at **`<= 0.0156`** in the `60 um` box ($r \le 25\text{ }\mu\text{m}$) and **`<= 0.0260`** in the `200 um` box ($r \le 50\text{ }\mu\text{m}$), though the first-step (`60 s`) error at $r = 50\text{ }\mu\text{m}$ is a large relative overestimate (`0.0402` vs `0.0060`, `6.7x`).
   - For weakly retarded **ColE2** ($D_{\text{eff}} = 1.96 \times 10^{-11}\text{ m}^2/\text{s}$), diffusive equilibration is fast ($\tau_{\text{diff}} \sim 1$–$30\text{ s}$), keeping $|\Delta\theta| \le 0.125$ (`60 um`) / $\le 0.093$ (`200 um`) and cumulative 300 s $|\Delta P_{\text{kill}}| \le 0.0066$ across all probed directions and distances.
   - Under the shipped corrinoid-competed BtuB affinity ($[B_{12}] = 10^{-3}\text{ mol/m}^3$), cumulative 300 s $|\Delta P_{\text{kill}}|$ is **`<= 0.0005`** across both geometries.
3. **Geometry effect**: Moving from `H = 50 um` (`60um_compact`) to `H = 100 um` (`200um_mucus`) leaves ColE1 downstream plume concentrations identical (the retarded plume stays well inside the channel walls over 300 s) while reducing ColE2 concentrations and transient errors by ~20–25% due to reduced vertical image confinement.
