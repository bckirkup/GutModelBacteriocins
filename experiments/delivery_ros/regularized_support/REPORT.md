# Population-scale acceptance gate for the regularized delivery sink (PR #335)

Binary: `/home/ubuntu/campaign-src-335/build-serial/gut_ibm`
(`cbb68fed032700d5e94fd44e6a6d6b7bb4720ad2fc3047cc4aa179f29f5f589e`), built from
commit `2a16cc7` of `devin/1787623566-regularized-delivery-sink` in a detached
worktree so the CI branch could keep moving.

Arms are the #334 far-field campaign configs. Only `grid_dx` and
`metabolism.delivery_far_field_radius` vary: 80 founders, anatomic placement,
Robin oxygen boundary, `oxygen.k_ROS = 0` (the retired ambient-ROS mortality),
6 h horizon, seed 1001. Nothing is fitted.

## Result

| arm | grid | radius | final N | peak N | divisions | funded O2 | funded C | O2 withheld (mol) | C withheld (mol) | O2 retries | C retries |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| P_ff10_res2 | 2 um | 10 um | 201 | 263 | 313 | 0.988 | 0.9967 | 1.29e-14 | 0 | 63 | 0 |
| P_ff10_res4 | 4 um | 10 um | 207 | 256 | 317 | 0.9885 | 0.9966 | 1.05e-14 | 0 | 59 | 0 |
| P_ff10_res6 | 6 um | 10 um | 204 | 258 | 322 | 0.9883 | 0.9966 | 1.15e-14 | 0 | 61 | 0 |
| P_ff0_res2 | 2 um | 0 | 76 | 156 | 103 | 0.620 | 0.6865 | 8.39e-13 | 3.43e-14 | 392 | 419 |
| P_ff0_res4 | 4 um | 0 | 123 | 232 | 212 | 0.9009 | 0.9941 | 3.44e-13 | 4.48e-16 | 257 | 44 |
| P_ff0_res6 | 6 um | 0 | 171 | 254 | 263 | 0.9479 | 0.9972 | 2.04e-13 | 0 | 169 | 0 |

All arms reached `horizon_reached`; oxygen reaction clips are exactly zero
everywhere and carbon clips are zero except the 1.15e-20 already present in the
`P_ff0_res2` control.

- **Acceptance criterion met.** With the radius on, final N is 201/207/204
  across a 3x grid change (1.03x spread, i.e. inside run-to-run stochastic
  variation) and funded oxygen fraction is 0.988 at every resolution
  (1.001x spread). With it off, the same quantities span 2.25x and 1.53x.
- **Both halves are asserted.** The radius-off arms reproduce the pre-#335
  campaign exactly (76 / 123 / 171 final N, funded fraction 0.620 / 0.901 /
  0.948), which confirms the radius-zero path is untouched and that the contrast
  is not an inert code path — the failure mode that made #334 look solved.
- **Mechanism confirmed.** The retry ledger from #333 is the discriminator:
  withheld oxygen falls 65x at 2 um (8.39e-13 -> 1.29e-14) and carbon retries
  fall 419 -> 0. What remains with the radius on is grid-independent
  (1.05e-14 - 1.29e-14 mol, 59-63 events), i.e. a genuine local depletion limit
  rather than a `dx^3` bookkeeping limit.
- **The artifact was suppressing growth, not causing it.** Carrying capacity
  with the radius on (~204) exceeds every radius-off arm including the coarsest,
  so the previously reported 76-vs-171 range was a floor set by discretisation.

## What this does not establish

- The 6 h horizon is short relative to mucus turnover; this is a supply-limited
  carrying-capacity measurement, not a persistence measurement.
- ROS hazard is off in all six arms. Nothing here bears on the lysis
  coefficient, which is still an in-vitro-derived prior pending comparison to
  in-vivo displacement kinetics.
- The default radius remains 0.0, so nothing shipped changes until the default
  is flipped. The case for flipping it to 10 um now rests on this campaign plus
  the physical argument that single-cell depletion at 10 um is ~0.008% of C_inf,
  so the averaging radius is not doing biological work.
- Any carbon-limited carrying capacity previously reported (#314) was measured
  under the single-voxel deposition and moves with this change.

## Correction

The default `delivery_far_field_radius` is `10 µm` since #336, not `0.0`.
