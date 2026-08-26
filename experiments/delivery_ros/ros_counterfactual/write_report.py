import json
from pathlib import Path

metrics_path = Path('/home/ubuntu/gutibm-campaign/ros-counterfactual/metrics.json')
report_path = Path('/home/ubuntu/gutibm-campaign/ros-counterfactual/REPORT.md')
data = json.loads(metrics_path.read_text())
arms = data['arms']

def f(v):
    if isinstance(v, float): return f'{v:.8g}'
    return str(v)

def row(vals):
    return '| ' + ' | '.join(f(v) for v in vals) + ' |'

def traj_at(traj, times, key):
    by = {round(float(x['time_s'])): x for x in traj}
    return [by.get(t, {}).get(key, '—') for t in times]

lines = ['# ROS counterfactual campaign report', '',
         ('Valid arms are the two Group A arms and four corrected Group B arms. '
          'Full per-step trajectories, SOS component ledgers, and all nonzero '
          'clip intervals are in `metrics.json`; `metrics.csv` contains the '
          'per-arm summary rows.'), '',
         '## Build and policy', '',
         '- Commit: `abf393750453673ebb1f4c4ab54f447797a647fd`',
         '- Binary SHA-256: `859fc77a8e5280b2eb0271ae0603cbf4965adf891ff5edc2ba4fd2e3817221d6`',
         '- Stop policy: relative clip threshold `1e-6` against cumulative agent growth uptake; funded-minus-realized checked at every ledger record.',
         '- `delivery_reduction_*`: absent from every emitted HDF5 output.', '',
         '## Run summary', '',
         '| Arm | dx (µm) | k_ROS | founders | final N | divisions | termination | time (s) | wall (s) |',
         '|---|---:|---:|---:|---:|---:|---|---:|---:|']
for x in arms:
    lines.append(row([x['arm'],x['grid_dx_um'],x['oxygen_k_ROS'],x['configured_founders'],x['final_N'],x['cumulative_divisions'],x['termination_cause'],x['termination_time_s'],x['wrapper_wall_seconds']]))
lines += ['', '## Population loss closure', '', '| Arm | lysis | colicin | CDI | washout-trapped | boundary outflow | channel sum | observed N decrement | divisions − losses | closure |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
for x in arms:
    p=x['population_loss']
    lines.append(row([x['arm'],p['mortality_lysis'],p['mortality_colicin'],p['mortality_cdi'],p['outflow_washout_trapped'],p['outflow_boundary'],p['gross_loss_channel_sum'],p['net_N_decrease'],p['adjusted_N_decrease_for_divisions_and_immigrations'],p['gross_channel_sum_matches_adjusted_decrease']]))
lines += ['', 'Loss channels include mortality and transport. The channel sum closes against `initial N + divisions + immigrations − final N`; therefore boundary/washout entries are transported out, while lysis/colicin/CDI entries are deaths.', '', '## Oxygen and carbon ledgers', '', '| Arm | O₂ growth funded/demanded | O₂ growth fraction | O₂ maintenance funded/demanded | O₂ maint shortfall | O₂ agent/flora removal | agent share | C growth funded/demanded | C maintenance funded/demanded | C maint shortfall | C clips (relative) |', '|---|---|---:|---|---:|---|---:|---|---|---:|---|']
for x in arms:
    o,c=x['oxygen'],x['carbon']
    lines.append(row([x['arm'],f'{f(o["funded_growth"])} / {f(o["demanded_growth"])}',f(o['funded_growth_fraction']),f'{f(o["funded_maintenance"])} / {f(o["demanded_maintenance"])}',f(o['maintenance_shortfall']),f'{f(o["agent_total_removal"])} / {f(o["flora_removal"])}',f(o['agent_share']),f'{f(c["funded_growth"])} / {f(c["demanded_growth"])}',f'{f(c["funded_maintenance"])} / {f(c["demanded_maintenance"])}',f(c['maintenance_shortfall']),f'{f(c["reaction_clips"])} ({f(c["reaction_clip_relative_to_growth_uptake"])} final; max {f(x["max_clip_relative_to_cumulative_growth_uptake"]["carbon"])})']))
lines += ['', 'All funded-minus-realized final residuals are zero at reported precision. The largest positive cumulative floating residual seen during the stepwise check was `4.04e-28 mol` (round-off).', '', '## Clips', '', '| Arm | O₂ cumulative clip | O₂ relative | C cumulative clip | C relative | nonzero interval records |', '|---|---:|---:|---:|---:|---:|']
for x in arms:
    lines.append(row([x['arm'],x['oxygen']['reaction_clips'],f(x['max_clip_relative_to_cumulative_growth_uptake']['oxygen']),x['carbon']['reaction_clips'],f(x['max_clip_relative_to_cumulative_growth_uptake']['carbon']),len(x['nonzero_clip_intervals'])]))
lines += ['', 'For `B_ros0_res2`, the first nonzero carbon interval is recorded in `metrics.json` at step 120 / 7200 s; maximum interval-relative carbon clip is below `1e-6`, so the run continued. Other arms have no nonzero clips.', '', '## SOS hazard components', '', '| Arm | cumulative basal | post-division | nuclease cross-induction | ROS | total | SOS inductions | final interval components (basal, post, nuclease, ROS) | final per-live-agent interval components |', '|---|---:|---:|---:|---:|---:|---:|---|---|']
for x in arms:
    s=x['sos']; z=s['interval_and_cumulative'][-1]; cc=z['cumulative_components']; ii=z['interval_components']; pp=z['interval_components_per_live_agent']
    lines.append(row([x['arm'],cc['sos_basal_rate'],cc['sos_post_division_rate'],cc['sos_nuclease_cross_induction_rate'],cc['sos_ros_rate'],z['cumulative_total'],s['cumulative_inductions'],tuple(ii.values()),tuple(pp.values())]))
lines += ['', 'Each SOS component sum equals the emitted total (`cumulative_component_sum_matches_total = true` for all arms). Full interval and cumulative records include per-live-agent means in JSON.', '', '## Fermentation, acetate, grids', '', '| Arm | mean fermentation | final fermentation | final acetate mean | acid inhibition | final O₂ mean/min | final C mean/min |', '|---|---:|---:|---:|---|---|---|']
for x in arms:
    g=x['final_grid']; lines.append(row([x['arm'],x['fermentation_fraction_mean'],x['fermentation_fraction_final'],x['final_acetate_mean'],x['acid_inhibition_enabled'],f'{g["oxygen"]["domain_mean"]} / {g["oxygen"]["domain_minimum"]}',f'{g["carbon"]["domain_mean"]} / {g["carbon"]["domain_minimum"]}']))
lines += ['', 'Final-grid agent-cell concentrations and cell/domain-mean ratios are retained per species in `metrics.json` under `final_grid`.', '', '## Sampled trajectories', '', 'Samples use exact saved agent/summary times where available. `—` means the run ended before that sample.', '']
times=[0,600,1200,3600,7200,10800,14400,18000,21600]
for x in arms:
    lines += [f'### {x["arm"]}', '', '| time (s) | N | divisions interval | mean mu_realized (s⁻¹) | mean biomass (mol) | bacteriostatic | washout-trapped |', '|---:|---:|---:|---:|---:|---:|---:|']
    nby={round(float(a['time_s'])):a for a in x['N_trajectory']}
    dby={round(float(a['time_s'])):a for a in x['raw_rows']}
    mby={round(float(a['time_s'])):a for a in x['mean_mu_realized_trajectory']}
    bby={round(float(a['time_s'])):a for a in x['mean_biomass_trajectory']}
    sby={round(float(a['time_s'])):a for a in x['live_agent_stocks']}
    for t in times:
        n=nby.get(t,{}); d=dby.get(t,{}); m=mby.get(t,{}); b=bby.get(t,{}); s=sby.get(t,{})
        if n or d or m or b or s:
            lines.append(row([t,n.get('N','—'),d.get('divisions','—'),m.get('mean_mu_realized','—'),b.get('mean_biomass','—'),s.get('bacteriostatic_live_agents','—'),s.get('washout_trapped_live_agents','—')]))
    lines.append('')
lines += ['## Invalid superseded arm', '', '- `B_ros0_res2` 100-founder output is excluded from `arms` metrics.', '- Preserved at `invalid_B_ros0_res2_100_founders/` and recorded in `metrics.json` under `invalid_arms`.', '- Reason: first dispatch used 80 strain-1 founders plus 20 strain-2 founders instead of 60 + 20 = 80 total.', '']
with Path('/home/ubuntu/gutibm-campaign/ros-counterfactual/REPORT.md').open(
    'w'
) as report_file:
    report_file.write('\n'.join(lines))
print(report_path)
