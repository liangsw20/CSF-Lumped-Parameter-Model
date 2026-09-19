# CSF_Model — multi-compartment lumped-parameter model of intracranial fluid dynamics

Supplementary code for the manuscript. It contains the Simulink model, the
parameter-computation scripts, the batch runners for every scenario reported in
Tables 2 and 3, and the scripts that draw the figures.

---

## 1. Requirements

* **MATLAB R2022a with Simulink** (the model was built and run in R2022a).
  Later releases normally open it as well, but the results reported in the
  manuscript were produced with R2022a.
* No additional toolbox is required (base MATLAB + Simulink only).
* No compilation step: the model is interpreted.

## 2. Layout

```
CSF_Model_Code/
├── README.md                 this file
├── run_all_scenarios.m       single entry point: runs everything in order
└── SimulinkModel/
    ├── CSF_Model.slx         the model: 9 volume-carrying compartments,
    │                         2 volume-less junction nodes and 3 constant-pressure
    │                         outlets (see Methods)
    ├── ABP_Data/             input arterial waveforms (radial pulse, Figshare
    │                          Pulse Wave Database), one file per scenario
    │   ├── 120.txt           120/80 mmHg  -> B, C-series, A3-series, S-series
    │   ├── 132.txt           132/84 mmHg  -> H_HM
    │   ├── 144.txt           144/88 mmHg  -> H_G1
    │   └── 164.txt           164/92 mmHg  -> H_G2, H_G2A, H+C, HA+C
    ├── MATLAB_script/
    │   ├── bp_compute.m          computes every model parameter for a scenario
    │   ├── Brain_Parameters.m    thin wrapper: writes bp_compute's output into
    │   │                         the base workspace, which is what the model reads
    │   ├── ABP_Get.m             reads one ABP_Data file, segments the cardiac
    │   │                         cycles, ensemble-averages and renormalises them
    │   ├── ABP_to_Central_GTF.m  radial -> central aortic waveform (12 harmonics,
    │   │                         200 Hz, Karamanoglu generalised transfer function)
    │   └── per_cen_TF.m          the transfer-function implementation used above
    └── Batch/
        ├── run_base_final.m      baseline B (must run first, see §3)
        ├── run_h_series.m        A1 hypertension
        ├── run_c_series.m        A2 communicating outflow impairment (14 scenarios)
        ├── run_cref_series.m     A2-ref single-outflow reference
        ├── run_a3_series.m       A3 combined pathology (H+C, HA+C)
        ├── run_s_series.m        Table 3 one-at-a-time sensitivity
        ├── base_ref.m            reads the B reference values back from
        │                         Results/H_series/B.mat (so no runner carries
        │                         hand-copied reference numbers)
        ├── collate_all.m         builds Results/collate_all.csv, one row per
        │                         scenario, from the saved .mat files
        ├── check_codes.m         checks all scenario codes against the
        │                         parameters they must produce (seconds, no simulation)
        ├── check_symbols.m       checks that the model's logged labels match the
        │                         symbol map below (seconds)
        ├── sym_rename_map.txt    manuscript symbol <-> code identifier map
        ├── compare_runs.m        difference table between any two saved runs
        ├── compare_a3_series.m   A3 interaction analysis (optional)
        ├── draw_baseline_timeseries.m   Fig. 3
        ├── draw_fig_hypertension.m      Fig. 4
        ├── draw_fig_comm_hydro.m        Fig. 5
        └── draw_fig_obst_hydro.m        Fig. 6
```

`Results/` is created automatically by the runners; it is not part of this
package because the per-scenario dumps are 100–350 MB each. Re-running the
included scripts reproduces them.

## 3. How to run

From the folder that contains `SimulinkModel/`:

```matlab
% everything, in the required order (baseline first)
run_all_scenarios
```

or step by step:

```matlab
matlab -batch "run_base_final"     % B          ~6 min
matlab -batch "run_h_series"       % A1         ~25 min
matlab -batch "run_c_series"       % A2         ~100 min
matlab -batch "run_cref_series"    % A2-ref     ~10 min
matlab -batch "run_a3_series"      % A3         ~50 min
matlab -batch "run_s_series"       % Table 3    ~110 min
matlab -batch "collate_all"        % master table
```

**The baseline must be run first.** Every other runner calls `base_ref.m`, which
reads the reference values from `Results/H_series/B.mat`; without it the runners
stop with an explicit error rather than silently using stale numbers.

Before a long batch, two checks take seconds each and catch the usual
mistakes (a broken parameter parser, a renamed signal):

```matlab
matlab -batch "check_codes"        % all scenario codes -> parameter values
matlab -batch "check_symbols"      % model's logged labels vs sym_rename_map.txt
```

### Numerical settings

* Solver `ode4` (fixed-step fourth-order Runge–Kutta), fixed step **0.01 s**.
* Initial conditions come from the `*0` and `*_free` parameters in `bp_compute.m`.
* **StopTime** is chosen per scenario, not fixed: each runner estimates the
  settling time from the outflow resistance and the compliance
  (`tau = 0.78 * R_eff * C` with `C = k_m / P_eq`) and takes `6*tau`, rounded up,
  with a floor of 3000 s. The estimate uses the *cranial* compliance only and
  therefore **under-estimates** the true settling time by roughly a factor of two
  (the spinal compartment carries 69 % of the craniospinal compliance); the
  runners print the estimate, and any scenario can be re-run with an explicit
  override, e.g.

  ```matlab
  matlab -batch "C_CODES={'C_OVS5'}; C_STOP_OVERRIDE=12000; run_c_series"
  ```

  The convergence of every run is reported as 12 block-means of the cycle-mean
  state plus the ratio of the last two differences (it should be < 1 and
  decaying).
* Results are always the **mean over the last 400 cardiac cycles** (0.8 s each),
  i.e. after the transient has decayed; pulsatile quantities are the
  peak-to-peak amplitude over the last 100 cycles.

## 4. Scenario codes

The codes are those of Table 2 of the manuscript.

| Group | Codes | Meaning |
|---|---|---|
| A0 | `B` | baseline, 120/80 mmHg, PVI 25 |
| A0-I | `B_I1.5`, `B_I25` | infusion tests (see note below) |
| A1 | `H_HM`, `H_G1`, `H_G2`, `H_G2A` | high-normal / grade 1 / grade 2 chronic hypertension, and acute grade 2 |
| A2 | `C_O…`, `C_V…`, `C_S…`, `C_OV…`, `C_OS…`, `C_VS…`, `C_OVS…` with suffix `2.5` or `5` | all fourteen combinations of the three outflow pathways at 250 % or 500 % of baseline resistance. `O` = olfactory lymphatic (`R_CP`), `V` = arachnoid granulation (`R_AG`), `S` = spinal lymphatic (`R_SPI`) |
| A2-ref | `C_Ref` | single-outflow reference: arachnoid granulations alone, at the unchanged total outflow resistance |
| A3 | `H+C`, `HA+C` | grade 2 hypertension (chronic / acute) + all three pathways at 500 % |
| Table 3 | `S_*` | one-at-a-time sensitivity, e.g. `S_PVI_30` (PVI 30 mL), `S_RBBB_4` (`R_BBB` ×4), `S_SWAP_AGMAX` (exchange `R_AG` and `R_CP`), `S_VALVE_OFF` / `S_VALVE_HIGH` (hold `R_IEG` at its low / high value for the whole run) |

A code is a `+`-joined list, so `H_G2+C_OVS5+Ven_E0.25` sets the vasculature
from the hypertension code, the outflow resistances from the `C_` code and the
ventricular elastance from the `Ven_` code. The Table 2 shorthands `H+C` and
`HA+C` are expanded inside `bp_compute.m`.

*The infusion scenarios (`B_I1.5`, `B_I25`) require the injection input of the
model to be driven; the corresponding runner is not part of this package.*

## 5. Symbol map

The manuscript's symbols and the code identifiers are not always spelled the
same way. `Batch/sym_rename_map.txt` is the complete table (`old`, `new`,
kind, note); the entries that matter most when reading the code are:

| Manuscript | Code | Meaning |
|---|---|---|
| $Q_{GS}$ | `Q_GS` | PAS → ECS flow across the IEG |
| $\overline{Q_{GS}}$ | `meanQ_GS` | its one-cycle moving average (Fig. 3e) |
| $Q_{ECS}$ | `Q_ECS` | ECS → PNS outflow |
| $Q_{OLS}$ | `Q_OLS` | cribriform plate → olfactory lymphatic system |
| $Q_{SSS}$ | `Q_SSS` | arachnoid granulations → superior sagittal sinus |
| $Q_{SLS}$ | `Q_SLS` | spinal nerve roots → spinal lymphatic system |
| $Q_{PCS}$ | `Q_PCS` | PAS → PVS, the paravascular route |
| $Q_{PNS}$ | `Q_PNS` | PNS → SASE |
| $P_{ECS}$ | `P_ECS` | intracranial pressure (ICP) |
| $V_T$ | `V_T` | total intracranial volume |
| — | `R_IEG` | the switched IEG resistance actually used, i.e. `R_IEGL` or `R_IEGH` |
| — | `R_FM` | craniospinal (foramen magnum) resistance between the ventricular outlet junction and the spinal SAS |

Quantities of the same name appear with identical labels in the model, in the
saved `.mat` datasets and in `collate_all.csv`, so a column of the master table
can be traced back to the model block that produces it.

## 6. Reproducing the figures

The figure scripts read the saved results, so the corresponding scenarios must
have been run first:

```matlab
matlab -batch "draw_baseline_timeseries"   % Fig. 3   (needs B)
matlab -batch "draw_fig_hypertension"      % Fig. 4   (needs the A1 group)
matlab -batch "draw_fig_comm_hydro"        % Fig. 5   (needs collate_all.csv)
matlab -batch "draw_fig_obst_hydro"        % Fig. 6   (needs collate_all.csv)
```

Figures are written next to the results they read (PDF, and PNG where the
script exports both).
