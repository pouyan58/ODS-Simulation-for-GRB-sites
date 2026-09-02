# GRB ODS Spot Beam Simulator

## Purpose and scope

The GRB ODS Spot Beam Simulator is a MATLAB engineering model for estimating
inner and outer operational data sharing (ODS) angles around a protected
GOES Rebroadcast (GRB) receive-antenna boresight. It supplies twelve standard
cases: four receive locations, each pointed toward GOES East, GOES West, or
the configured GOES Backup longitude.

Angles are measured from the protected GRB antenna boresight toward the
instantaneous LEO satellite direction. They are not measured from zenith or
from the local horizon.

The model combines LEO orbit geometry, ten independently steered spot-beam
positions per satellite, carrier reuse and burst traffic, exact overlap with
the 10.9 MHz GRB passband, aggregate received interference, a GRB receive-
antenna pattern, atmospheric and clutter loss, command prediction, and
per-beam ODS actions. It searches for the smallest modeled inner/outer angle
pair meeting the configured interference-to-noise criterion.

This package is a public engineering sensitivity model. It is not a
coordination determination, spectrum authorization, site certification, or
representation of an operator's actual constellation, scheduler, traffic,
antenna mask, geographic exclusion policy, or ODS performance.

## Public baseline

The supplied baseline uses:

- a nominal 4.5 m GRB reflector, 120 K noise temperature, and 15.2 dB/K
  minimum G/T;
- 1686.6 MHz center frequency and 10.9 MHz necessary bandwidth;
- an ITU-R S.580-6/APEREC015 receive-pattern envelope with a -27 dBi
  absolute far-angle sensitivity assumption;
- a -6 dB aggregate I/N criterion at 99.9% modeled availability;
- 10 spot-beam steering positions per LEO satellite;
- a 60-degree satellite off-nadir field of regard;
- a 2-degree satellite-view half-power radius per spot;
- 30 dBW/MHz maximum boresight EIRP density per active spot;
- four contiguous 5 MHz carriers from 1675 to 1695 MHz;
- at most one active spot on each GRB-overlapping carrier per satellite;
- aggregation of every eligible scheduled satellite and beam in linear
  power;
- a 30 km geographic avoidance radius at Fairmont and Wallops;
- a 50 km geographic avoidance radius at Sioux Falls and Suitland;
- outer-zone retasking to a feasible beam center at least 180 km from the
  protected receiver;
- protected-carrier shutdown inside the inner ODS angle;
- geometric LEO visibility down to the local horizon; and
- no minimum operational outer-angle floor.

The 60-degree value is a steering field, not an individual spot beamwidth.
The ten stored steering positions form a snapshot lattice and do not assert
continuous tiling of the complete field of regard.

## Requirements

- MATLAB R2025b or a compatible later release.
- Satellite Communications Toolbox for satellite-scenario and Walker-
  constellation functions.
- The four supplied site-specific ITU-R P.618 lookup tables in `data/`.
- MathWorks ITU digital maps only when regenerating P.618 tables.

Normal simulations use the included atmospheric tables, need no internet
connection, and do not need the ITU digital-map files.

## Quick start

Add the extracted package folder to the MATLAB path:

```matlab
packageFolder = "C:\path\to\GRB_ODS_Spot_Beam_Simulator";
addpath(packageFolder)
```

Display the twelve choices and run one case interactively:

```matlab
results = runGrbOdsSpotBeamSimulator();
```

Run a case by menu number, `SiteId`, or combined label:

```matlab
results = runGrbOdsSpotBeamSimulator(10);
results = runGrbOdsSpotBeamSimulator("wallops_va_goes_east");
results = runGrbOdsSpotBeamSimulator("Wallops, VA | GOES East");
```

The selected primary angles are available as:

```matlab
outerAngleDeg = results.practical.outerAngleDeg;
innerAngleDeg = results.practical.innerAngleDeg;
```

List cases without running the simulation:

```matlab
sites = listGrbOdsSpotBeamCases();
```

Run all twelve public cases and reuse matching saved results:

```matlab
summary = runAllGrbOdsSpotBeamCases(ResumeCompleted=true);
```

Set `ResumeCompleted=false` to force recalculation. A saved result is reused
only when its release ID, site, ZA radius, duration, time step, Monte Carlo
count, and -27 dBi far-angle setting match the requested configuration.

## Supplied site and boresight catalog

The simulator reads `data/GRB_Receiver_Sites.csv`.

| No. | Site | Latitude | Longitude | GOES boresight | Feed height | Default ZA |
|---:|---|---:|---:|---|---:|---:|
| 1 | Fairmont, WV | 39.43361 | -80.1928 | East | 5.0 m | 30 km |
| 2 | Fairmont, WV | 39.43361 | -80.1928 | West | 5.0 m | 30 km |
| 3 | Fairmont, WV | 39.43361 | -80.1928 | Backup | 5.0 m | 30 km |
| 4 | Sioux Falls, SD | 43.73513 | -96.6255 | West | 7.4 m | 50 km |
| 5 | Sioux Falls, SD | 43.73513 | -96.6255 | Backup | 7.4 m | 50 km |
| 6 | Sioux Falls, SD | 43.73513 | -96.6255 | East | 7.4 m | 50 km |
| 7 | Suitland, MD | 38.85217 | -76.9367 | East | 24.0 m | 50 km |
| 8 | Suitland, MD | 38.85217 | -76.9367 | West | 24.0 m | 50 km |
| 9 | Suitland, MD | 38.85217 | -76.9367 | Backup | 24.0 m | 50 km |
| 10 | Wallops, VA | 37.94651 | -75.4621 | East | 5.0 m | 30 km |
| 11 | Wallops, VA | 37.94651 | -75.4621 | West | 5.0 m | 30 km |
| 12 | Wallops, VA | 37.94651 | -75.4621 | Backup | 5.0 m | 30 km |

GOES East is modeled at 75.2 degrees W, GOES West at 137.0 degrees W, and
GOES Backup at 104.7 degrees W. The local antenna-boresight vector is
calculated from the site coordinates and assigned GEO longitude.

Feed height is height above local ground, not surveyed geodetic height. The
baseline also uses it as a geodetic-height proxy. Supply `GeodeticHeightM`
when a surveyed antenna-reference height is available, and regenerate the
P.618 lookup when the propagation geometry changes materially.

Required custom-catalog columns are `Site`, `LatitudeDeg`, `LongitudeDeg`,
`Satellite`, and `FeedHeightM`.

## Verified all-case results

Every supplied result uses two simulated hours, 2-second samples, four
Monte Carlo realizations, and the public assumptions described here.

| Site | Boresight | ZA | Primary outer | Primary inner | Robust availability |
|---|---|---:|---:|---:|---:|
| Fairmont, WV | GOES East | 30 km | 40.5 deg | 1.0 deg | 99.917% |
| Fairmont, WV | GOES West | 30 km | 39.0 deg | 2.0 deg | 99.917% |
| Fairmont, WV | GOES Backup | 30 km | 27.5 deg | 4.0 deg | 99.917% |
| Sioux Falls, SD | GOES West | 50 km | 41.5 deg | 8.0 deg | 99.917% |
| Sioux Falls, SD | GOES Backup | 50 km | 45.0 deg | 5.0 deg | 99.917% |
| Sioux Falls, SD | GOES East | 50 km | 39.0 deg | 3.0 deg | 99.917% |
| Suitland, MD | GOES East | 50 km | 46.0 deg | 5.0 deg | 99.917% |
| Suitland, MD | GOES West | 50 km | 44.5 deg | 6.0 deg | 99.917% |
| Suitland, MD | GOES Backup | 50 km | 34.5 deg | 20.5 deg | 99.944% |
| Wallops, VA | GOES East | 30 km | 34.0 deg | 3.5 deg | 99.917% |
| Wallops, VA | GOES West | 30 km | 41.0 deg | 2.0 deg | 99.917% |
| Wallops, VA | GOES Backup | 30 km | 27.5 deg | 9.5 deg | 99.917% |

The complete machine-readable table is
`grb_ods_spot_beam_all_12_results.csv`. Each result also retains a stricter
zero-modeled-exceedance diagnostic. That diagnostic is not the primary
public answer.

These values are conditional model outputs, not recommended operational
limits. The Suitland/GOES Backup inner angle is especially sensitive to
boresight geometry, scheduler states, aggregate traffic, and the assumed
retasking behavior.

## GRB receiver and antenna model

| Parameter | Value |
|---|---:|
| Center frequency | 1686.6 MHz |
| Necessary bandwidth | 10.9 MHz |
| Nominal reflector diameter | 4.5 m |
| System noise temperature | 120 K |
| Minimum G/T | 15.2 dB/K |
| Derived peak gain | approximately 35.99 dBi |
| Derived aperture efficiency | approximately 0.6282 |
| Primary protection criterion | -6 dB I/N at 99.9% availability |
| Polarization coupling | worst-case co-polar |
| Implementation loss | 1 dB |

Thermal noise at the antenna-terminal reference plane is:

```text
N_W   = k T B
N_dBm = 10 log10(N_W) + 30
```

The model protects one worst-case co-polar receiver chain. It does not add
independent interference power across the two GRB circular-polarization
chains.

The receive-pattern envelope follows ITU-R S.580-6 with the attached
APEREC015 small-aperture extension. Peak gain is derived from G/T and noise
temperature. Beyond the S.580 far-angle breakpoint, the public baseline
uses -27 dBi absolute receive gain. This is a sensitivity assumption, not a
measured GRB site pattern. Measured azimuth/elevation patterns should replace
it for operational work.

## Constellation, spot beams, spectrum, and traffic

| Shell | Altitude | Inclination | Satellites | Planes | Phasing |
|---|---:|---:|---:|---:|---:|
| DTC-53 | 340 km | 53 deg | 325 | 13 | 1 |
| DTC-43 | 355 km | 43 deg | 325 | 13 | 1 |

Every satellite has ten candidate steering positions inside a 60-degree
off-nadir field of regard. The lattice is stored in
`evidence/beam_packing_layout.csv` and scaled to each shell altitude. Each
spot has a separately configured 2-degree satellite-view half-power radius.

The per-spot transmit pattern uses the ITU-R S.1528-0 LEO reference
envelope. Scan loss follows projected-aperture cosine loss, with modeled
power control maintaining the same 30 dBW/MHz peak EIRP-density ceiling at
every steering position. An optional off-axis-mask credit is configurable
but is zero in the public baseline.

Four contiguous 5 MHz carriers cover 1675-1695 MHz. Static four-color reuse
maps the ten steering positions cyclically across those carriers. The GRB
passband spans 1681.15-1692.05 MHz, producing these rectangular overlaps:

| Carrier | Nominal range | GRB overlap |
|---:|---|---:|
| 1 | 1675-1680 MHz | 0 MHz |
| 2 | 1680-1685 MHz | 3.85 MHz |
| 3 | 1685-1690 MHz | 5.00 MHz |
| 4 | 1690-1695 MHz | 2.05 MHz |

Each beam's received interference is integrated over its own overlap width.
The scheduler permits at most one active spot on each overlapping carrier
per satellite, so up to three beams from one satellite may contribute at the
same time.

Each carrier uses an independent two-state burst scheduler with 20% mean
duty cycle, 30-second mean on time, 120-second mean off time, 60% mean
scheduled loading, and bounded loading variation. Every visible satellite
remains a potential contributor through its active spot's main lobe or
sidelobes. The aggregate is not limited to a fixed number of closest or
strongest satellites.

## Geographic avoidance and ODS actions

No active overlapping-carrier beam center or served user is permitted
inside the configured geographic avoidance radius. The standard site
mapping is 30 km for Fairmont and Wallops and 50 km for Sioux Falls and
Suitland.

ODS state is based on angular separation from the protected GRB boresight:

- outside the outer angle: normal scheduled operation;
- between inner and outer angles: overlapping-carrier spots are physically
  retasked to a feasible beam center at least 180 km from the protected
  receiver, and gain toward the receiver is recalculated; and
- inside the inner angle: the affected carrier is shut down, contributing
  zero power to the protected receiver in the numerical model.

The model applies 15 seconds of command latency plus 1 second of clock
allowance. It explicitly predicts satellite direction 16 seconds ahead and
uses the smaller of current and predicted separation for entry. Measured
separation controls exit. A 0.5-degree hysteresis limits rapid state changes.

Change ZA at runtime without editing source code:

```matlab
results = runGrbOdsSpotBeamSimulator( ...
    "wallops_va_goes_east",ZaRadiusKm=40);
```

A nonstandard run is written to a ZA-tagged folder such as
`results_wallops_va_goes_east_za40km`, preserving the supplied public result.
ZA must remain below the 150 km calibration endpoint and the 180 km outer-
retasking distance. `configureGrbOdsZaRadius` updates every linked field.

## Propagation and interference calculation

### Free-space spreading and PFD cap

For each active spot, boresight EIRP density is reduced by transmit-pattern
offset, scheduled loading, and bounded uncertainty. Free-space spreading is
applied once as power-flux density:

```text
PFD_dBW/m2/MHz = EIRP_dBW/MHz - 10 log10(4 pi R^2)
```

The configured -80 dBW/m2/MHz PFD cap is applied independently to each
beam after range loss.

### Receive effective area

For receive gain `G_rx` at the LEO-to-GRB-boresight separation:

```text
A_e,dB(m2) = G_rx,dBi + 20 log10(lambda) - 10 log10(4 pi)
```

### Atmospheric attenuation

The supplied lookup tables were generated at 1686.6 MHz with the nominal
4.5 m/0.6282-efficiency antenna using ITU-R P.618 through MATLAB Satellite
Communications Toolbox and a 5% annual exceedance percentage. Gaseous,
cloud, rain, and scintillation contributions are combined in
`TotalAtmosphericLossDb` and interpolated by elevation with PCHIP. Tables
cover 5 through 90 degrees. Visible links from 0 through 5 degrees use the
5-degree value.

### Local clutter

The model uses the ITU-R P.2108 open/rural height-gain form with a 10 m
representative clutter height. For feed height below 10 m:

```text
L_clutter = -(21.8 + 6.2 log10(f_GHz)) log10(h_rx / 10 m)
```

Clutter loss is zero at or above the representative height. Feed height,
terrain class, and clutter height are provisional inputs, not site-survey
results.

### Received interference and aggregation

The implemented per-beam relationship is equivalent to:

```text
P_r,dBm = PFD_dBW/m2/MHz
          + 10 log10(overlap bandwidth / 1 MHz)
          + A_e,dB(m2)
          + 30
          + power uncertainty
          + 10 log10(scheduled loading)
          - ODS action
          - polarization loss
          - implementation loss
          - atmospheric loss
          - clutter loss
```

All eligible post-action powers are converted to milliwatts and summed:

```text
I_aggregate,mW = sum(10^(P_r,dBm/10))
I/N_dB         = 10 log10(I_aggregate,mW) - N_dBm
```

The provisional NRAO-informed calibration applies up to -18 dB additional
distance-ramped satellite sidelobe correction between the configured ZA
boundary and 150 km. It transfers 1990-1995 MHz test evidence to the
hypothetical 1675-1695 MHz case and must be replaced by band- and operator-
specific measurements for operational use.

## Monte Carlo criteria and angle search

| Parameter | Value |
|---|---:|
| Start time | 2026-08-01 00:00 UTC |
| Duration | 2 hours |
| Time step | 2 seconds |
| Samples per realization | 3,601 |
| Monte Carlo realizations | 4 |
| Random seed | 16751695 |
| Link-power sigma | 3 dB |
| Ephemeris-angle sigma | 0.05 deg |
| Receiver-pointing sigma | 0.10 deg |
| Scheduled-load sigma | 0.15 |
| Across-run robustness percentile | 95th |

Positive link-power uncertainty is clipped at the configured 30 dBW/MHz
beam ceiling. Angular errors are clamped to visible-sky geometry.

For each candidate pair and realization, the model calculates the fraction
of time samples above -6 dB I/N and the realization's 99.9th-percentile
I/N. Across four realizations, it uses the 95th percentile of both metrics.
Primary feasibility requires:

```text
robust exceedance <= 0.1 percent
robust 99.9th-percentile I/N <= -6 dB
```

Thus the reported robust availability is at least 99.9%. A separate strict
diagnostic requires zero modeled samples above -6 dB in every realization.
Neither interpretation proves performance over every orbit, traffic state,
uncertainty draw, or year.

Only pairs with `inner <= outer` are evaluated. Search stages are:

1. coarse candidates at 0, 0.5, 1, 2, 3, and 5-degree increments through
   145 degrees, clipped to visible-sky separation;
2. 1-degree refinement within 7 degrees of coarse strict and primary
   candidates; and
3. 0.5-degree refinement within 1.5 degrees of the refined candidates.

Feasible pairs are ranked by `outer + inner`; ties prefer smaller outer and
then smaller inner angle. Service impact is reported but is not an added
optimization weight in the public baseline.

## Runtime options

`runGrbOdsSpotBeamSimulator` accepts:

| Option | Default | Purpose |
|---|---|---|
| `CatalogFile` | supplied CSV | Use another compatible site catalog |
| `OutputRoot` | package folder | Isolate results under another root |
| `GeodeticHeightM` | feed-height proxy | Supply surveyed reference height |
| `ZaRadiusKm` | site mapping | Override geographic avoidance radius |
| `DurationSeconds` | 7200 | Change simulated duration |
| `TimeStepSec` | 2 | Change time resolution |
| `MonteCarloRuns` | 4 | Change realization count |
| `WriteDetailedTimeSeries` | `true` | Write selected time-series CSV |
| `RegenerateAtmosphericLookup` | `false` | Rebuild the P.618 table |
| `ItuDigitalMapsFolder` | empty | Map folder required for regeneration |
| `DryRun` | `false` | Validate and display configuration only |

Example short smoke test in a separate output folder:

```matlab
results = runGrbOdsSpotBeamSimulator( ...
    "wallops_va_goes_east", ...
    OutputRoot="C:\temp\grb_ods_smoke", ...
    DurationSeconds=300,TimeStepSec=10,MonteCarloRuns=1, ...
    WriteDetailedTimeSeries=false);
```

Run all cases with one common ZA radius:

```matlab
summary = runAllGrbOdsSpotBeamCases(ZaRadiusKm=40);
```

Or provide one radius per catalog row:

```matlab
zaByCase = [30 30 30 50 50 50 50 50 50 30 30 30];
summary = runAllGrbOdsSpotBeamCases(ZaRadiusKm=zaByCase);
```

## Changing configuration

The central public configuration is
`config/defaultGrbOdsSpotBeamConfig.m`.

| Area | Fields or source |
|---|---|
| Site, feed height, default ZA | `data/GRB_Receiver_Sites.csv` |
| GEO longitude mapping | `loadGrbReceiverSites.m` |
| GRB receiver and criterion | `cfg.receiver.*` |
| Receive antenna pattern | `cfg.antenna.*` |
| LEO shells | `cfg.constellation.shells` |
| Spot geometry and S.1528 pattern | `cfg.beams.*` and packing CSV |
| Frequency overlap and traffic | `cfg.resources.*` |
| EIRP, PFD cap, and fixed losses | `cfg.rf.*` |
| Geographic zone | runtime `ZaRadiusKm` or `configureGrbOdsZaRadius` |
| P.618 and P.2108 propagation | `cfg.propagation.*` |
| Aggregate contributor selection | `cfg.aggregation.*` |
| Retasking, shutdown, and prediction | `cfg.ods.*` |
| NRAO-informed transfer calibration | `cfg.calibration.*` |
| Monte Carlo uncertainty | `cfg.uncertainty.*` |
| Duration, time step, and seed | `cfg.simulation.*` |
| Angle grids and ranking | `cfg.search.*` |

Regenerate the P.618 table after changing site coordinates, receiver
frequency, antenna diameter, aperture efficiency, geodetic height, or
annual percentage. A different beam count requires a matching layout CSV
and intentional changes to reuse mapping, validation, documentation, and
verification. A new GOES assignment name requires a longitude mapping in
`loadGrbReceiverSites.m`.

Use a separate `OutputRoot` for sensitivity studies to avoid replacing the
supplied public results. When changing linked traffic fields, update their
compatibility aliases and derived fields consistently.

## Output files

Each standard case writes to `results_<site_id>`:

- `ods_results.mat`: complete configuration, geometry summary, traffic
  summary, search table, primary/strict metrics, and representative series;
- `configuration.json`: resolved configuration;
- `grb_ods_angle_summary.csv`: primary and strict case summary;
- `primary_selected_angles.csv`: selected 99.9%-availability metrics;
- `practical_selected_angles.csv`: engine-level primary metrics;
- `selected_angles.csv`: strict zero-exceedance diagnostic;
- `selected_time_series.csv`: aggregate I/N and contributor counts;
- `angle_search.csv` and `angle_search.png`: candidate search;
- `aggregate_in_over_n.png` and its CDF: strict diagnostic plots; and
- `practical_aggregate_in_over_n.png` and its CDF: primary result plots.

The all-case runner writes `grb_ods_spot_beam_all_12_results.csv` in the
package root.

## Verification

Run the public-package checks after extraction:

```matlab
verification = verifyGrbOdsSpotBeamSimulator();
```

The verifier checks all catalog/configuration mappings, default and custom
ZA behavior, the GRB receiver and three-carrier overlap, the fixed ten-beam
assumptions, the four atmospheric tables, saved metadata, required files,
angle validity, and the -6 dB/99.9%-availability criterion.

## Included evidence and references

- `evidence/grb_receiver_profile.csv`
- `evidence/beam_packing_layout.csv`
- `evidence/nrao_dtc_calibration.csv`
- `references/R-REC-S.580-6-200401.pdf`
- `references/APEREC015V01.pdf`

See `NOTICE.md` for source attribution, model-status statements, and
redistribution cautions. Users must independently validate constellation,
traffic, antenna, propagation, receiver, regulatory, and ODS inputs before
using any result for planning or coordination.
