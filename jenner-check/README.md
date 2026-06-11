# Jenner compatibility tests

This directory was added by a pull request from the
[Jenner](https://jenneranalytics.com) project. Each `tNNN_*` subdirectory
contains a runnable bundle built from the SCORECI and PAIRBINCI macros in
this repository, paired with the published example datasets from
`tests/v_scoreci.sas` and `tests/v_pairbinci.sas`.

## What's in here

```
jenner-check/
├── README.md                       # this file
├── run_jenner.sh                   # mac/linux runner (curl)
├── run_jenner.bat                  # windows runner
├── run_jenner.sas                  # SAS-side runner
├── t001_scoreci_stratified_rd/     # stratified MN + SCAS risk difference (CMH-coherent)
├── t002_scoreci_unstratified_rd/   # Newcombe 1998 single-stratum examples
├── t003_scoreci_relative_risk/     # Gart & Nam 1988 unstratified RR
├── t004_scoreci_odds_ratio/        # Gart 1985 stratified OR + Gart-Nam RR (SCAS)
├── t005_scoreci_poisson_rates/     # cisapride meta-analysis, binomial + Poisson
└── t006_pairbinci_paired_proportions/  # Fagerland 2014 / Laud 2025 paired examples
```

Each bundle contains:

```
tNNN_*/
├── script.sas       # the macro + example calls, runnable end to end
├── autoexec.sas     # options applied before the script
├── expected.json    # fields a re-run is expected to reproduce
├── expected/        # captured log, listing and artifact links
└── meta.json        # provenance + the exact list of adaptations applied
```

`script.sas` embeds the macro source with a small set of mechanical,
semantics-preserving syntax adaptations so it runs on the Jenner engine;
`meta.json` in each bundle lists every adaptation. The confidence limits
reproduce the published values cited in the validation scripts
(for example the cisapride SCAS interval (0.245523, 0.370330) and the
Fagerland Table V Tango interval (-0.517, -0.026)).

## How to run

```bash
cd jenner-check/
./run_jenner.sh --all          # run every bundle against the Jenner API
./run_jenner.sh t001_scoreci_stratified_rd   # run one
```

No SAS installation, license, or local install is needed — the runner
POSTs each bundle to `https://api.jenneranalytics.com/v1/run` and prints
the status, log and listing. You can also paste any `script.sas` into the
hosted workspace at [jenneranalytics.com](https://jenneranalytics.com).
