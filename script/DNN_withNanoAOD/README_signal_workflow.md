# DNN_withNanoAOD Workflow

This directory contains two separate model lines:

1. `ttbb` B-jet assignment
2. signal classification

Do not mix them. The BJA files are kept as-is and are not part of the simplified signal workflow below.

## Current Contract

For the simplified signal workflow, the main chain is now:

1. extract raw per-sample features
2. mix datasets
3. assign labels during mix only
4. reweight during mix if requested
5. validate the mixed dataset contract before training/testing

This means:

- `signalExtractionInputFromNanoAOD.py` stays classification-agnostic
- `extractSignalFromSampleJson.py` only dispatches extraction jobs
- `prepareSignalDatasets.py` is the only place that assigns `label`
- `simpleSignalModel.py train/test` validates the mixed dataset structure against the mix report before use

ROOT deployment scripts may still expose `--class-schema` for score naming or checkpoint compatibility, but the extract -> mix -> train/test/analyze data contract is now centralized in the mix stage.

## Signal Files

Main workflow:

- `signalExtractionInputFromNanoAOD.py`
  ROOT -> single-sample signal input `.npy`

- `prepareSignalDatasets.py`
  many per-sample `.npy` -> `OS_train.npy`, `OS_val.npy`, `OS_test.npy`, `SRtest.npy`
  Main option:
  - `--weighting-mode balanced` : recommended main workflow
  - `--weighting-mode raw` : split only, keep original weights
  - `--class-schema ...` : choose how channels are mapped into training categories during mix
  Main responsibilities:
  - assign labels from `channelId`
  - build `channel_map.json`
  - write `mix_report.json`
  - validate the feature / field naming contract for the mixed outputs

- `simpleSignalModel.py`
  - `train` : train from `OS_train.npy` and `OS_val.npy`
  - `test` : run prediction on `OS_test.npy` or `SRtest.npy`
  Main responsibilities:
  - validate dataset field names, feature names, label-channel consistency, and OS/SS split expectations
  - save model checkpoint with the dataset contract embedded

- `plotSignalScoreDistributions.py`
  dataset + prediction + channel map -> score plots

- `plotInputFeatureDistributions.py`
  dataset + channel map -> feature shape plots

- `analyzeSimpleSignalFeatures.py`
  dataset + model + channel map -> permutation / gradient / optional covariance analysis

ROOT deployment branch:

- `predictSignalRootToNpy.py`
  ROOT -> signal prediction `.npy`
  This is the preferred name for ROOT inference.

- `writeSignalPredictionToRoot.py`
  ROOT + signal prediction `.npy` -> new ROOT with signal branches
  This is the preferred name for writing signal scores back to ROOT.

BJA deployment branch:

- `predictBJARootToNpy.py`
  ROOT -> BJA prediction `.npy`

- `writeBJAPredictionToRoot.py`
  ROOT + BJA prediction `.npy` -> new ROOT with BJA branches

Schema definition:

- `signal_class_config.py`

## Default Label Configuration

The mix stage uses the `--class-schema` selected in `prepareSignalDatasets.py`.
If not provided, it defaults to the default schema from `signal_class_config.py`.

The current default schema is:

```bash
three_class_omit_ttz_merged_tttt
```

Meaning:

- `TTHH`
  - `TTHH_DL`
  - `TTHH_SL`
- `tt_b`
  - `ttbar` = `ttbarDL + ttbarSL`
  - `TTBB` = `TTBB_DL + TTBB_SL`
  - `TTHBB`
- `ttX_like`
  - `TTW`
  - `TTHnonBB`
  - `TTTT`
- `TTZ`
  - removed from all OS splits
  - kept only in `SRtest`

For balancing inside the default three-class setup, the effective balance groups are:

- category `0`: `TTHH_DL`, `TTHH_SL`
- category `1`: `ttbar`, `TTBB`, `TTHBB`
- category `2`: `ttX_like_merged = TTW + TTHnonBB + TTTT`

Your requested four-class setup is still available as:

```bash
four_class_merged
```

Meaning:

- `TTHH`
  - `TTHH_DL`
  - `TTHH_SL`
- `tt_b`
  - `ttbar` = `ttbarDL + ttbarSL`
  - `TTBB` = `TTBB_DL + TTBB_SL`
  - `TTHBB`
- `ttX_like`
  - `TTW`
  - `TTZ` = `TTZ_low + TTZ_high`
  - `TTHnonBB`
- `TTTT`

For balancing inside each category, the effective balance groups are:

- category `0`: `TTHH_DL`, `TTHH_SL`
- category `1`: `ttbar`, `TTBB`, `TTHBB`
- category `2`: `TTW`, `TTZ`, `TTHnonBB`
- category `3`: `TTTT`

The underlying raw channels inside merged groups such as `ttbarDL/ttbarSL` and `TTZ_low/TTZ_high` do not become separate training labels and are not balanced separately.

Additional schema now available:

```bash
four_class_omit_ttz
```

Meaning:

- same four training categories as `four_class_merged`
- `TTZ_low + TTZ_high -> TTZ -> ttX_like`
- but `TTZ` is excluded from all OS splits:
  - not in `OS_train`
  - not in `OS_val`
  - not in `OS_test`
- `TTZ` is still kept in `SRtest`

This is useful when you want to suppress the OS-specific `TTZ -> Z->mumu` behavior from shaping the ttX training category, while still retaining TTZ in the same-sign validation region.

Additional binary sub-training schemas:

```bash
binary_tthh_vs_ttb
```

- `TTHH`
- `tt_b = ttbar + TTBB + TTHBB`

```bash
binary_tthh_vs_ttw_tthnonbb
```

- `TTHH`
- `ttW_like = TTW + TTHnonBB`

```bash
binary_tthh_vs_tttt
```

- `TTHH`
- `TTTT`

When using one of these reduced schemas, you can still pass the full extracted input directory to `prepareSignalDatasets.py`.
Channels that are not part of the selected schema will be skipped automatically and recorded in `mix_report.json`.

Preferred approach for sub-training:

- keep one main mixed dataset
- project to a sub-task only at `train/test/analyze` time with `--label-subset`

Examples:

```bash
--label-subset 0 1
```

- `TTHH` vs `tt_b`

```bash
--label-subset 0 2
```

- `TTHH` vs `ttX_like`

```bash
--label-subset 0 3
```

- `TTHH` vs `TTTT`

This avoids regenerating many redundant dataset folders when the mixed dataset itself does not change.

You can also sub-train on reduced input groups without regenerating datasets.
`simpleSignalModel.py train/test` accepts:

```bash
--feature-groups 1 2 6 7
```

Group mapping:

- `1`: muon four-momentum
- `2`: jet four-momentum
- `3`: jet-pair `dr` + invariant mass
- `4`: dimuon `dr` + mass
- `5`: muon-jet `dr`
- `6`: extra muon information (`miniIso`, `jetRelIso`, `jetDF`, `promptMVA`)
- `7`: extra jet information (`btag`, jet/bjet multiplicities and summed `pt`, jet centrality)

For batch ablation training, use:

```bash
python3 runSimpleSignalAblations.py \
  --train signal_split/OS_train.npy \
  --val signal_split/OS_val.npy \
  --test signal_split/OS_test.npy \
  --dataset-report signal_split/mix_report.json \
  --output-root ablation_runs \
  --label-subset 0 1
```

This runs the full batch pipeline for each preset:

- train the masked-input model
- test it and save `predictions.npy`, confusion matrix, and globally normalized score plots
- run `analyzeSimpleSignalFeatures.py`
- build a combined comparison package under `ablation_runs/comparison`

By default, the batch runner uses `--score-normalization global`, which keeps the real event weights and only applies one overall normalization to the plotted sample. If you want raw weighted yields, use `--score-normalization raw`. If you want per-class shape comparison with the true signed weights, use `--score-normalization shape`. If you specifically want a pure visualization-oriented unit-area plot that ignores weight signs, use `--score-normalization shape_abs`.

## Minimal Signal Workflow

### 1. Extract per-sample inputs

```bash
python3 signalExtractionInputFromNanoAOD.py \
  --input sample.root \
  --output sample_signal.npy \
  --era RunIII2024Summer24NanoAODv15 \
  --sample-name TTHH_DL_2B2W_batch1
```

Output fields at this stage are only:

- all extracted feature columns
- `weight`
- `isSameSign`

No classification label is written here.

Batch wrapper:

```bash
python3 extractSignalFromSampleJson.py \
  --sample-json source_cleanup/json/samples/SR_medium_muon/RunIII2024Summer24NanoAODv15_corrected_temp.json \
  --output-dir extracted_signal_inputs \
  --era RunIII2024Summer24NanoAODv15
```

### 2. Build mixed datasets

Recommended main workflow:

```bash
python3 prepareSignalDatasets.py \
  --inputs sample1.npy sample2.npy sample3.npy \
  --output-dir signal_split \
  --weighting-mode balanced \
  --class-schema four_class_merged
```

If you want to remove TTZ from all OS splits while keeping it in `SRtest`, use:

```bash
python3 prepareSignalDatasets.py \
  --inputs sample1.npy sample2.npy sample3.npy \
  --output-dir signal_split \
  --weighting-mode balanced \
  --class-schema four_class_omit_ttz
```

`--inputs` now accepts either:

- explicit `.npy` files
- directories containing `.npy` files
- parent directories whose subdirectories contain per-channel `.npy` files

For example, this now works directly:

```bash
python3 prepareSignalDatasets.py \
  --inputs threeBJet_signalExtraction_input/* \
  --output-dir four_class_split_dataset \
  --weighting-mode balanced
```

If the shell expands `*` into channel directories such as `TTHH_DL/`, `ttbarSL/`, etc., the script will recursively collect the `.npy` files inside them.

Main outputs:

- `signal_split/OS_train.npy`
- `signal_split/OS_val.npy`
- `signal_split/OS_test.npy`
- `signal_split/SRtest.npy`
- `signal_split/channel_map.json`
- `signal_split/mix_report.json`

`mix_report.json` records:

- one shared `common_input_fields` block
- input file list
- per-input field override only if a file differs from the common field layout
- per-input event counts
- OS / SS counts
- input weight sums
- channel -> canonical channel -> class -> label mapping
- split summaries
- reweighting statistics before / after balancing

Optional debug files are only written if you pass:

```bash
--save-debug-files
```

### 3. Train

```bash
python3 simpleSignalModel.py train \
  --train signal_split/OS_train.npy \
  --val signal_split/OS_val.npy \
  --output-dir train_run \
  --training-weight-mode input \
  --dataset-report signal_split/mix_report.json
```

Binary sub-training on top of the same mixed dataset:

```bash
python3 simpleSignalModel.py train \
  --train signal_split/OS_train.npy \
  --val signal_split/OS_val.npy \
  --output-dir train_run_tthh_vs_ttb \
  --training-weight-mode input \
  --dataset-report signal_split/mix_report.json \
  --label-subset 0 1
```

Training automatically resolves `signal_split/mix_report.json` and validates:

- feature names
- mixed field names
- label range
- `channelId -> label` consistency
- OS-only expectation for train / val

If the report is not next to the dataset, pass it explicitly:

```bash
--dataset-report path/to/mix_report.json
```

Main outputs:

- `train_run/model.pt`
- `train_run/loss.png`

Optional extra:

```bash
--save-weight-summary
```

### 4. Test on dataset

OS:

```bash
python3 simpleSignalModel.py test \
  --test signal_split/OS_test.npy \
  --model train_run/model.pt \
  --output-dir test_os
```

SS:

```bash
python3 simpleSignalModel.py test \
  --test signal_split/SRtest.npy \
  --model train_run/model.pt \
  --output-dir test_sr
```

Testing also validates the dataset contract against:

- `mix_report.json`
- the contract embedded inside `model.pt`

Main outputs:

- `test_os/predictions.npy`
- `test_os/confusion_matrix_*.png`

Optional basic score plots:

```bash
--save-basic-plots
```

### 5. Plot detailed score distributions

```bash
python3 plotSignalScoreDistributions.py \
  --dataset signal_split/OS_test.npy \
  --predictions test_os/predictions.npy \
  --channel-map signal_split/channel_map.json \
  --output-dir score_plots_os \
  --class-schema three_class
```

Use this when you want:

- category-wise score comparison
- channel-wise score comparison
- category-internal channel breakdown

For a sub-trained model, the plotting task is inferred automatically from `predictions.npy`.
You can also force it explicitly, for example:

```bash
--label-subset 0 1
```

### 6. Plot input feature shapes

```bash
python3 plotInputFeatureDistributions.py \
  --dataset signal_split/OS_train.npy \
  --channel-map signal_split/channel_map.json \
  --output-dir feature_plots_train \
  --class-schema three_class
```

Use this when you want:

- category feature shape comparison
- channel breakdown inside each category

For a sub-trained task, add:

```bash
--dataset-report signal_split/mix_report.json --label-subset 0 1
```

Additional diagnostic for one specific channel:

```bash
python3 plotInputFeatureDistributions.py \
  --dataset signal_split/OS_test.npy \
  --channel-map signal_split/channel_map.json \
  --predictions test_os/predictions.npy \
  --focus-channel TTHH_DL \
  --group-by prediction \
  --output-dir feature_plots_os_by_prediction \
  --class-schema four_class_merged
```

This produces per-feature plots inside one chosen channel, with separate histograms for events predicted as:

- signal
- ttbar / tt_b
- ttX
- TTTT

This is useful for diagnostics such as:

- `TTHH_DL` events predicted as `TTHH`
- `TTHH_DL` events predicted as `tt_b`
- `TTHH_DL` events predicted as `ttX_like`
- `TTHH_DL` events predicted as `TTTT`

This also works for sub-trained models. The script uses the prediction file's score fields to infer the active task if you do not pass `--label-subset`.

### 7. Analyze feature importance

```bash
python3 analyzeSimpleSignalFeatures.py \
  --dataset signal_split/OS_val.npy \
  --model train_run/model.pt \
  --output-dir feature_analysis \
  --dataset-report signal_split/mix_report.json
```

For a sub-trained model, you can use the same dataset and either:

- omit `--label-subset`, in which case the analysis follows the task stored in the checkpoint
- or pass the same subset explicitly, for example:

```bash
--label-subset 0 1
```

Default outputs include:

- `gradients/gradient_abs_by_target.png`
- `gradients/gradient_importance_summary.csv`
- `permutation/permutation_importance.png`
- `permutation/permutation_importance_summary.csv`
- `group_permutation/group_permutation_importance.png`
- `group_permutation/group_permutation_importance_summary.csv`
- `channel_permutation/*.png`
- `channel_permutation/channel_permutation_heatmap_normalized.png`
- `covariance/covariance_overall.png`
- `covariance/correlation_overall.png`
- `covariance/covariance_<class>.png`
- `covariance/correlation_<class>.png`

Optional extras:

```bash
--gradient-views full
```

The grouped permutation view partitions the 61 inputs into 7 physics-motivated blocks:

- muon four-momentum
- jet four-momentum
- jet-pair `dr` + invariant mass
- dimuon `dr` + mass
- muon-jet `dr`
- extra muon information (`miniIso`, `jetRelIso`, `jetDF`, `promptMVA`)
- extra jet information (`btag`, jet/bjet multiplicities and summed `pt`, jet centrality)

If you already trained and analyzed multiple runs and want only the aggregate comparison package for score shapes, gradients, and permutation studies, use:

```bash
python source_cleanup/script/DNN_withNanoAOD/compareSimpleSignalRuns.py \
  --runs runs_01_full runs_01_g1267 runs_01_g126 runs_01_g127 \
  --labels full g1267 g126 g127 \
  --output-dir compare_01 \
  --top-n 15
```

This writes:

- `run_summary.csv` with unweighted and `abs(weight)`-weighted accuracy
- `score_summary.csv` with per-truth-class score means for each run
- `scores/*.png` with overlaid score histograms across runs
- `gradient_comparison.csv` and `gradient_top_features.png`
- `permutation_comparison.csv` and `permutation_top_features.png`
- `group_permutation_comparison.csv` and `group_permutation_all_groups.png`

## Gradient Importance Meaning

The script uses gradient of the model output node with respect to the input features.
In practice this is the gradient of the selected logit, not the softmax probability.

Main views:

- `by_target`
  - fixed target node
  - example: for every event, compute `|d logit_TTTT / d input|`
  - this is the main gradient view to use

- `overall`
  - for each event, use the model's own predicted node
  - example: if the event is predicted as `TTTT`, use `|d logit_TTTT / d input|`
  - this describes what the model uses in its actual decisions

- `true_class`
  - restrict to events whose true label is a given class
  - then use that class node
  - example: on true `TTTT` events, compute `|d logit_TTTT / d input|`

If you only want one clear gradient diagnostic, use `by_target`.

## ROOT Deployment Branch

Use these names to avoid mixing ROOT inference with dataset testing:

### ROOT -> signal prediction `.npy`

```bash
python3 predictSignalRootToNpy.py \
  --input input.root \
  --model train_run/model.pt \
  --output signal_prediction.npy \
  --class-schema three_class
```

### signal prediction `.npy` -> new ROOT

```bash
python3 writeSignalPredictionToRoot.py \
  --input input.root \
  --prediction signal_prediction.npy \
  --output output_with_signal.root \
  --class-schema three_class
```

These are the preferred names.
For signal ROOT deployment, prefer the two explicit names above.

For batch deployment from a sample json, use:

```bash
python3 runSignalRootInferenceFromSampleJson.py \
  --sample-json path/to/sample.json \
  --model four_class_model/model.pt \
  --output-dir signal_inference_outputs \
  --prefix Signal
```

This creates a namespaced output tree:

- `signal_inference_outputs/_predictions/<json_stem>/...`
- `signal_inference_outputs/data/...`
- `signal_inference_outputs/mc/...`

The wrapper preserves the `data/...` or `mc/...` directory suffix found in the sample json paths, so different periods and channels stay separated in the final scored ROOT outputs. The namespace is only used for the intermediate prediction cache under `_predictions/`. You can override that namespace with `--namespace`.

## BJA Reminder

Preferred BJA entry names:

- `predictBJARootToNpy.py`
- `writeBJAPredictionToRoot.py`

Internal BJA implementation files kept:

- `RNNBJAExtraction_staged.py`
- `ttbb_BJA_model.py`
