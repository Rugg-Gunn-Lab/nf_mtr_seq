# Dev branch
# Currently this will only work for paired mode

This branch uses STARsolo instead of STAR + custom counting script. This means that we do not double count clustered features, instead only keeping unique ones, this resolves some strange observations in UMAPs.

Added the `--soloTE` param, when provided this will ALSO run soloTE as well as standard mapping, so you will get 2 output matrices, one for TE features and one for non-TE (standard output) features. If `--soloTE` is provided, the user must also ensure soloTE has been cloned locally, the TE bed file has been generated (via SoloTE_RepeatMasker_to_BED.py helper function) and then both the local soloTE dir and TE bed file must be provided via `--soloTE_dir` and `--soloTE_bed` respectively. 