## Methods and scripts for "The role of germline TE insertions in pediatric cancer predisposition"

### Software used within

-  [xTea](https://github.com/parklab/xTea)
-  [AEon](https://github.com/GenomicRisk/aeon)
-  [goleft/indexcov](https://github.com/brentp/goleft/blob/master/indexcov/README.md)
-  bcftools, python (pyranges)

-  `run_xtea.sh` to generate xTea calls.
-  `postprocess_xTea/postprocess_xTEA.sh` to annotate and merge xTea calls within 35bp
-  `sample_filtering_steps/run_nextflow_meiqc.sh` to get ancestry, depth, and indexcov estimates
-  For ancestry-matched control cohort generation see `matched_control_generation/`


