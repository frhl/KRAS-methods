#!/usr/bin/env bash
#
#$ -N hypergeom
#$ -wd /well/lindgren/flassen/projects/kras/MassSpectrometry
#$ -o logs/hypergeom.log
#$ -e logs/hypergeom.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe


set_up_rpy() {
  module load Anaconda3/2020.07
  module load java/1.8.0_latest
  source "/apps/eb/skylake/software/Anaconda3/2020.07/etc/profile.d/conda.sh"
  conda activate /well/lindgren/users/mmq446/conda/skylake/envs/rpy
}

readonly rscript="workflows/KRAS/211001_volcano_hypergeometric_enrichment.R"

set_up_rpy
Rscript "${rscript}"

