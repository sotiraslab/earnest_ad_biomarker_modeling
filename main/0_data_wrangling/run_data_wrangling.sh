#!/bin/bash

echo ''
echo 'WRANGLING FOR ADNI'
echo '------------------'

SCRIPT="adni_build_table.R"
Rscript -e "source('${SCRIPT}', echo=T)"

echo ''
echo 'WRANGLING FOR OASIS'
echo '------------------'

SCRIPT="oasis_compute_centiloid.R"
Rscript -e "source('${SCRIPT}', echo=T)"

SCRIPT="oasis_build_table.R"
Rscript -e "source('${SCRIPT}', echo=T)"

echo ''
echo 'HARMONIZATION'
echo '------------------'

SCRIPT="harmonization.R"
Rscript -e "source('${SCRIPT}', echo=T)"

echo ''
echo 'FEATURE ENGINEERING'
echo '------------------'

SCRIPT="feature_engineering.R"
Rscript -e "source('${SCRIPT}', echo=T)"
