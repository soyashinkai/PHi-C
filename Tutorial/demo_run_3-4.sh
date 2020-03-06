#!/bin/bash

# Name
NAME=Bonev_ES_observed_KR_chr8_42100-44500kb_res25kb
RES=25000

# Hyperparameters for optimization
ALPHA1=0.002
ALPHA2=0.0001
STEP1=2000
STEP2=5000
ITERATION=100
INIT_K_BACKBONE=0.3

# Number of optimization sample
SAMPLE=1

# Parameters for plot
PLT_MAX_LOG_C=0
PLT_MIN_LOG_C=-3
PLT_MAX_K_BACKBONE=0.5
PLT_MAX_K=0.010
PLT_K_DIS_BINS=100
PLT_MAX_K_DIS=1500

#---------------------------------------------------------------------------------------------------
# Run python codes
python3 3_optimization.py ${NAME} ${SAMPLE} ${ALPHA1} ${ALPHA2} ${STEP1} ${STEP2} ${ITERATION} ${INIT_K_BACKBONE}

python3 4_validation.py ${NAME} ${RES} ${SAMPLE} ${PLT_MAX_LOG_C} ${PLT_MIN_LOG_C} ${PLT_MAX_K_BACKBONE} ${PLT_MAX_K} ${PLT_K_DIS_BINS} ${PLT_MAX_K_DIS}
