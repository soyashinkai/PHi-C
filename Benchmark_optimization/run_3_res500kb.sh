#!/bin/bash

# Name
NAME="Bonev_ES_observed_KR_chr06_50-100Mb_res500kb"

# Hyperparameters for optimization
ALPHA1=0.002
ALPHA2=0.0001
STEP1=2000
STEP2=5000
ITERATION=100
INIT_K_BACKBONE=0.5

# Number of optimization sample
SAMPLE=1

#---------------------------------------------------------------------------------------------------
# Run python codes
(time python3 3_optimization_1.0.1.py $NAME $SAMPLE $ALPHA1 $ALPHA2 $STEP1 $STEP2 $ITERATION $INIT_K_BACKBONE) 2> run-time_res500kb.txt
