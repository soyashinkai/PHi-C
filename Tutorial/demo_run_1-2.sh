#!/bin/bash

# Name
NAME="Bonev_ES_observed_KR_chr8_42100-44500kb_res25kb"
HiCFILE=$NAME".txt"

# Parameters for conversion and normalization
START=42100000
END=44500000
RES=25000
OFFSET=2
PLT_MIN_LOG_C=-3

#---------------------------------------------------------------------------------------------------
# Run python codes
python3 1_conversion.py $HiCFILE $START $END $RES

python3 2_normalization.py $NAME $RES $OFFSET $PLT_MIN_LOG_C
