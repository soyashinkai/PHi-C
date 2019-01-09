#!/bin/bash

# Name
NAME="Bonev_ES_observed_KR_chr8_42100-44500kb_res25kb"
HiCFILE=$NAME".txt"

# Parameters for conversion and normalization
START=42100000
END=44500000
RES=25000
OFFSET=2

#---------------------------------------------------------------------------------------------------
# Run python codes
python 1_conversion.py $HiCFILE $START $END $RES

python 2_normalization.py $NAME $RES $OFFSET
