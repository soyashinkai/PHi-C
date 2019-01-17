#!/bin/bash

# Name
KNUM="00"
KFILE="Bonev_ES_observed_KR_chr8_42100-44500kb_res25kb/optimized_data/"$KNUM"_K.txt"

# Parameters
FRAME=1000
SAMPLE=100

#---------------------------------------------------------------------------------------------------
# Run python codes
python3 5_4d_simulation.py $KFILE $FRAME

python3 6_conformation.py $KFILE $SAMPLE
