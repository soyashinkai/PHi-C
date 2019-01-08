#!/bin/bash

KNUM="00"

# Parameters
FRAME=1000
SAMPLE=100

KFILE="Bonev_ES_observed_KR_chr8_42100-44500kb_res25kb/optimized_data/"$KNUM"_K.txt"

# Run python codes
python 5_4d_simulation.py $KFILE $FRAME

python 6_conformation.py $KFILE $SAMPLE
