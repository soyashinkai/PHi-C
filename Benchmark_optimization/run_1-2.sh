#!/bin/bash

# --------------------------------------------------------------------------------------------------
# Name
NAME="Bonev_ES_observed_KR_chr06_50-100Mb_res100kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
RES=100000
START=50000000
END=99900000
OFFSET=2
PLT_MIN_LOG_C=-3
# Run python codes
python 1_conversion_1.0.1.py $HiCFILE $START $END $RES
python 2_normalization_1.0.1.py $NAME $RES $OFFSET $PLT_MIN_LOG_C
# --------------------------------------------------------------------------------------------------
# Name
NAME="Bonev_ES_observed_KR_chr06_50-100Mb_res250kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
RES=250000
START=50000000
END=99750000
OFFSET=2
PLT_MIN_LOG_C=-3
# Run python codes
python 1_conversion_1.0.1.py $HiCFILE $START $END $RES
python 2_normalization_1.0.1.py $NAME $RES $OFFSET $PLT_MIN_LOG_C
# --------------------------------------------------------------------------------------------------
# Name
NAME="Bonev_ES_observed_KR_chr06_50-100Mb_res500kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
RES=500000
START=50000000
END=99500000
OFFSET=2
PLT_MIN_LOG_C=-3
# Run python codes
python 1_conversion_1.0.1.py $HiCFILE $START $END $RES
python 2_normalization_1.0.1.py $NAME $RES $OFFSET $PLT_MIN_LOG_C
# --------------------------------------------------------------------------------------------------
# Name
NAME="Bonev_ES_observed_KR_chr06_50-100Mb_res1000kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
RES=1000000
START=50000000
END=99000000
OFFSET=2
PLT_MIN_LOG_C=-3
# Run python codes
python 1_conversion_1.0.1.py $HiCFILE $START $END $RES
python 2_normalization_1.0.1.py $NAME $RES $OFFSET $PLT_MIN_LOG_C
# --------------------------------------------------------------------------------------------------
