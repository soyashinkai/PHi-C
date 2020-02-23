#!/bin/bash

# Juicertools and URL
JUICERTOOLS="juicer_tools_1.13.02.jar"
URL="http://hicfiles.s3.amazonaws.com/external/bonev/ES_mapq30.hic"

# --------------------------------------------------------------------------------------------------
# Name
NAME="observed_KR_ES_chr06_50-100Mb_res1000kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
CHR=6
RES=1000000
START=50000000
END=99000000
# Run juicertools
java -jar $JUICERTOOLS dump observed KR $URL $CHR:$START:$END $CHR:$START:$END BP $RES $HiCFILE
# --------------------------------------------------------------------------------------------------
# Name
NAME="observed_KR_ES_chr06_50-100Mb_res500kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
CHR=6
RES=500000
START=50000000
END=99500000
# Run juicertools
java -jar $JUICERTOOLS dump observed KR $URL $CHR:$START:$END $CHR:$START:$END BP $RES $HiCFILE
# --------------------------------------------------------------------------------------------------
# Name
NAME="observed_KR_ES_chr06_50-100Mb_res250kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
CHR=6
RES=250000
START=50000000
END=99750000
# Run juicertools
java -jar $JUICERTOOLS dump observed KR $URL $CHR:$START:$END $CHR:$START:$END BP $RES $HiCFILE
# --------------------------------------------------------------------------------------------------
# Name
NAME="observed_KR_ES_chr06_50-100Mb_res100kb"
HiCFILE=$NAME".txt"
# Parameters for conversion and normalization
CHR=6
RES=100000
START=50000000
END=99900000
# Run juicertools
java -jar $JUICERTOOLS dump observed KR $URL $CHR:$START:$END $CHR:$START:$END BP $RES $HiCFILE
# --------------------------------------------------------------------------------------------------
