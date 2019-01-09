# PHi-C

PHi-C consists of **Python3** codes for deciphering Hi-C data into polymer dynamics simulation.
The input is a contact matrix data generated from a _hic_ file through [Juicer](https://github.com/aidenlab/juicer).
PHi-C assumes that a genomic region of interest is modeled by a polymer network model including repulsive and attractive interactions between monomers.
Instead of finding optimized 3D conformations, we can obtain optimal parameters of the polymer network model.
Through the optimization procedure, an optimized contact matrix is reconstructed.
Finally, we can carry out polymer dynamics simulations of the polymer network model equipped with the optimal interactions.

![overview](https://github.com/soyashinkai/PHi-C/blob/master/images/overview.png)

### Citation

If you use PHi-C, please cite:

Soya Shinkai, Masaki Nakagawa, Takeshi Sugawara, Yuichi Togashi, Hiroshi Ochiai,
Ryuichiro Nakato, Yuichi Taniguchi, and Shuichi Onami.
**"PHi-C: deciphering Hi-C data into polymer dynamics."**

## Requirements on Python libraries

PHi-C codes require the following Python libraries:
- os
- sys
- numpy
- matplotlib
- scipy
- numba


-----------------
## Usage

### 1. Conversion of a sparce matrix format into a dense contact matrix

Here, _NAME.txt_ as an ipunt is in sparse matrix format produced from [“dump” command of Juicebox](https://github.com/aidenlab/juicer/wiki/Data-Extraction).

    python3 1_conversion.py NAME.txt START END RES

A directory named _NAME_ is made,
and the output _contact_matrix.txt_ is generated in the directory _NAME_.

### 2. Normalization of the dense contact matrix

    python3 2_normalization.py NAME RES OFFSET


### 3. Optimization

    python3 3_optimization.py NAME SAMPLE ALPHA1 ALPHA2 STEP1 STEP2 ITERATION INIT_K_BACKBONE

### 4. Validattion of the optimized contact matrix data 

    python3 4_validation.py NAME RES SAMPLE PLT_MIN_LOG_C PLT_MAX_K_BACKBONE PLT_MAX_K PLT_K_DIS_BINS PLT_MAX_K_DIS

### 5. 4D simulation of the optimal polymer model

    python3 5_4d_simulation.py K-file FRAME

### 6. Normalize the dense contact matrix

    python3 6_conformation.py K-file SAMPLE

