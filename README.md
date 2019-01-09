# PHi-C

PHi-C consists of **Python3** codes for deciphering Hi-C data into polymer dynamics simulation.
The input is a contact matrix data generated from a _hic_ file through [Juicer](https://github.com/aidenlab/juicer).
PHi-C assumes that a genomic region of interest is modeled by a polymer network model including repulsive and attractive interactions between monomers.
Instead of finding optimized 3D conformations, we can obtain optimal parameters of the polymer network model.
Through the optimization procedure, an optimized contact matrix is reconstructed.
Finally, we can carry out polymer dynamics simulations of the polymer network model equipped with the optimal interactions.

![overview](/images/overview.png)

### Citation

If you use PHi-C, please cite:

Soya Shinkai, Masaki Nakagawa, Takeshi Sugawara, Yuichi Togashi, Hiroshi Ochiai,
Ryuichiro Nakato, Yuichi Taniguchi, and Shuichi Onami.
**"PHi-C: deciphering Hi-C data into polymer dynamics."**

## Requirements

PHi-C codes require the following Python libraries:
- os
- sys
- numpy
- matplotlib
- scipy
- numba

To visualize the simulated polymer dynamics and conformations, [VMD](https://www.ks.uiuc.edu/Research/vmd/) is needed.


## Quick Start
Move to the directory [_Tutorial_](/Tutorial):

    cd Tutorial

Then, run the following scripts:

    ./demo_run_1-2.sh
    ./demo_run_3-4.sh
    ./demo_run_5-6.sh

It will take a few seconds, 20 minutes or less, and a few minutes, respectively.

You will obtain the same output data and figures in the directory [_Demos_](/Demos).
We used Hi-C data for mouse embryo stem cells (chr8: 42,100-44,525 kb; bin size: 25 kb) with KR normalization by [Bonev et al.](https://doi.org/10.1016/j.cell.2017.09.043).

### Visualization of polymer conformations

To visualize the simulated polymer dynamics, run VMD, firstly read _polymer_N97.psf_, and then read _dynamics_00_K.xyz_ on VMD.

To visualize the sampled polymer conformations, run VMD, firstly read _polymer_N97.psf_, and then read _conformations_00_K.xyz_ on VMD.


-------------------------
## Usage
PHi-C consists of the following six Python codes:
- 1_conversion.py
- 2_normalization.py
- 3_optimization.py
- 4_validation.py
- 5_4d_simulation.py
- 6_conformations.py

### 1. Conversion of a sparce matrix format into a dense contact matrix
    python3 1_conversion.py NAME.txt START END RES

Here, _NAME.txt_ as an ipunt is in sparse matrix format produced from [“dump” command of Juicebox](https://github.com/aidenlab/juicer/wiki/Data-Extraction).

A directory named _NAME_ is made, and the output _contact_matrix.txt_ is generated in the directory _NAME_.

- START:
- END:
- RES:

### 2. Normalization of the dense contact matrix

    python3 2_normalization.py NAME RES OFFSET

- RES: 
- OFFSET: 

### 3. Optimization

    python3 3_optimization.py NAME SAMPLE ALPHA1 ALPHA2 STEP1 STEP2 ITERATION INIT_K_BACKBONE

- SAMPLE:
- ALPHA1:
- ALPHA2:
- STEP1:
- STEP2:
- ITERATION:
- INIT_K_BACKBONE

### 4. Validation of the optimized contact matrix data 

    python3 4_validation.py NAME RES SAMPLE PLT_MIN_LOG_C PLT_MAX_K_BACKBONE PLT_MAX_K PLT_K_DIS_BINS PLT_MAX_K_DIS

- PLT_MIN_LOG_C
- PLT_MAX_K_BACKBONE
- PLT_MAX_K
- PLT_K_DIS_BINS
- PLT_MAX_K_DIS

### 5. 4D simulation of the optimal polymer dynamics

    python3 5_4d_simulation.py KFILE FRAME

- KFILE:
- FRAME:

### 6. Sampling the optimal polymer conformations

    python3 6_conformation.py KFILE SAMPLE

- SAMPLE:
