# PHi-C

PHi-C consists of **Python3** codes for deciphering Hi-C data into polymer dynamics simulations.
The input is a contact matrix data generated from a _hic_ file through [Juicer](https://github.com/aidenlab/juicer).
PHi-C assumes that a genomic region of interest at an appropriate resolution can be modeled using a polymer network model including attractive and repulsive interactions between monomers.
Instead of finding optimized 3D conformations, PHi-C's optimization procedure provides optimal interaction parameters of the polymer network model.
We can then reconstruct an optimized contact matrix.
Finally, we can carry out polymer dynamics simulations of the polymer network model equipped with the optimal interaction parameters.

![overview](/images/overview.png)


### Citation

If you use PHi-C, please cite:

Soya Shinkai, Masaki Nakagawa, Takeshi Sugawara, Yuichi Togashi, Hiroshi Ochiai,
Ryuichiro Nakato, Yuichi Taniguchi, and Shuichi Onami.
**"PHi-C: deciphering Hi-C data into polymer dynamics."**
*bioRxiv* [https://doi.org/10.1101/574962](https://doi.org/10.1101/574962) (2019).

## Requirements

PHi-C codes require the following Python libraries:

-   os
-   sys
-   numpy
-   matplotlib
-   scipy
-   numba

To visualize the simulated polymer dynamics and conformations, [VMD](https://www.ks.uiuc.edu/Research/vmd/) is needed.

## Quick Start

Move to the directory [_Tutorial_](/Tutorial):

    cd Tutorial

Then, run the following scripts:

    ./demo_run_1-2.sh
    ./demo_run_3-4.sh
    ./demo_run_5-6.sh

It will take a few seconds, 20 minutes or less, and a few minutes, respectively.

[](You will obtain the same output data and figures in the directory [_Demos_](/Demos).)
We used Hi-C data for mouse embryo stem cells (chr8: 42,100-44,525 kb; bin size: 25 kb) with KR normalization by [Bonev et al.](https://doi.org/10.1016/j.cell.2017.09.043).

### Visualization of polymer conformations

To visualize the simulated polymer dynamics, run VMD, firstly read _polymer_N97.psf_, and then read _dynamics_00_K.xyz_ on VMD.

To visualize the sampled polymer conformations, run VMD, firstly read _polymer_N97.psf_, and then read _conformations_00_K.xyz_ on VMD.

* * *

## Usage

PHi-C consists of the following six Python codes:

-   1_conversion.py
-   2_normalization.py
-   3_optimization.py
-   4_validation.py
-   5_4d_simulation.py
-   6_conformations.py

### 1. Conversion of a sparce matrix format into a dense contact matrix

Here, _NAME.txt_ as an example ipunt is in sparse matrix format produced from [“dump” command of Juicebox](https://github.com/aidenlab/juicer/wiki/Data-Extraction).

    python3 1_conversion.py NAME.txt START END RES

The command converts to dense matrix format data, _contact_matrix.txt_, at the newly made directory _NAME_:  
_./NAME/contact_matrix.txt_.

The other three arguments of the command represent the followings:

-   START: the start genomic coordinate of the input Hi-C data,
-   END: the end genomic coordinate of the input Hi-C data,
-   RES: the bin size or resolution of the input Hi-C data.

### 2. Normalization of the dense contact matrix

    python3 2_normalization.py NAME RES OFFSET PLT_MIN_LOG_C

The command normalizes the Hi-C matrix data, _./NAME/contact_matrix.txt_, so that the diagonal elements satisfy _C<sub>ii</sub>_ = 1 as probability, with an interpolation if needed.  
The output six files are the followings:  
_./NAME/normalized_contact_matrix.txt_  
_./NAME/normalized_contact_probability.txt_  
_./NAME/normalized_Cij.svg_  
_./NAME/normalized_Cij_log.svg_  
_./NAME/contact_probability.svg_  
_./NAME/interpolation.log_

The other two arguments of the command represent the followings:

-   RES: the bin size or resolution of the input Hi-C data,
-   OFFSET: the offset value for ND contact probability _P(s)_ if needed,
-   PLT_MIN_LOG_C: the lower limit to plot _./NAME/normalized_Cij_log.svg_.

### 3. Optimization

    python3 3_optimization.py NAME SAMPLE ALPHA1 ALPHA2 STEP1 STEP2 ITERATION INIT_K_BACKBONE

The command carries out the PHi-C optimization.  
The log data of the optimization are stored as _./NAME/optimized_data/optimization.log_.  
The optimized interaction matrix data of the polymer network model with _SAMPLE-INDEX_ are output as
_./NAME/optimized_data/{SAMPLE-INDEX}\_K.txt_.

The other seven arguments of the command represent the followings:

-   SAMPLE: the number of samples to obtain optimized output,
-   ALPHA1: the learning rate of the optimization for _k<sub>i, i+1</sub>_,
-   ALPHA2: the learning rate of the optimization for _k<sub>i, j</sub>_,
-   STEP1: the number of the optimization steps for _k<sub>i, i+1</sub>_,
-   STEP2: the number of the optimization steps for _k<sub>i, j</sub>_,
-   ITERATION: the number to iterate (STEP1 + STEP2) optimization steps,
-   INIT_K_BACKBONE: initial values of \_k<sub>i, i+1</sub>_.

### 4. Validation of the optimized contact matrix data

    python3 4_validation.py NAME RES SAMPLE PLT_MIN_LOG_C PLT_MAX_K_BACKBONE PLT_MAX_K PLT_K_DIS_BINS PLT_MAX_K_DIS

The output eight files are the followings:  
_./NAME/cost_correlation.txt_  
_./NAME/optimized_data/{SAMPLE-INDEX}\_C.svg_  
_./NAME/optimized_data/{SAMPLE-INDEX}\_C_log.svg_  
_./NAME/optimized_data/{SAMPLE-INDEX}\_contact_probabilities.svg_  
_./NAME/optimized_data/{SAMPLE-INDEX}\_Correlation.svg_  
_./NAME/optimized_data/{SAMPLE-INDEX}\_K.svg_  
_./NAME/optimized_data/{SAMPLE-INDEX}\_K_distribution.svg_  
_./NAME/optimized_data/{SAMPLE-INDEX}\_k_polymer_backbone.svg_

The other five arguments of the command represent the followings:

-   PLT_MIN_LOG_C: the lower limit to plot _{SAMPLE-INDEX}\_C_log.svg_, _{SAMPLE-INDEX}\_contact_probabilities.svg_ and _{SAMPLE-INDEX}\_Correlation.svg_
-   PLT_MAX_K_BACKBONE: the upper limit to plot _{SAMPLE-INDEX}\_k_polymer_backbone.svg_,
-   PLT_MAX_K: the upper and lower limit to plot _{SAMPLE-INDEX}\_K.svg\_,
-   PLT_K_DIS_BINS: the number of the bins of _{SAMPLE-INDEX}\_K_distribution.svg_,
-   PLT_MAX_K_DIS: the upper limit to plot _{SAMPLE-INDEX}\_K_distribution.svg_.

### 5. 4D simulation of the optimal polymer dynamics

    python3 5_4d_simulation.py KFILE FRAME

The output two files are the followings:  
_./polymer_N{NUMBER-OF-BEADS}.psf_  
_./dynamics\_{INPUT-KFILE}.xyz_

The other two arguments of the command represent the followings:

-   KFILE: the optimized interaction matrix data of the polymer network as input,
-   FRAME: the number of frames to visualize simulated polymer dynamics.

### 6. Sampling the optimal polymer conformations

    python3 6_conformation.py KFILE SAMPLE

The output two files are the followings:  
_./polymer_N{NUMBER-OF-BEADS}.psf_  
_./conformations\_{INPUT-KFILE}.xyz_

The other argument of the command represents the followings:

-   KFILE: the optimized interaction matrix data of the polymer network as input,
-   SAMPLE: the number of polymer conformations to caluclate.
