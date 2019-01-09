# PHi-C

PHi-C consists of Python3 codes for deciphering Hi-C data into polymer dynamics simulation.
The input is a contact matrix data generated from a _hic_ file through [Juicer](https://github.com/aidenlab/juicer).

## Citation

If you use PHi-C, please cite:

Soya Shinkai, Masaki Nakagawa, Takeshi Sugawara, Yuichi Togashi, Hiroshi Ochiai,
Ryuichiro Nakato, Yuichi Taniguchi, and Shuichi Onami.
**"PHi-C: deciphering Hi-C data into polymer dynamics."**

# Requirements

# Usage
Here, _FILENAME.txt_ as an ipunt is in sparse matrix format produced from [“dump” command of Juicebox](https://github.com/aidenlab/juicer/wiki/Data-Extraction).

    python3 1_conversion.py FILENAME.txt START END RESOLUTION

A directory named _FILENAME_ is made,
and the output _contact_matrix.txt_ is generated in the directory _FILENAME_.

    python3
