import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# --------------------------------------------------------------------------------------------------


def main():
    decay1000kb = np.loadtxt(
        "Bonev_ES_observed_KR_chr06_50-100Mb_res1000kb/optimized_data/decay_cost.txt")
    decay500kb = np.loadtxt(
        "Bonev_ES_observed_KR_chr06_50-100Mb_res500kb/optimized_data/decay_cost.txt")
    decay250kb = np.loadtxt(
        "Bonev_ES_observed_KR_chr06_50-100Mb_res250kb/optimized_data/decay_cost.txt")
    decay100kb = np.loadtxt(
        "Bonev_ES_observed_KR_chr06_50-100Mb_res100kb/optimized_data/decay_cost.txt")
    # ----------------------------------------------------------------------------------------------
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 24
    # ----------------------------------------------------------------------------------------------
    FILE_FIG = "cost_decay.svg"
    plt.figure(figsize=(8, 6))
    plt.xlim(0, 7e5)
    plt.ylim(0, 220)
    plt.xlabel(r"$\mathrm{\mathbf{Iteration}}$")
    plt.ylabel(r"$\mathrm{\mathbf{Cost}}$")

    plt.plot(decay100kb[:, 0], decay100kb[:, 1], linewidth=3,
             label="500x500 (bin size: 100 kb, $r$ = 0.972)", color=cm.viridis(0.0))
    plt.plot(decay250kb[:, 0], decay250kb[:, 1], linewidth=3,
             label="200x200 (bin size: 250 kb, $r$ = 0.984)", color=cm.viridis(0.333))
    plt.plot(decay500kb[:, 0], decay500kb[:, 1], linewidth=3,
             label="100x100 (bin size: 500 kb, $r$ = 0.990)", color=cm.viridis(0.666))
    plt.plot(decay1000kb[:, 0], decay1000kb[:, 1], linewidth=3,
             label="50x50 (bin size: 1 Mb, $r$ = 0.995)", color=cm.viridis(1.0))

    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(FILE_FIG)
    plt.close()
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
