import numpy as np
import matplotlib.pyplot as plt
import sys
# --------------------------------------------------------------------------------------------------
argv = sys.argv
argc = len(argv)
if (argc != 4):
    print("Usage: python " + argv[0] + " DIR RES OFFSET")
    exit()
DIR = argv[1]
RES = int(argv[2])
OFFSET = float(argv[3])
# --------------------------------------------------------------------------------------------------
FILE_READ = DIR + "/contact_matrix.txt"
FILE_OUT_MATRIX = DIR + "/normalized_contact_matrix.txt"
FILE_OUT_PROBABILITY = DIR + "/normalized_contact_probability.txt"
# --------------------------------------------------------------------------------------------------


def Calc_Contact_Probability(C):
    N = C.shape[0]
    P = np.zeros(N)
    for n in range(0, N):
        count = 0
        for m in range(0, N - n):
            if C[m, m + n] > 0:
                P[n] += C[m, m + n]
                count += 1
            else:
                print("C[%d, %d]\tis a 0 value element." % (m, m + n))
        if count > 0:
            P[n] /= count
        else:
            P[n] = OFFSET
            print("Interpolated: P[n=%d] = %f" % (n, OFFSET))
    return P
# --------------------------------------------------------------------------------------------------


def Interpolate_Contact_Matrix(C, P):
    N = C.shape[0]
    C_interpolated = np.zeros((N, N))
    for i in range(0, N):
        for j in range(i, N):
            if C[i, j] > 0:
                C_interpolated[i, j] = C[i, j]
            else:
                C_interpolated[i, j] = P[j - i]
            C_interpolated[j, i] = C_interpolated[i, j]
    return C_interpolated
# --------------------------------------------------------------------------------------------------


def Calc_Normalized_Contact_Matrix(C_interpolated, P):
    N = C_interpolated.shape[0]
    C_normalized = np.zeros((N, N))
    for i in range(0, N):
        for j in range(i, N):
            if j == i:
                C_normalized[i, j] = 1.0
            else:
                C_normalized[i, j] = C_interpolated[i, j] / P[0]
                if C_normalized[i, j] > 1:
                    C_normalized[i, j] = P[j - i] / P[0]
            C_normalized[j, i] = C_normalized[i, j]
    return C_normalized
# --------------------------------------------------------------------------------------------------


def Plot_Figs(C_normalized, P, P_normalized):
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 24
    # ----------------------------------------------------------------------------------------------
    FILE_FIG = DIR + "/normalized_Cij_log.svg"
    plt.figure(figsize=(8, 8))
    plt.imshow(np.log10(C_normalized), cmap="inferno_r", clim=(-3, 0))
    plt.colorbar(ticks=[-3, 0], shrink=0.6, orientation="horizontal")
    plt.tick_params(labelbottom=0, labelleft=0, color="white")
    plt.tight_layout()
    plt.savefig(FILE_FIG)
    plt.clf()
    # ----------------------------------------------------------------------------------------------
    FILE_FIG = DIR + "/normalized_Cij.svg"
    plt.figure(figsize=(8, 8))
    plt.imshow(C_normalized, cmap="magma_r", clim=(0, 1))
    plt.colorbar(ticks=[0, 1], shrink=0.6, orientation="horizontal")
    plt.tick_params(labelbottom=0, labelleft=0, color="white")
    plt.tight_layout()
    plt.savefig(FILE_FIG)
    plt.clf()
    # ----------------------------------------------------------------------------------------------
    FILE_FIG = DIR + "/contact_probabilities.svg"
    plt.figure(figsize=(6, 6))
    plt.style.use("ggplot")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    N = C_normalized.shape[0]
    s = np.zeros(N)
    for n in range(0, N):
        s[n] = n * RES
    plt.plot(s[1:], P[1:], label="JUICER")
    plt.plot(s[1:], P_normalized[1:], label="Normalized s.t. Cii = 1")
    plt.legend(handlelength=1, loc="lower left")
    plt.tight_layout()
    plt.savefig(FILE_FIG)
    plt.clf()
    return s
# --------------------------------------------------------------------------------------------------


def main():
    C = np.loadtxt(FILE_READ)
    P = Calc_Contact_Probability(C)
    C_interpolated = Interpolate_Contact_Matrix(C, P)
    C_normalized = Calc_Normalized_Contact_Matrix(C_interpolated, P)
    P_normalized = Calc_Contact_Probability(C_normalized)
    np.savetxt(FILE_OUT_MATRIX, C_normalized, fmt="%f")
    s = Plot_Figs(C_normalized, P, P_normalized)
    np.savetxt(FILE_OUT_PROBABILITY, np.c_[
               s, P_normalized], fmt=["%d", "%f"])
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
