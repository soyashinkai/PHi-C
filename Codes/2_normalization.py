import numpy as np
import matplotlib.pyplot as plt
import sys
# --------------------------------------------------------------------------------------------------
argv = sys.argv
argc = len(argv)
if (argc != 6):
    print("Usage: python " + argv[0] +
          " DIR RES OFFSET PLT_MAX_LOG_C PLT_MIN_LOG_C")
    exit()
DIR = argv[1]
RES = int(argv[2])
OFFSET = float(argv[3])
PLT_MAX_LOG_C = float(argv[4])
PLT_MIN_LOG_C = float(argv[5])
# --------------------------------------------------------------------------------------------------
FILE_READ = DIR + "/contact_matrix.txt"
FILE_OUT_MATRIX = DIR + "/normalized_contact_matrix.txt"
FILE_OUT_PROBABILITY = DIR + "/normalized_contact_probability.txt"
FILE_LOG = DIR + "/interpolation.log"
# --------------------------------------------------------------------------------------------------


def Calc_Contact_Probability(C):
    N = C.shape[0]
    P = np.zeros((N, 2))
    fp = open(FILE_LOG, "w")
    for n in range(0, N):
        P[n, 0] = RES * n
        cnt = 0
        for m in range(0, N - n):
            if C[m, m + n] > 0:
                P[n, 1] += C[m, m + n]
                cnt += 1
            else:
                print("C[%d, %d]\tis a 0 value element." % (m, m + n), file=fp)
        if cnt > 0:
            P[n, 1] /= cnt
        else:
            P[n, 1] = OFFSET
            print("Interpolated: P[n=%d] = %f" % (n, OFFSET), file=fp)
    fp.close()
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
                C_interpolated[i, j] = P[j - i, 1]
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
                C_normalized[i, j] = C_interpolated[i, j] / P[0, 1]
                if C_normalized[i, j] > 1:
                    C_normalized[i, j] = P[j - i, 1] / P[0, 1]
            C_normalized[j, i] = C_normalized[i, j]
    return C_normalized
# --------------------------------------------------------------------------------------------------


def Calc_Contact_Normalized_Probability(C_normalized):
    N = C_normalized.shape[0]
    P_normalized = np.zeros((N, 2))
    for n in range(0, N):
        P_normalized[n, 0] = RES * n
        cnt = 0
        for m in range(0, N - n):
            if C_normalized[m, m + n] > 0:
                P_normalized[n, 1] += C_normalized[m, m + n]
                cnt += 1
        if cnt > 0:
            P_normalized[n, 1] /= cnt
    return P_normalized
# --------------------------------------------------------------------------------------------------


def Plot_Figs(C_normalized, P, P_normalized):
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 24
    # ----------------------------------------------------------------------------------------------
    FILE_FIG = DIR + "/normalized_Cij_log.svg"
    plt.figure(figsize=(8, 8))
    plt.imshow(np.log10(C_normalized), cmap="inferno_r",
               clim=(PLT_MIN_LOG_C, PLT_MAX_LOG_C))
    plt.colorbar(ticks=[PLT_MIN_LOG_C, PLT_MAX_LOG_C],
                 shrink=0.6, orientation="horizontal")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(FILE_FIG)
    plt.close()
    # ----------------------------------------------------------------------------------------------
    FILE_FIG = DIR + "/normalized_Cij.svg"
    plt.figure(figsize=(8, 8))
    plt.imshow(C_normalized, cmap="magma_r", clim=(0, 1))
    plt.colorbar(ticks=[0, 1], shrink=0.6, orientation="horizontal")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(FILE_FIG)
    plt.close()
    # ----------------------------------------------------------------------------------------------
    FILE_FIG = DIR + "/contact_probabilities.svg"
    plt.figure(figsize=(6, 6))
    plt.style.use("ggplot")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    plt.plot(P[1:, 0], P[1:, 1], label="JUICER")
    plt.plot(P_normalized[1:, 0], P_normalized[1:, 1],
             label="Normalized s.t. Cii = 1")
    plt.legend(handlelength=1, loc="lower left")
    plt.tight_layout()
    plt.savefig(FILE_FIG)
    plt.close()
# --------------------------------------------------------------------------------------------------


def main():
    C = np.loadtxt(FILE_READ)
    P = Calc_Contact_Probability(C)
    C_interpolated = Interpolate_Contact_Matrix(C, P)
    C_normalized = Calc_Normalized_Contact_Matrix(C_interpolated, P)
    np.savetxt(FILE_OUT_MATRIX, C_normalized, fmt="%e")
    # ----------------------------------------------------------------------------------------------
    P_normalized = Calc_Contact_Normalized_Probability(C_normalized)
    np.savetxt(FILE_OUT_PROBABILITY, P_normalized, fmt="%d\t%e")
    # ----------------------------------------------------------------------------------------------
    Plot_Figs(C_normalized, P, P_normalized)
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
