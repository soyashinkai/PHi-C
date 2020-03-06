import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys
import numba
# --------------------------------------------------------------------------------------------------
argv = sys.argv
argc = len(argv)
if (argc != 10):
    print("Usage: python " + argv[0] +
          " DIR RES SAMPLE PLT_MAX_LOG_C PLT_MIN_LOG_C PLT_MAX_K_BACKBONE PLT_MAX_K PLT_K_DIS_BINS PLT_MAX_K_DIS")
    exit()
DIR = argv[1]
RES = int(argv[2])
SAMPLE = int(argv[3])
PLT_MAX_LOG_C = float(argv[4])
PLT_MIN_LOG_C = float(argv[5])
PLT_MAX_K_BACKBONE = float(argv[6])
PLT_MAX_K = float(argv[7])
PLT_K_DIS_BINS = int(argv[8])
PLT_MAX_K_DIS = int(argv[9])
# --------------------------------------------------------------------------------------------------


def Read_Normalized_C():
    FILE_READ = DIR + "/normalized_contact_matrix.txt"
    C = np.loadtxt(FILE_READ)
    N = C.shape[0]
    return C, N
# --------------------------------------------------------------------------------------------------


@numba.jit
def Convert_K_into_C(K, N):
    # K to Laplacian matrix
    d = np.sum(K, axis=0)
    D = np.diag(d)
    L = D - K
    # Eigenvalues and eigenvectors
    lam, Q = np.linalg.eigh(L)
    inv_lam = 1 / lam   # inverse of the eigenvalues
    inv_lam[0] = 0
    inv_Lam = np.diag(inv_lam)
    # L to M
    M = np.dot(Q, np.dot(inv_Lam, Q.T))
    # M to Σ^2
    M_diag = np.diag(np.diag(M))
    Ones = np.ones((N, N))
    A = np.dot(M_diag, Ones)
    Sigma2 = (A + A.T - 2 * M) / 3
    # Σ^2 to C
    C = (1 + Sigma2)**(-1.5)
    return C
# --------------------------------------------------------------------------------------------------


@numba.jit
def Calc_Diff_Cost(log10_C_reconstructed, log10_C_normalized):
    Diff = log10_C_reconstructed - log10_C_normalized
    Cost = np.sqrt(np.trace(np.dot(Diff.T, Diff)))
    return Diff, Cost
# --------------------------------------------------------------------------------------------------


def Calc_Correlation(A, B, N):
    list_X = []
    list_Y = []
    # ----------------------------------------------------------------------------------------------
    for i in range(N):
        for j in range(i, N):
            list_X.append(A[i, j])
            list_Y.append(B[i, j])
    # ----------------------------------------------------------------------------------------------
    X = np.array(list_X)
    Y = np.array(list_Y)
    # ----------------------------------------------------------------------------------------------
    r, p = scipy.stats.pearsonr(X, Y)
    return r, X, Y
# --------------------------------------------------------------------------------------------------


def Calc_Contact_Probability(C):
    N = C.shape[0]
    P = np.zeros((N, 2))
    for n in range(0, N):
        P[n, 0] = RES * n
        cnt = 0
        for m in range(0, N - n):
            P[n, 1] += C[m, m + n]
            cnt += 1
        P[n, 1] /= cnt
    return P
# --------------------------------------------------------------------------------------------------


def main():
    # ----------------------------------------------------------------------------------------------
    DIR_OPTIMIZATION = DIR + "/optimized_data"
    FILE_COST = DIR + "/cost_correlation.txt"
    # ----------------------------------------------------------------------------------------------
    C_normalized, N = Read_Normalized_C()
    log10_C_normalized = np.log10(C_normalized)
    # ----------------------------------------------------------------------------------------------
    stream = open(FILE_COST, "w")
    print("sample\tCost\tCorrelation for logCij", file=stream)
    for sample in range(SAMPLE):
        # ------------------------------------------------------------------------------------------
        # READ matrix K
        FILE_READ = DIR_OPTIMIZATION + "/{0:02d}_K.txt".format(sample)
        K = np.loadtxt(FILE_READ)
        # ------------------------------------------------------------------------------------------
        # CALC cost and correlation
        C_optimized = Convert_K_into_C(K, N)
        log10_C_optimized = np.log10(C_optimized)
        Diff, Cost = Calc_Diff_Cost(log10_C_optimized, log10_C_normalized)
        r, Optimized_log, Normalized_log = Calc_Correlation(
            log10_C_optimized, log10_C_normalized, N)
        print("%d\t%f\t%f" % (sample, Cost / N, r), file=stream)
        # ------------------------------------------------------------------------------------------
        P_normalized = Calc_Contact_Probability(C_normalized)
        P_optimized = Calc_Contact_Probability(C_optimized)
        # ------------------------------------------------------------------------------------------
        C = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                if i > j:
                    C[i, j] = C_optimized[i, j]
                else:
                    C[i, j] = C_normalized[i, j]
        # ------------------------------------------------------------------------------------------
        list_K_wo_BACKBONE = []
        for i in range(0, N - 2):
            for j in range(i + 2, N):
                list_K_wo_BACKBONE.append(K[i, j])
        K_wo_BACKBONE = np.array(list_K_wo_BACKBONE)
        # ------------------------------------------------------------------------------------------
        # OUTPUT K[i, i+1]
        FILE_OUT = DIR_OPTIMIZATION + \
            "/{0:02d}_k_polymer_backbone.txt".format(sample)
        fp = open(FILE_OUT, "w")
        for i in range(N - 1):
            print("%d\t%f" % (i, K[i, i + 1]), file=fp)
        fp.close()
        # ------------------------------------------------------------------------------------------
        plt.style.use("default")
        plt.rcParams["font.family"] = "Arial"
        plt.rcParams["font.size"] = 24
        # ------------------------------------------------------------------------------------------
        FILE_OUT = DIR_OPTIMIZATION + "/{0:02d}_Correlation.svg".format(sample)
        plt.figure(figsize=(5, 5))
        plt.axes().set_aspect("equal")
        x = np.linspace(PLT_MIN_LOG_C, PLT_MAX_LOG_C)
        plt.plot(x, x, linestyle="dashed", color="gray", linewidth=3)
        plt.scatter(Optimized_log, Normalized_log, color="blue", alpha=0.5)
        plt.xlim(PLT_MIN_LOG_C, PLT_MAX_LOG_C + 0.1)
        plt.ylim(PLT_MIN_LOG_C, PLT_MAX_LOG_C + 0.1)
        plt.xticks([PLT_MIN_LOG_C, PLT_MIN_LOG_C / 2, PLT_MAX_LOG_C])
        plt.yticks([PLT_MIN_LOG_C, PLT_MIN_LOG_C / 2, PLT_MAX_LOG_C])
        plt.tight_layout()
        plt.savefig(FILE_OUT)
        plt.close()
        # ------------------------------------------------------------------------------------------
        FILE_READ = DIR_OPTIMIZATION + \
            "/{0:02d}_k_polymer_backbone.txt".format(sample)
        k_profile = np.loadtxt(FILE_READ)
        FILE_OUT = DIR_OPTIMIZATION + \
            "/{0:02d}_k_polymer_backbone.svg".format(sample)
        plt.figure(figsize=(10, 2))
        plt.bar(k_profile[:, 0], k_profile[:, 1],
                width=1.0,
                color="#FF0000",
                linewidth=0)
        plt.ylim(0, PLT_MAX_K_BACKBONE)
        plt.tight_layout()
        plt.savefig(FILE_OUT)
        plt.close()
        # ------------------------------------------------------------------------------------------
        FILE_OUT = DIR_OPTIMIZATION + "/{0:02d}_C.svg".format(sample)
        plt.figure(figsize=(10, 13))
        plt.imshow(C, cmap="magma_r", clim=(0, 1))
        plt.colorbar(ticks=[0, 1], orientation="horizontal", shrink=0.6)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(FILE_OUT)
        plt.close()
        # ------------------------------------------------------------------------------------------
        FILE_OUT = DIR_OPTIMIZATION + "/{0:02d}_C_log.svg".format(sample)
        plt.figure(figsize=(10, 13))
        plt.imshow(np.log10(C), cmap="inferno_r",
                   clim=(PLT_MIN_LOG_C, PLT_MAX_LOG_C))
        plt.colorbar(ticks=[PLT_MIN_LOG_C, PLT_MAX_LOG_C],
                     orientation="horizontal", shrink=0.6)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(FILE_OUT)
        plt.close()
        # ------------------------------------------------------------------------------------------
        FILE_OUT = DIR_OPTIMIZATION + "/{0:02d}_K.svg".format(sample)
        plt.figure(figsize=(10, 13))
        plt.imshow(K, cmap="bwr", clim=(-PLT_MAX_K, PLT_MAX_K))
        plt.colorbar(ticks=[-PLT_MAX_K, 0, PLT_MAX_K],
                     orientation="horizontal", shrink=0.6)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(FILE_OUT)
        plt.close()
        # ------------------------------------------------------------------------------------------
        FILE_OUT = DIR_OPTIMIZATION + \
            "/{0:02d}_contact_probabilities.svg".format(sample)
        plt.figure(figsize=(5, 5))
        plt.ylim(10**PLT_MIN_LOG_C, 10**PLT_MAX_LOG_C)
        plt.xscale("log")
        plt.yscale("log")
        plt.plot(P_normalized[1:, 0], P_normalized[1:, 1],
                 label="Hi-C", linewidth=2)
        plt.plot(P_optimized[1:, 0], P_optimized[1:, 1],
                 label="Optimized", linewidth=2)
        plt.legend(handlelength=1, loc="upper right", fontsize=20)
        plt.tight_layout()
        plt.savefig(FILE_OUT)
        plt.close()
        # ------------------------------------------------------------------------------------------
        FILE_OUT = DIR_OPTIMIZATION + \
            "/{0:02d}_K_distribution.svg".format(sample)
        plt.figure(figsize=(10, 10))
        plt.hist(K_wo_BACKBONE, bins=PLT_K_DIS_BINS,
                 range=(-PLT_MAX_K, PLT_MAX_K), density=True)
        plt.xlim(-PLT_MAX_K, PLT_MAX_K)
        plt.ylim(0, PLT_MAX_K_DIS)
        plt.xticks([-PLT_MAX_K, -PLT_MAX_K / 2, 0, PLT_MAX_K / 2, PLT_MAX_K])
        plt.tight_layout()
        plt.savefig(FILE_OUT)
        plt.close()
    stream.close()
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
