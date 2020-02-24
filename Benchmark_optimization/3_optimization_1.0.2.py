import numpy as np
import os
import sys
import numba
# --------------------------------------------------------------------------------------------------
argv = sys.argv
argc = len(argv)
if (argc != 9):
    print("Usage: python " +
          argv[0] + " DIR SAMPLE ALPHA1 ALPHA2 STEP1 STEP2 ITERATION INIT_K_BACKBONE")
    exit()
DIR = argv[1]
SAMPLE = int(argv[2])
# Hyperparameters
ALPHA1 = float(argv[3])  # for (i, i+1)
ALPHA2 = float(argv[4])  # for (i, j)
STEP1 = int(argv[5])  # for (i, i+1)
STEP2 = int(argv[6])  # for (i, j)
ITERATION = int(argv[7])
INIT_K_BACKBONE = float(argv[8])
# --------------------------------------------------------------------------------------------------
SEED = 3939
np.random.seed(SEED)
# --------------------------------------------------------------------------------------------------


def Read_Normalized_C():
    FILE_READ = DIR + "/normalized_contact_matrix.txt"
    C = np.loadtxt(FILE_READ)
    N = C.shape[0]
    return C, N
# --------------------------------------------------------------------------------------------------


def Init_K(N):
    K = np.zeros((N, N))
    for i in range(N - 1):
        j = i + 1
        K[i, j] = K[j, i] = INIT_K_BACKBONE
    return K
# --------------------------------------------------------------------------------------------------


@numba.jit
def Check_PSD_of_L(K):
    # K to Laplacian matrix
    d = np.sum(K, axis=0)
    D = np.diag(d)
    L = D - K
    # Eigenvalues and eigenvectors
    lam, Q = np.linalg.eigh(L)
    if ~np.all(lam[1:] > 0):
        print("Error: The Laplacian matrix does not the positive-semidefiniteness!")
        exit()
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


def Random_Int_Pair_Adjacent(N):
    i = np.random.randint(N - 1)
    j = i + 1
    return i, j
# --------------------------------------------------------------------------------------------------


def Random_Int_Pair_Different(N):
    while True:
        i = np.random.randint(N)
        j = np.random.randint(N)
        if i != j:
            break
    return i, j
# --------------------------------------------------------------------------------------------------


def main():
    DIR_OPT = DIR + "/optimized_data"
    os.makedirs(DIR_OPT, exist_ok=True)
    # ----------------------------------------------------------------------------------------------
    C_normalized, N = Read_Normalized_C()
    log10_C_normalized = np.log10(C_normalized)
    K = Init_K(N)
    C_reconstructed = Convert_K_into_C(K, N)
    log10_C_reconstructed = np.log10(C_reconstructed)
    Diff, Cost = Calc_Diff_Cost(log10_C_reconstructed, log10_C_normalized)
    # ----------------------------------------------------------------------------------------------
    FILE_LOG = DIR_OPT + "/optimization.log"
    fp_log = open(FILE_LOG, "w")
    FILE_DECAY_COST = DIR_OPT + "/decay_cost.txt"
    fp_decay = open(FILE_DECAY_COST, "w")
    cnt = 0
    # ----------------------------------------------------------------------------------------------
    for sample in range(SAMPLE):
        for iteration in range(ITERATION):
            # --------------------------------------------------------------------------------------
            print("sample\titeration\tstep1\tCost\tUpdate K[i,j]", file=fp_log)
            for step1 in range(STEP1):
                cnt += 1
                tmp_K = K.copy()

                i, j = Random_Int_Pair_Adjacent(N)

                delta = np.random.uniform(0, ALPHA1)
                if Diff[i, j] < 0:
                    tmp_K[i, j] += delta
                    tmp_K[j, i] += delta
                else:
                    tmp_K[i, j] -= delta
                    tmp_K[j, i] -= delta

                C_reconstructed = Convert_K_into_C(tmp_K, N)
                log10_C_reconstructed = np.log10(C_reconstructed)
                tmp_Diff, tmp_Cost = Calc_Diff_Cost(
                    log10_C_reconstructed, log10_C_normalized)

                if tmp_Cost < Cost:
                    K = tmp_K.copy()
                    Diff = tmp_Diff.copy()
                    Cost = tmp_Cost
                    print("%d\t%d\t%04d\t%f\tK[%d,%d]" %
                          (sample, iteration, step1, Cost, i, j), file=fp_log)
                    print("%d\t%f" % (cnt, Cost), file=fp_decay)
                    Check_PSD_of_L(K)
            # --------------------------------------------------------------------------------------
            print("sample\titeration\tstep2\tCost\tUpdate K[i,j]", file=fp_log)
            for step2 in range(STEP2):
                cnt += 1
                tmp_K = K.copy()

                i, j = Random_Int_Pair_Different(N)

                delta = np.random.uniform(0, ALPHA2)
                if Diff[i, j] < 0:
                    tmp_K[i, j] += delta
                    tmp_K[j, i] += delta
                else:
                    tmp_K[i, j] -= delta
                    tmp_K[j, i] -= delta

                C_reconstructed = Convert_K_into_C(tmp_K, N)
                log10_C_reconstructed = np.log10(C_reconstructed)
                tmp_Diff, tmp_Cost = Calc_Diff_Cost(
                    log10_C_reconstructed, log10_C_normalized)

                if tmp_Cost < Cost:
                    K = tmp_K.copy()
                    Diff = tmp_Diff.copy()
                    Cost = tmp_Cost
                    print("%d\t%d\t%04d\t%f\tK[%d,%d]" %
                          (sample, iteration, step2, Cost, i, j), file=fp_log)
                    print("%d\t%f" % (cnt, Cost), file=fp_decay)
                    Check_PSD_of_L(K)
        # ------------------------------------------------------------------------------------------
        FILE_OUT = DIR_OPT + "/iterated-sample{0:02d}_K.txt".format(sample)
        np.savetxt(FILE_OUT, K, fmt="%e")
    # ----------------------------------------------------------------------------------------------
    fp_log.close()
    fp_decay.close()
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
