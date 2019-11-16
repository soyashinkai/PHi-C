import numpy as np
import scipy.stats
import sys
import numba
# --------------------------------------------------------------------------------------------------
argv = sys.argv
argc = len(argv)
if (argc != 2):
    print("Usage: python " + argv[0] + " RES")
    exit()
RES = argv[1]
DIR = "Bonev_ES_observed_KR_chr06_50-100Mb_res{0:s}kb".format(RES)
FILE_OUT = "run-time_res{0:s}kb.txt".format(RES)
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


def Calc_Correlation(A, B, N):
    NN = int(N * (N + 1) / 2)
    X = np.zeros(NN)
    Y = np.zeros(NN)
    n = 0
    for i in range(N):
        for j in range(i, N):
            X[n] = A[i, j]
            Y[n] = B[i, j]
            n += 1
    r, p = scipy.stats.pearsonr(X, Y)
    return r
# --------------------------------------------------------------------------------------------------


def main():
    # ----------------------------------------------------------------------------------------------
    DIR_OPT = DIR + "/optimized_data"
    fp = open(FILE_OUT, "a")
    # ----------------------------------------------------------------------------------------------
    C_normalized, N = Read_Normalized_C()
    log10_C_normalized = np.log10(C_normalized)
    # ----------------------------------------------------------------------------------------------
    # READ matrix K
    FILE_READ = DIR_OPT + "/iterated-sample00_K.txt"
    K = np.loadtxt(FILE_READ)
    # ----------------------------------------------------------------------------------------------
    # CALC cost and correlation
    C_optimized = Convert_K_into_C(K, N)
    log10_C_optimized = np.log10(C_optimized)
    r = Calc_Correlation(log10_C_optimized, log10_C_normalized, N)
    print("Correlation: %f" % r, file=fp)
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
