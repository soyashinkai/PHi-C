import numpy as np
import sys
import os

# --------------------------------------------------------------------------------------------------
# Set command line arguments
argv = sys.argv
argc = len(argv)
if (argc != 3):
    print("Usage: python " + argv[0] + " K-file SAMPLE")
    exit()
FILE_READ = argv[1]
SAMPLE = int(argv[2])
FILE_NAME = os.path.basename(FILE_READ)
NAME, EXT = os.path.splitext(FILE_NAME)
FILE_OUT = "conformations_" + NAME + ".xyz"
# --------------------------------------------------------------------------------------------------
SEED = 3932
np.random.seed(SEED)
# --------------------------------------------------------------------------------------------------


def Read_K():
    K = np.loadtxt(FILE_READ)
    N = K.shape[0]
    return K, N
# --------------------------------------------------------------------------------------------------


def Transform_K_into_L(K):
    d = np.sum(K, axis=0)
    D = np.diag(d)
    L = D - K
    return L
# --------------------------------------------------------------------------------------------------


def Equilibrium_Conformation_of_Normal_Coordinates(lam, P):
    Xx = np.zeros(P)
    Xy = np.zeros(P)
    Xz = np.zeros(P)
    for p in range(1, P):
        Xx[p] = np.sqrt(1 / 3 / lam[p]) * np.random.randn()
        Xy[p] = np.sqrt(1 / 3 / lam[p]) * np.random.randn()
        Xz[p] = np.sqrt(1 / 3 / lam[p]) * np.random.randn()
    return Xx, Xy, Xz
# --------------------------------------------------------------------------------------------------


def Convert_X_to_R(Xx, Xy, Xz, Q):
    Rx = np.dot(Q, Xx)
    Ry = np.dot(Q, Xy)
    Rz = np.dot(Q, Xz)
    return Rx, Ry, Rz
# --------------------------------------------------------------------------------------------------


def Write_Psfdata(N):
    FILE_PSF = "polymer_N" + str(N) + ".psf"
    fp = open(FILE_PSF, "w")
    print("PSF\n\n       1 !NTITLE\n REMARKS %s\n" % NAME, file=fp)
    print(" %7d !NATOM" % N, file=fp)
    for n in range(N):
        print(" %7d A    %04d GLY  CA   CT1    0.070000       12.0110           0" % (
            n + 1, n + 1), file=fp)
    print("\n %7d !NBOND: bonds" % (N - 1), file=fp)
    j = 0
    for i in range(N):
        if i % N != 0:
            print(" %7d %7d" % (i, i + 1), end="", file=fp)
            j += 1
            if j % 4 == 0:
                print("\n", end="", file=fp)
    print("\n", end="", file=fp)
    fp.close()
# --------------------------------------------------------------------------------------------------


def main():
    K, N = Read_K()
    P = N
    Write_Psfdata(N)
    # ----------------------------------------------------------------------------------------------
    L = Transform_K_into_L(K)
    lam, Q = np.linalg.eigh(L)
    # ----------------------------------------------------------------------------------------------
    fp = open(FILE_OUT, "w")
    for sample in range(SAMPLE):
        Xx, Xy, Xz = Equilibrium_Conformation_of_Normal_Coordinates(lam, P)
        Rx, Ry, Rz = Convert_X_to_R(Xx, Xy, Xz, Q)
        # ------------------------------------------------------------------------------------------
        print("%d" % N, file=fp)
        print("sample = %d" % sample, file=fp)
        for n in range(N):
            print("CA\t%f\t%f\t%f" % (Rx[n], Ry[n], Rz[n]), file=fp)
        # ------------------------------------------------------------------------------------------
    fp.close()
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
