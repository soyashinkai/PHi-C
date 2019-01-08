import numpy as np
import sys
import os

#---------------------------------------------------------------------------------------------------
# Set command line arguments
argv = sys.argv
argc = len(argv)
if (argc != 3):
    print("Usage: python " + argv[0] + " K-file FRAME")
    exit()
FILE_READ = argv[1]
FRAME = int(argv[2])
FILE_NAME = os.path.basename(FILE_READ)
NAME, EXT = os.path.splitext(FILE_NAME)
FILE_OUT = "dynamics_" + NAME + ".xyz"
#---------------------------------------------------------------------------------------------------
# Parameters
EPS = 0.0001  # = kB T Δt / (γ σ^2): a non-dimensional parameter
INTERVAL = 1000
#---------------------------------------------------------------------------------------------------
# Physical constants
NOISE = np.sqrt(2 * EPS)
F_Coefficient = 3 * EPS
#---------------------------------------------------------------------------------------------------
SEED = 2980
np.random.seed(SEED)
#---------------------------------------------------------------------------------------------------


def Read_K():
    K = np.loadtxt(FILE_READ)
    N = K.shape[0]
    return K, N
#---------------------------------------------------------------------------------------------------


def Transform_K_into_L(K):
    d = np.sum(K, axis=0)
    D = np.diag(d)
    L = D - K
    return L
#---------------------------------------------------------------------------------------------------


def Equilibrium_Conformation_of_Normal_Coordinates(lam, P):
    Xx = np.zeros(P)
    Xy = np.zeros(P)
    Xz = np.zeros(P)
    for p in range(1, P):
        Xx[p] = np.sqrt(1 / 3 / lam[p]) * np.random.randn()
        Xy[p] = np.sqrt(1 / 3 / lam[p]) * np.random.randn()
        Xz[p] = np.sqrt(1 / 3 / lam[p]) * np.random.randn()
    return Xx, Xy, Xz
#---------------------------------------------------------------------------------------------------


def Convert_X_to_R(Xx, Xy, Xz, Q):
    Rx = np.dot(Q, Xx)
    Ry = np.dot(Q, Xy)
    Rz = np.dot(Q, Xz)
    return Rx, Ry, Rz
#---------------------------------------------------------------------------------------------------


def Integrate_Polymer_Network(x, y, z, L, N):
    noise_x = NOISE * np.random.randn(N)
    noise_y = NOISE * np.random.randn(N)
    noise_z = NOISE * np.random.randn(N)

    force_x = - F_Coefficient * np.dot(L, x)
    force_y = - F_Coefficient * np.dot(L, y)
    force_z = - F_Coefficient * np.dot(L, z)

    x_dt = x + force_x + noise_x
    y_dt = y + force_y + noise_y
    z_dt = z + force_z + noise_z

    force_x = - F_Coefficient * np.dot(L, x_dt)
    force_y = - F_Coefficient * np.dot(L, y_dt)
    force_z = - F_Coefficient * np.dot(L, z_dt)

    x_2dt = x_dt + force_x + noise_x
    y_2dt = y_dt + force_y + noise_y
    z_2dt = z_dt + force_z + noise_z

    X = 0.5 * (x + x_2dt)
    Y = 0.5 * (y + y_2dt)
    Z = 0.5 * (z + z_2dt)
    return X, Y, Z
#---------------------------------------------------------------------------------------------------


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
#---------------------------------------------------------------------------------------------------


def main():
    K, N = Read_K()
    P = N
    Write_Psfdata(N)
    #-----------------------------------------------------------------------------------------------
    L = Transform_K_into_L(K)
    lam, Q = np.linalg.eigh(L)
    Xx, Xy, Xz = Equilibrium_Conformation_of_Normal_Coordinates(lam, P)
    Rx, Ry, Rz = Convert_X_to_R(Xx, Xy, Xz, Q)
    #-----------------------------------------------------------------------------------------------
    fp = open(FILE_OUT, "w")
    for frame in range(FRAME + 1):
        #---------------------------------------------------------------------------------------
        print("%d" % N, file=fp)
        print("frame = %d" % frame, file=fp)
        for n in range(N):
            print("CA\t%f\t%f\t%f" % (Rx[n], Ry[n], Rz[n]), file=fp)
        #---------------------------------------------------------------------------------------
        for step in range(INTERVAL):
            Rx, Ry, Rz = Integrate_Polymer_Network(Rx, Ry, Rz, L, N)
    fp.close()
#---------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
