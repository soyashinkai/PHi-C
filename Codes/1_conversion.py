import numpy as np
import os
import sys
# --------------------------------------------------------------------------------------------------
argv = sys.argv
argc = len(argv)
if (argc != 3):
    print("Usage: python " + argv[0] + " Hi-C-FILE RES")
    exit()
FILE_READ = argv[1]
RES = int(argv[2])
DIR, EXT = os.path.splitext(FILE_READ)
os.makedirs(DIR, exist_ok=True)
FILE_OUT = DIR + "/contact_matrix.txt"
# --------------------------------------------------------------------------------------------------


def Read_Data():
    data = np.loadtxt(FILE_READ, delimiter="\t")
    data[np.isnan(data)] = 0
    return data
# --------------------------------------------------------------------------------------------------


def Convert_Data_into_Matrix(data):
    Len = data.shape[0]
    START = data[0, 1]
    END = data[Len - 1, 1]
    N = int((END - START) / RES) + 1
    matrix = np.zeros((N, N))
    M = data.shape[0]
    for m in range(M):
        i = int((data[m, 0] - START) / RES)
        j = int((data[m, 1] - START) / RES)
        matrix[i, j] = data[m, 2]
        matrix[j, i] = data[m, 2]
    return matrix
# --------------------------------------------------------------------------------------------------


def main():
    data = Read_Data()
    matrix = Convert_Data_into_Matrix(data)
    np.savetxt(FILE_OUT, matrix, fmt="%f")
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
