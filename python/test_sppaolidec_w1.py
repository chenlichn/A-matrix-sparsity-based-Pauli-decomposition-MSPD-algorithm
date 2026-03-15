#Python test environment: Python 3.12 (based on the PyCharmMiscProject)
import h5py
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import time
from sppaolidec_w1 import sppaolidec_w1

def load_mat73(mat_path, var_name):
    """Read MATLAB 7.3+ .mat files (HDF5 format)"""
    with h5py.File(mat_path, 'r') as f:
        data = f[var_name][:]
        data = data.T
    return data

def matrix_bandwidth(mat, symmetry_desc):
    # Convert to sparse matrix if the input is a dense matrix
    if not sp.issparse(mat):
        mat_sparse = sp.csr_matrix(mat)
    else:
        mat_sparse = mat.copy()

    n = mat_sparse.shape[0]
    lower_bw = 0
    upper_bw = 0

    # Get row and column indices of non-zero elements
    rows, cols = mat_sparse.nonzero()
    for r, c in zip(rows, cols):
        if r > c:
            lower_bw = max(lower_bw, r - c)
        elif r < c:
            upper_bw = max(upper_bw, c - r)

    # Determine final bandwidth based on matrix symmetry
    if symmetry_desc == 'symmetric':
        final_bw = upper_bw
    else:
        final_bw = max(upper_bw, lower_bw)

    final_bw = int(final_bw)
    # Ensure the minimum bandwidth is 1
    final_bw = max(final_bw, 1)

    return int(lower_bw), int(upper_bw), final_bw

if __name__ == "__main__":
    # 1. Load the matrix
    load_path = r'Path to store K1 matrix\K1.mat' # K1 matrix is stored in v7.3 format and must be a full matrix
    K1 = load_mat73(load_path, 'K1')

    # 2. New: Check matrix symmetry in main program (migrated from sppaoli function)
    # Principle: Use Frobenius norm to check difference between matrix and its transpose; <1e-13 is considered symmetric
    symmetry_norm = la.norm(K1 - K1.T, 'fro')  # Calculate Frobenius norm of the difference between the matrix and its transpose
    if symmetry_norm < 1e-13:
        symmetry_desc = 'symmetric'  # Condition true: symmetric matrix
    else:
        symmetry_desc = 'asymmetric'  # Condition false: asymmetric matrix
    print(f'Matrix symmetry judgment result: {symmetry_desc} (Frobenius norm difference: {symmetry_norm:.6e})')

    # 3. Calculate matrix bandwidth
    Lb, Ub, bandwidth_python = matrix_bandwidth(K1, symmetry_desc)


    # 4. Initialize parameters
    total_time = 0.0
    loop_count = 10
    h = bandwidth_python
    type_flag = 'dev'
    is_symmetric = symmetry_desc
    iscorrect_flag = 'off'


    # 5. Run decomposition in a loop
    for i in range(loop_count):
        start_time = time.perf_counter()

        a, DD = sppaolidec_w1(
            h=h,
            A=K1,
            type_=type_flag,
            is_symmetric=is_symmetric,
            iscorrect=iscorrect_flag
        )
        elapsed_time = time.perf_counter() - start_time
        print(f'Runtime for iteration {i + 1}: {elapsed_time:.6f} seconds')
        total_time += elapsed_time

    # 6. Output statistical results
    average_time = total_time / loop_count
    print('------------------------')
    print(f'Average runtime: {average_time:.6f} seconds')