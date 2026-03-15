import numpy as np
import scipy.sparse as sp
import scipy.linalg as la

try:
    import paoli_s
    import paoli_is
except ImportError as e:
    raise ImportError(f"C extension module import failed: {e}") from e

def dec2quat_vector(dec_num, n):
    """Convert a decimal number to a quaternary vector"""
    quat = []
    for _ in range(n):
        quat.append(dec_num % 4)
        dec_num = dec_num // 4
    quat = quat + [0] * (n - len(quat))
    return np.array(quat[::-1], dtype=int)

def sppaolidec_w1(h, A, type_, is_symmetric, iscorrect=None):

    # ===================== 1. Required Parameter Validation =====================
    nargin = len(locals())
    if nargin < 4:
        raise ValueError("Error: At least 4 parameters are required: h, A, type, is_symmetric!")

    # ===================== 2. Check input matrix A: MUST be a full matrix (not sparse) =====================
    if sp.issparse(A):
        raise ValueError("Error: Input matrix A must be a FULL matrix, sparse matrix is not allowed!")

    # ===================== 3. h Parameter Validation =====================
    if isinstance(h, (np.integer, np.floating)):
        h = h.item()

    if isinstance(h, str):
        if h != 'off':
            raise ValueError("Error: The 1st parameter h only supports 'off' as a string!")
    elif isinstance(h, (int, float)):
        if not np.isreal(h) or h <= 0 or np.mod(h, 1) != 0:
            raise ValueError("Error: The 1st parameter h must be a positive integer (e.g., 1,5,10) or 'off'!")
        h = int(h)
    else:
        raise TypeError("Error: The 1st parameter h only supports positive integers or the string 'off'!")

    # ===================== 4. type Parameter Validation =====================
    valid_type = {'qua', 'dec', 'off'}
    if type_ not in valid_type:
        raise ValueError("Error: The 3rd parameter (type) must be 'qua', 'dec' or 'off'!")

    # ===================== 5. Symmetry Parameter Validation =====================
    if is_symmetric == "symmetric":
        is_symmetric = True
    elif is_symmetric == "asymmetric":
        is_symmetric = False
    else:
        raise ValueError("Error: is_symmetric must be 'symmetric' or 'asymmetric'!")

    # ===================== 6. Pre-matrix Preprocessing =====================
    if sp.issparse(A):
        A = A.toarray()
    A = np.array(A, dtype=np.float64, copy=False)

    # ===================== 7. Type Conversion =====================
    type_num = 0
    if type_ == 'dec':
        type_num = 1
    elif type_ == 'off':
        type_num = 2

    # ===================== 8. Call C Extension Module =====================
    if is_symmetric:
        a, D = paoli_s.sppaolidec_w1(h, A, type_num)
    else:
        a, D = paoli_is.sppaolidec_w1(h, A, type_num)

    # ===================== 9. Original validation logic (unchanged, no additional parameter checks) =====================
    if type_ == 'off':
        pass
    elif iscorrect == 'open':
        PL = [
            np.eye(2, dtype=np.complex128),
            np.array([[1, 0], [0, -1]], dtype=np.complex128),
            np.array([[0, 1], [1, 0]], dtype=np.complex128),
            np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
        ]

        N = A.shape[0]
        n = int(np.log2(N))
        Ki = np.zeros((N, N), dtype=np.complex128)
        Nd = len(a)

        for i in range(Nd):
            if type_ == 'qua':
                q = D[:, i] + 1
            elif type_ == 'dec':
                q = dec2quat_vector(int(D[i]), n) + 1

            Gi = PL[q[0] - 1]
            for k in range(1, n):
                Gi = np.kron(Gi, PL[q[k] - 1])

            Ki += Gi * a[i]

        er = la.norm(A - Ki.real, 1)
        if er <= 1e-12:
            print(f"Pauli decomposition is correct, error (1-norm): {er:.16f}")
        else:
            print(f"Pauli decomposition is incorrect, error (1-norm): {er:.16f}")

    return a, D