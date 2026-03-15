- - - # Code Description for Matrix Sparsity-Based Pauli Decomposition (Python Version)

      **Version**: v1.0

      **Authors**: Feng Wu, Chen Li, Xuanlong Wu, Yansong Guo, Li Zhu, Yuxiang Yang

      **Date**: March 15, 2026

      ------

      ## 1. Introduction

      ### 1.1 Document Purpose

      This document is designed to help users quickly understand, use, and maintain the Python implementation of the **Matrix-Sparsity-Based Pauli Decomposition (MSPD)** algorithm. It specifies the function, operation procedures, and dependent environment of the code, to reduce the cost of deployment, usage, and maintenance.

      The MSPD algorithm implemented in this toolkit is based on the method proposed in the literature: Feng Wu, Chen Li, Xuanlong Wu, Yansong Guo, Xu Guo, *A matrix sparsity-based Pauli decomposition algorithm*, 2026.

      ### 1.2 Project Background

      Efficient numerical simulation is the core of modern engineering mechanics, yet the storage and computational demands of discrete matrices for large-scale mechanical problems have far exceeded the limits of classical computers. Quantum computing offers a promising path to break through this computational bottleneck, and **Pauli decomposition** is a critical step to convert mechanical matrices into linear combinations of unitary matrices that can be operated by quantum circuits.

      Existing Pauli decomposition methods are mostly designed for dense matrices, which fail to fully leverage the banded sparsity characteristic of mechanical matrices, leading to limited computational efficiency. The **Matrix-Sparsity-Based Pauli Decomposition (MSPD)** method implemented in this code, by exploring the sparsity of Pauli strings and the banded structure of mechanical matrices combined with the fast Hadamard transform (FHT), significantly reduces the computational complexity. It provides an efficient preprocessing tool for the integration of quantum computing and computational mechanics.

      ### 1.3 Scope of Application

      - Applicable to the Pauli decomposition of square matrices with a size of **power of 2** (N=2^n, where n is the number of qubits);
      - Supports **symmetric/asymmetric banded sparse matrices** (input must be in dense format), and also delivers excellent efficiency for dense matrices;
      - Outputs the coefficients of Pauli decomposition and the corresponding indices of Pauli strings, which can be directly used in the subsequent workflow of quantum computing;
      - The input matrix must be a full (dense) matrix; sparse matrices need to be converted to dense format before input.

      ## 2. System/Environment Requirements

      ### 2.1 Hardware Requirements

      - Memory: 8GB or above is recommended (larger memory is required for processing high-dimensional matrices);
      - CPU: Supports common instruction sets (AVX2 instruction set is recommended for optimal performance; BMI/LZCNT support is required for hardware acceleration).

      ### 2.2 Software Requirements

      - **Python Version**: Python 3.12 (compatibility with 3.8+ versions needs to be verified);
      - **Dependent Python Libraries**: `numpy` (numerical computation), `scipy` (sparse matrix and linear algebra operations), `h5py` (reading MATLAB v7.3 .mat files);
      - **C Compiler**:
        - Windows: Microsoft Visual C++ (MSVC) compatible with the Python version (recommended for the compilation script provided), or MinGW64 Compiler;
        - Linux: GCC compiler;
        - macOS: Clang compiler;
      - **Build Tool**: `setuptools` (for compiling C extension modules).

      ## 3. Installation and Configuration

      ### 3.1 Code File List

      Ensure that all the following files are in the same working directory:

      | Filename                       | File Type                 | Core Function Description                                    |
      | ------------------------------ | ------------------------- | ------------------------------------------------------------ |
      | `paoli_is.c`                   | C Source File             | C source code of Pauli decomposition for asymmetric matrices |
      | `paoli_s.c`                    | C Source File             | C source code of Pauli decomposition for symmetric matrices  |
      | `paoli_is.cp313-win_amd64.pyd` | Python Extension Module   | Precompiled asymmetric decomposition module for Windows (Python 3.13, x64) |
      | `paoli_s.cp313-win_amd64.pyd`  | Python Extension Module   | Precompiled symmetric decomposition module for Windows (Python 3.13, x64) |
      | `setup_paoli_is.py`            | Python Compilation Script | Compilation configuration script for the asymmetric matrix decomposition C extension |
      | `setup_paoli_s.py`             | Python Compilation Script | Compilation configuration script for the symmetric matrix decomposition C extension |
      | `sppaolidec_w1.py`             | Python Main Function      | Main interface function for MSPD, including parameter validation and result validation |
      | `test_sppaolidec_w1.py`        | Python Test Script        | Complete test script with matrix loading, symmetry check, timing, and result verification |

      > **Note**: The suffix of the precompiled Python extension module varies with operating system and Python version:
      >
      > - Windows: `.cpXX-win_amd64.pyd` (XX corresponds to the Python version, e.g., 313 for Python 3.13)
      > - Linux: `.cpython-XX-x86_64-linux-gnu.so`
      > - macOS (Intel): `.cpython-XX-darwin.so`
      > - macOS (Apple Silicon): `.cpython-XX-darwin-arm64.so`
      >
      > If there is no precompiled extension module matching your system and Python version, please compile the C extension first according to the steps in Section 3.3.

      ### 3.2 Install Dependent Libraries

      Execute the following command in the terminal/command prompt to install the required Python libraries:

      ```
      pip install numpy scipy h5py setuptools
      ```

      ### 3.3 Compile C Extension Modules

      If there is no matching precompiled extension module, or you need to recompile the C code, execute the following commands in the terminal/command prompt (ensure the C compiler is correctly configured and the working directory is the code directory):

      1. Compile the decomposition module for symmetric matrices:

      ```
      python setup_paoli_s.py build_ext --inplace
      ```

      1. Compile the decomposition module for asymmetric matrices:

      ```
      python setup_paoli_is.py build_ext --inplace
      ```

      After successful compilation, the corresponding Python extension module (`.pyd`/`.so` file) will be generated in the current directory.

      ## 4. Operation Guide

      ### 4.1 Program Startup

      1. Configure the Python environment (ensure the dependent libraries are installed and the C extension modules are compiled successfully);
      2. Use a Python IDE (such as PyCharm, VS Code) or terminal to open the code directory;
      3. Run the test script or call the main interface function in your own code.

      ### 4.2 Core Function Interface

      #### Main Interface Function: `sppaolidec_w1`

      ```
      a, D = sppaolidec_w1(h, A, type_, is_symmetric, iscorrect=None)
      ```

      **Input Parameter Description**:

      | Parameter Name | Type / Value Range                                    | Parameter Description                                        |
      | -------------- | ----------------------------------------------------- | ------------------------------------------------------------ |
      | `h`            | Positive integer / string `'off'`                     | **Required**. Matrix bandwidth: directly used if a positive integer is input; automatically calculated internally if `'off'` is input. |
      | `A`            | Numpy array (dense square matrix, size of power of 2) | Matrix to be decomposed (**only dense numpy array is supported**, both symmetric and asymmetric matrices are acceptable). |
      | `type_`        | `'qua'`/`'dec'`/`'off'`                               | Output type: `'qua'` = quaternary index, `'dec'` = decimal index, `'off'` = no index output. |
      | `is_symmetric` | `'symmetric'`/`'asymmetric'`                          | Matrix symmetry: `'symmetric'` for symmetric matrix, `'asymmetric'` for asymmetric matrix. |
      | `iscorrect`    | `'open'`/`'off'` (Optional)                           | Validation switch: `'open'` = verify the correctness of decomposition, `'off'`/`None` = disable validation (disabled by default). |

      **Output Parameter Description**:

      | Parameter Name | Type                        | Output Description                                           |
      | -------------- | --------------------------- | ------------------------------------------------------------ |
      | `a`            | Numpy array (vector)        | Vector of Pauli decomposition coefficients (arranged in non-zero order); may be complex for asymmetric matrix decomposition. |
      | `D`            | Numpy array (matrix/vector) | Indices of Pauli strings: n×length(a) quaternary matrix when `type_='qua'`; 1×length(a) decimal vector when `type_='dec'`; empty array when `type_='off'`. |

      ### 4.3 Execution Flow

      1. **Input Preparation**: Construct or load the matrix A to be decomposed (ensure the order is a power of 2, and the matrix is a dense numpy array);
      2. **Parameter Setting**: Set `h`, `type_`, and `is_symmetric` according to the matrix characteristics (note that `is_symmetric` must be passed as the specified string);
      3. **Call Decomposition**: Run the `sppaolidec_w1` function, and the code will automatically complete the following process:
         - Strict parameter validity check;
         - Call the corresponding C extension module according to the matrix symmetry;
         - Calculate/use the specified bandwidth;
         - Generate indices of Pauli strings;
         - Perform fast Hadamard transform;
         - Filter non-zero coefficients and return the results;
      4. **Result Validation**: If `iscorrect='open'` is enabled, the code will reconstruct the matrix from the decomposition results and output the 1-norm error of the decomposition.

      ## 5. Examples and Precautions

      ### 5.1 Functional Example: Run the Test Script `test_sppaolidec_w1.py`

      `test_sppaolidec_w1.py` is a complete test script, including MATLAB v7.3 matrix loading, automatic matrix symmetry judgment, bandwidth calculation, loop decomposition timing, and result statistics.

      **Operation Steps**:

      1. Open `test_sppaolidec_w1.py` in your Python IDE or editor;

      2. Modify the `load_path` variable to the actual path of your matrix file (the example uses `K1.mat`, which must be a dense matrix stored in MATLAB v7.3 format);

      3. Run the script directly in the IDE, or execute the following command in the terminal:

         ```
         python test_sppaolidec_w1.py
         ```

      **Expected Output**:

      - Matrix symmetry judgment result (including the Frobenius norm difference between the matrix and its transpose);
      - Runtime of each loop iteration;
      - Average runtime of multiple cycles;
      - If `iscorrect='open'` is enabled, the decomposition error (1-norm) and correctness verification result will be output.

      ### 5.2 Custom Call Examples

      #### Example 1: Decomposition of Symmetric Banded Matrix (Specified Bandwidth)

      ```
      import numpy as np
      import scipy.sparse as sp
      from sppaolidec_w1 import sppaolidec_w1
      
      # Construct a symmetric banded matrix of order 2^13 (bandwidth h=6)
      n_qubit = 13
      N = 2 ** n_qubit
      h = 6
      A = sp.csr_matrix((N, N), dtype=np.float64)
      
      for i in range(N):
          j_start = max(0, i - h)
          j_end = min(N, i + h + 1)
          A[i, j_start:j_end] = np.random.randn(j_end - j_start)
      
      # Convert to symmetric matrix and dense format
      A = (A + A.T) / 2
      A = A.toarray()  # Must be converted to a dense numpy array
      
      # Call decomposition (specify bandwidth h=6, decimal output, symmetric matrix, enable validation)
      a, D = sppaolidec_w1(h=6, A=A, type_='dec', is_symmetric='symmetric', iscorrect='open')
      ```

      #### Example 2: Decomposition of Asymmetric Matrix (Internal Automatic Bandwidth Calculation)

      ```
      import numpy as np
      from sppaolidec_w1 import sppaolidec_w1
      
      # Construct a random asymmetric matrix of order 2^4
      n_qubit = 4
      N = 2 ** n_qubit
      A = np.random.rand(N, N)  # Directly a dense numpy array
      
      # Call decomposition (internal bandwidth calculation, quaternary output, asymmetric matrix, disable validation)
      a, D = sppaolidec_w1(h='off', A=A, type_='qua', is_symmetric='asymmetric', iscorrect='off')
      ```

      ### 5.3 Precautions

      1. **Matrix Order Limitation**: The order of the input matrix A must be a power of 2 (e.g., 2, 4, 8, 16…), otherwise an error will be reported in the C extension module. If the matrix order does not meet the requirement, zero padding or truncation must be performed first. The maximum tested matrix order is up to 2^15.
      2. **Input Matrix Format Requirement**: The input matrix A must be a **dense numpy array**. Direct input of `scipy` sparse matrices is not supported, otherwise a parameter error will be triggered.
      3. **Parameter Format Specification**:
         - `is_symmetric` must be passed as the string `'symmetric'` or `'asymmetric'`, boolean values or other strings are not allowed;
         - The parameter name for output type is `type_` (with a trailing underscore), not `type`, to avoid conflict with Python's built-in keyword.
      4. **C Extension Compilation Issues**:
         - If compilation fails on Windows, ensure that the Microsoft Visual C++ Build Tools are installed and match the Python version;
         - For Linux/macOS, ensure that GCC/Clang is installed, and modify the compilation parameters in the setup script if needed;
         - If the compiled module cannot be imported, check whether the generated `.pyd`/`.so` file is in the current working directory, and whether the Python version matches the compilation environment.
      5. **Time Consumption Reminder for Validation**: When the matrix scale is large (e.g., n_qubit > 12), enabling `iscorrect='open'` will lead to a significant increase in time consumption due to the Kronecker product operation and matrix reconstruction. It is recommended to enable validation only for small matrix testing.
      6. **Bandwidth Parameter Optimization**: If the exact bandwidth of the matrix is known, it is recommended to directly pass the integer value of `h`, which can skip the internal bandwidth calculation step in the C module and further improve the decomposition efficiency.
      7. **MAT File Reading**: The built-in `load_mat73` function only supports MATLAB v7.3 .mat files (HDF5 format). For earlier versions of .mat files, use `scipy.io.loadmat` to read the matrix instead.

      ------

      ## Document Update Log

      - **v1.0** (March 15, 2026): Initial version of the Python code description, covering core function description, installation and compilation guide, usage examples, and precautions.