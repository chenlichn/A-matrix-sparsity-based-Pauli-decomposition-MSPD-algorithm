- - # Code Description for Matrix-Sparsity-Based Pauli Decomposition

    **Version**: v1.0

    **Authors**: Feng Wu, Chen Li, Xuanlong Wu, Yansong Guo, Li Zhu, Yuxiang Yang

    **Date**: March 15, 2026

    ------

    ## 1. Introduction

    ### 1.1 Document Purpose

    This document aims to help users quickly understand, use, and maintain the code for the **Matrix-Sparsity-Based Pauli Decomposition (MSPD)** algorithm. It clarifies the function, operation method, and dependent environment of the code, to reduce the cost of use and maintenance.

    The MSPD algorithm implemented in this toolkit is based on the method proposed in the following literature: Feng Wu, Chen Li, Xuanlong Wu, Yansong Guo, Xu Guo, *A matrix sparsity-based Pauli decomposition algorithm*, 2026.

    ### 1.2 Project Background

    Efficient numerical simulation is the core of modern engineering mechanics, yet the storage and computational requirements of discrete matrices for large-scale mechanical problems have far exceeded the limits of classical computers. Quantum computing brings hope for breaking through the computational bottleneck, and **Pauli decomposition** is an essential step to convert mechanical matrices into combinations of unitary matrices operable by quantum circuits.

    Existing Pauli decomposition methods are mostly designed for dense matrices, which fail to fully exploit the banded sparsity characteristic of mechanical matrices, resulting in limited computational efficiency. The **Matrix-Sparsity-Based Pauli Decomposition (MSPD)** algorithm implemented in this code, by exploring the sparsity of Pauli strings and the banded structure of mechanical matrices combined with the fast Hadamard transform (FHT), significantly reduces the computational complexity. It provides an efficient preprocessing tool for the integration of quantum computing and computational mechanics.

    ### 1.3 Scope of Application
  
    - Applicable to the Pauli decomposition of square matrices with a size of **power of 2** (N=2^n, where n is the number of qubits);
    - Supports **symmetric/asymmetric banded sparse matrices** (input must be in dense format), and also delivers good efficiency for dense matrices;
    - Outputs the coefficients of Pauli decomposition and the corresponding indices of Pauli strings, which can be directly used in the subsequent workflow of quantum computing;
    - The input matrix must be a full matrix; sparse matrices need to be converted via the `full` function before input.
  
    ## 2. System/Environment Requirements
  
    ### 2.1 Hardware Requirements

    - Memory: 8GB or above is recommended (larger memory is required for processing high-dimensional matrices);
    - CPU: Supports common instruction sets (BMI/LZCNT support is required for hardware-accelerated instructions).
  
    ### 2.2 Software Requirements

    - **MATLAB Version**: MATLAB R2023a (compatibility needs to be verified for other versions);
    - **C Compiler**: MinGW64 Compiler (C) (for compiling MEX functions, configuration via `mex -setup` in MATLAB is required).
  
    ## 3. Installation and Configuration
  
    ### 3.1 Code File List
  
    Ensure that the following files are in the same directory:
  
    | Filename            | File Type       | Core Function Description                                    |
    | ------------------- | --------------- | ------------------------------------------------------------ |
    | `dec2quat_vector.m` | MATLAB Function | Convert decimal numbers to fixed-length quaternary vectors   |
    | `paoli_is.c`        | C Source File   | MEX source code for Pauli decomposition of asymmetric matrices |
    | `paoli_s.c`         | C Source File   | MEX source code for Pauli decomposition of symmetric matrices |
    | `paoli_is.mexw64`   | MEX Executable  | Compiled decomposition function for asymmetric matrices (Windows) |
    | `paoli_s.mexw64`    | MEX Executable  | Compiled decomposition function for symmetric matrices (Windows) |
    | `sppaoli.m`         | MATLAB Function | Main interface function for Pauli decomposition              |
    | `test_example.m`    | MATLAB Script   | Main test file (with examples and timing)                    |
  
    > **Note**: The suffix of MEX files is related to the operating system:
    >
    > 
    >
    > - Windows: `.mexw64`
    > - Linux: `.mexa64`
    > - macOS (Intel): `.mexmaci64`
    > - macOS (Apple Silicon): `.mexmaca64`
    >
    > 
    >
    > If there is no precompiled MEX file for your system in the directory, please compile first according to the steps in Section 3.2.
  
    ### 3.2 Compile MEX Functions
  
    If there is no MEX file corresponding to your system, or you need to recompile, go to the code directory in the MATLAB command line and execute the following commands:
  
    ```
    % Compile the decomposition function for symmetric matrices
    mex paoli_s.c
    
    % Compile the decomposition function for asymmetric matrices
    mex paoli_is.c
    ```
  
    If a "target specific option mismatch" error occurs during compilation (related to hardware instruction sets), you can use the following commands to enable native CPU instruction set support:
  
    ```
    mex CFLAGS="$CFLAGS -march=native" paoli_s.c
    mex CFLAGS="$CFLAGS -march=native" paoli_is.c
    ```

    After successful compilation, MEX files such as `paoli_s.mexw64` and `pauli_is.mexw64` (for Windows system) will be generated in the directory.
  
    ### 3.3 Configure MATLAB Path
  
    Add the code directory to the MATLAB search path:
  
    ```
    % Temporarily add (valid for the current session)
    addpath('Path of the code directory');
    
    % Permanently add (path saving is required)
    savepath;
    ```
  
    ## 4. Operation Guide
  
    ### 4.1 Program Startup
  
    1. Open MATLAB R2023a;
    2. Enter the code directory via the `cd` command, or open the code directory in the "Current Folder" panel of MATLAB.

    ### 4.2 Core Function Interface

    #### Main Interface Function: `sppaoli`
  
    ```
    [a, D] = sppaoli(h, A, type, is_symmetric, iscorrect)
    ```
  
    **Input Parameter Description**:
  
    | Parameter Name | Type / Value Range                       | Parameter Description                                        |
    | -------------- | ---------------------------------------- | ------------------------------------------------------------ |
    | `h`            | Positive integer / string `'off'`        | **Required**. Matrix bandwidth: directly used if a positive integer is input; automatically calculated internally if `'off'` is input. |
    | `A`            | Dense square matrix (size of power of 2) | Matrix to be decomposed (**only dense format is supported**, both symmetric and asymmetric matrices are acceptable). |
    | `type`         | `'qua'`/`'dec'`/`'off'`                  | Output type: `'qua'` = quaternary index, `'dec'` = decimal index, `'off'` = no output. |
    | `is_symmetric` | `'symmetric'`/`'asymmetric'`             | Matrix symmetry: `'symmetric'` for symmetric matrix, `'asymmetric'` for asymmetric matrix. |
    | `iscorrect`    | `'open'`/`'off'` (Optional)              | Validation switch: `'open'` = verify the correctness of decomposition, `'off'` = disable validation (disabled by default). |

    **Output Parameter Description**:

    | Parameter Name | Type            | Output Description                                           |
    | -------------- | --------------- | ------------------------------------------------------------ |
    | `a`            | Vector          | Vector of Pauli decomposition coefficients (arranged in non-zero order); may be complex for asymmetric matrix decomposition. |
    | `D`            | Matrix / Vector | Indices of Pauli strings: n×length(a) quaternary matrix when `type='qua'`; 1×length(a) decimal vector when `type='dec'`. |
  
    ### 4.3 Execution Flow
  
    1. **Input Preparation**: Construct or load the matrix A to be decomposed (ensure the order is a power of 2 and the matrix is in dense format);
    2. **Parameter Setting**: Set `h`, `type`, and `is_symmetric` according to the matrix characteristics (note that `is_symmetric` must be passed as the specified string);
    3. **Call Decomposition**: Run the `sppaoli` function, and the code will automatically complete the following process:
       - Calculate/use the specified bandwidth;
       - Generate indices of Pauli strings;
       - Perform fast Hadamard transform;
       - Filter non-zero coefficients and output results;
    4. **Log Monitoring**: The command line window will output timing information for each step (such as bandwidth calculation time, core computing time, etc.).

    ## 5. Examples and Precautions

    ### 5.1 Functional Example: Run the Test Main File `test_example`
  
    `test_example` is a complete test script, including matrix loading, automatic symmetry judgment, bandwidth calculation, loop decomposition timing, and result verification functions.
  
    **Operation Steps**:
  
    1. Open `test_example.m` in MATLAB;
  
    2. Modify `load_path` in the script to the actual path of the matrix file (in the example, it is `K1.mat`, which needs to be a dense matrix in v7.3 format);
  
    3. Click the "Run" button, or enter the following in the command line:
  
       ```
       test_example
       ```
  
    **Expected Output**:
  
    - Matrix symmetry judgment result (including Frobenius norm difference);
    - Time consumption statistics for a single run;
    - Average running time of multiple loops;
    - If `iscorrect='open'` is enabled, the decomposition error (1-norm) and correctness verification result will be output.

    ### 5.2 Custom Call Examples
  
    #### Example 1: Decomposition of Symmetric Banded Matrix (Specified Bandwidth)
  
    ```
    % Construct a symmetric banded matrix of order 2^13 (bandwidth h=6)
    n = 2^13; 
    h = 6; 
    A = sparse(n,n);
    for i=1:n
        j = max(1,i-h):min(n,i+h);
        A(i,j) = randn(1,length(j));
    end
    A = (A+A')/2;
    A = full(A);  % Must be converted to a dense matrix
    
    % Call decomposition (specify bandwidth h=6, decimal output, symmetric matrix, enable validation)
    [a, D] = sppaoli(6, A, 'dec', 'symmetric', 'open');
    ```
  
    #### Example 2: Decomposition of Asymmetric Matrix (Internal Automatic Bandwidth Calculation)
  
    ```
    % Construct a random asymmetric matrix of order 2^4
    n = 4; 
    N = 2^n;
    A = rand(N, N);  % Directly a dense matrix
    
    % Call decomposition (internal bandwidth calculation, quaternary output, asymmetric matrix, disable validation)
    [a, D] = sppaoli('off', A, 'qua', 'asymmetric', 'off');
    ```
  
    ### 5.3 Precautions
  
    1. **Matrix Order Limitation**: The order of the matrix A to be decomposed must be a power of 2 (e.g., 2, 4, 8, 16…), otherwise an error will be reported. If the matrix order does not meet the requirement, zero padding or truncation is required first. The maximum tested matrix order is up to 2^15.
    2. **Input Matrix Format Requirement**: The matrix A to be decomposed must be a **dense (full) matrix**. Direct input of sparse matrices is not supported, otherwise a parameter error will be triggered.
    3. **Parameter Format Specification**: `is_symmetric` must be passed as the string `'symmetric'` or `'asymmetric'`. Boolean values (such as `true`/`false`) or other strings are not allowed, otherwise a parameter verification error will be triggered.
    4. **MEX Compilation Issues**: If compilation fails, first check whether the MinGW64 compiler is correctly configured via `mex -setup`; if an error related to hardware instruction sets is encountered, compile with the `-march=native` parameter.
    5. **Time Consumption Reminder for Validation Function**: When the matrix scale is large (e.g., n>12), enabling `iscorrect='open'` will lead to a significant increase in time consumption due to matrix reconstruction. It is recommended to enable validation only for small matrix testing.
    6. **Bandwidth Parameter Optimization**: If the exact bandwidth of the matrix is known, it is recommended to directly pass the value of `h`, which can skip the internal bandwidth calculation step and further improve the decomposition efficiency.
  
    ------
  
    ## Document Update Log
  
    - **v1.0** (March 12, ): Initial version, covering core function description and user guide.