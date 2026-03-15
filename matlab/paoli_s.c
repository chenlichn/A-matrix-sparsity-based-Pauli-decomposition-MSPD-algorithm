#include "mex.h"
#include "matrix.h"
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <intrin.h>
#include <stdio.h>

// Convert mwSize integer j's binary vector to specified quaternary number
mwSize binary_to_quaternary_mwsize(mwSize j) {
    if (j == 0) {
        return 0; // Special case: 0 converts to 0
    }

    // Step 1: Locate the highest binary bit quickly (platform-dependent instructions)
    unsigned long highest_bit;
#ifdef _WIN64
    // 64-bit platform: mwSize is 64-bit, use 64-bit leading zero count
    unsigned long leading_zeros = __lzcnt64(j);
    highest_bit = 63 - leading_zeros; // Max index for 64-bit is 63
#else
    // 32-bit platform: mwSize is 32-bit, use 32-bit leading zero count
    unsigned long leading_zeros = __lzcnt(j);
    highest_bit = 31 - leading_zeros; // Max index for 32-bit is 31
#endif

    // Step 2: Compute quaternary result from highest to lowest bit (shift instead of multiply)
    mwSize result = 0;
    for (unsigned long k = highest_bit; k != (unsigned long)-1; --k) {
        // Extract k-th bit (0 or 1)
        mwSize bit = (j >> k) & 1;
        // Equivalent to result = result * 4 + bit (shift left 2 bits = ¡Á4)
        result = (result << 2) + bit;
    }

    return result;
}

// Compute dot product of binary vectors of two integers
int binary_vector_dot_mwsize(mwSize i, mwSize j) {
    // 1. Bitwise AND: keep positions where both bits are 1
    mwSize common_ones = i & j;
    
    // 2. Select hardware-accelerated function based on platform
#ifdef _WIN64
    // 64-bit MATLAB: mwSize is 64-bit, use 64-bit popcount instruction
    return _mm_popcnt_u64((uint64_t)common_ones);
#else
    // 32-bit MATLAB: mwSize is 32-bit, use 32-bit popcount instruction
    return _mm_popcnt_u32((uint32_t)common_ones);
#endif
}

// Check if n is a power of 2
static int is_pow2(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Fast Hadamard Transform (in-place operation)
static void fht(double *x, int n) {
    for (int m = 1; m < n; m *= 2) {
        for (int i = 0; i < n; i += 2 * m) {
            for (int j = 0; j < m; j++) {
                double a = x[i + j];
                double b = x[i + j + m];
                x[i + j] = a + b;
                x[i + j + m] = a - b;
            }
        }
    }
}

// Compute lower and upper bandwidth of a dense matrix
static void get_bandwidth(const mxArray *mat, int *Lb, int *Ub) {
    mwSize m = mxGetM(mat);
    mwSize n = mxGetN(mat);

    *Ub = 0;

    // Process dense matrix only, get data pointer directly
    const double *pr = mxGetPr(mat);
    
    for (mwIndex j = 0; j < n; j++) {
        // Check if last row of current column has non-zero element
        if (fabs(pr[j*m + (m-1)]) > 1e-13) {
            *Ub = m - 1 - j;
            *Lb = *Ub;
            return;
        }

        // Check all non-zeros in current column to update upper bandwidth
        for (mwIndex i = 0; i < m; i++) {
            if (fabs(pr[j*m + i]) > 1e-13) {
                int di = i - j;
                if (di > *Ub) {
                    *Ub = di;
                }
            }
        }
    }

    // For symmetric matrix, lower bandwidth equals upper bandwidth
    *Lb = *Ub;
}

// Decimal to quaternary conversion
static void t2f(const uint32_t *D, mwSize m, mwSize n, uint8_t *M) {
    uint8_t shifts[16];
    for (mwSize j = 0; j < n; j++) {
        shifts[j] = 2 * (n - 1 - j);
    }
    
    for (mwSize j = 0; j < n; j++) {
        uint8_t shift = shifts[j];
        for (mwSize i = 0; i < m; i++) {
            M[j * m + i] = (uint8_t)((D[i] >> shift) & 0x3);
        }
    }
}

// Helper function to calculate elapsed time in milliseconds
static double get_elapsed_ms(clock_t start, clock_t end) {
    return ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
}

// MEX gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Input argument check
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("paoli_s:inputCount", "3 input arguments required: h, Kt, type");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("paoli_s:outputCount", "2 output arguments required: a and DD");
    }
    
    // Parse first input: h (supports numeric or 'off' string)
    int input_bandwidth = -1; // Default: -1 = compute internally
    if (mxIsChar(prhs[0])) {
        // Check if string is 'off'
        char *bandwidth_str = mxArrayToString(prhs[0]);
        if (strcmp(bandwidth_str, "off") == 0) {
            input_bandwidth = -1;
            mexPrintf("Bandwidth set to 'off', will compute internally\n");
        } else {
            mexErrMsgIdAndTxt("paoli_s:invalidBandwidthStr", "h only supports string 'off'");
        }
        mxFree(bandwidth_str);
    } else if (mxIsNumeric(prhs[0])) {
        // Parse numeric bandwidth
        input_bandwidth = (int)mxGetScalar(prhs[0]);
        if (input_bandwidth <= 0) {
            mexErrMsgIdAndTxt("paoli_s:invalidBandwidthNum", "h must be positive integer");
        }
        mexPrintf("Using input bandwidth value: %d\n", input_bandwidth);
    } else {
        mexErrMsgIdAndTxt("paoli_s:invalidBandwidthType", "h must be numeric or string 'off'");
    }
    
    // Parse second input: Kt matrix
    const mxArray *Kt = prhs[1];
    if (mxIsComplex(Kt)) {
        mexErrMsgIdAndTxt("paoli_s:complexInput", "Only real matrix supported");
    }
    mwSize N = mxGetM(Kt);
    mwSize N_cols = mxGetN(Kt);
    if (N_cols != N) {
        mexErrMsgIdAndTxt("paoli_s:notSquare", "Input must be square matrix");
    }
    
    // Parse third input: type (0/1/2 supported, 2=off)
    int type = (int)mxGetScalar(prhs[2]);
    if (type < 0 || type > 2) {
        mexErrMsgIdAndTxt("paoli_s:invalidType", "type must be 0, 1, or 2 (2=off)");
    }
    if (type == 2) {
        mexPrintf("skip non-zero filtering and output nothing\n");
    }

    // Compute n = log2(N)
    int n = 0;
    mwSize temp = N;
    while (temp > 1) {
        temp >>= 1;
        n++;
    }
    if ((1 << n) != N) {
        mexErrMsgIdAndTxt("paoli_s:notPowerOf2", "Matrix order must be power of 2");
    }

    // 1. Compute bandwidth
    clock_t t_start, t_end;
    int Lb, Ub;
    int h;
    if (input_bandwidth > 0) {
        // Use input bandwidth, skip internal calculation
        h = input_bandwidth;
    } else {
        // Compute bandwidth internally
        t_start = clock();
        get_bandwidth(Kt, &Lb, &Ub);
        t_end = clock();
        h = (abs(Lb) > abs(Ub)) ? abs(Lb) : abs(Ub);
        mexPrintf("[Timer] Bandwidth calculation time: %.8f ms, h=%d\n", get_elapsed_ms(t_start, t_end), h);
    }
    

    // 2. Compute n4 and Kn4
    mwSize *Kn4 = (mwSize *)mxMalloc(N * sizeof(mwSize));
    for (mwSize s = 0; s < N; s++) {
        Kn4[s] = binary_to_quaternary_mwsize(s);
    }


    // 3. Generate in_k array
    mwSize bi = 0;
    mwSize *in_k = (mwSize *)mxMalloc(N * sizeof(mwSize));
    
    for (int i = -1; i <= (int)n - 1; i++) {
        mwSize pow2i = (i >= 0) ? (1 << i) : 1;
        mwSize pow2j = (i >= 0) ? (1 << i) : 0;
        mwSize start_j = (pow2i > h) ? (pow2j - h + 1) : 1;
        mwSize end_j = pow2i;

        if (end_j < start_j) continue;
        
        for (mwSize j = start_j; j <= end_j; j++) {
            mwSize k = j + pow2j - 1;
            if (k >= N) continue;
            in_k[bi] = k;
            bi++;
            if (bi > N) goto end_in_k_loop;
        }
    }
end_in_k_loop:
    in_k = (mwSize *)mxRealloc(in_k, bi * sizeof(mwSize));


    // 4. Compute NS3 matrix
    mwSize *NS3 = (mwSize *)mxMalloc(bi * N * sizeof(mwSize));
    for (mwSize idx = 0; idx < bi; idx++) {
        mwSize i = in_k[idx];
        for (mwSize j = 0; j < N; j++) {
            NS3[idx * N + j] = binary_vector_dot_mwsize(i, j);
        }
    }


    // 5. Compute DSD matrix
    double *DSD = (double *)mxMalloc(N * bi * sizeof(double));
    const double *Kt_pr = mxGetPr(Kt);
    
    for (mwSize idx = 0; idx < bi; idx++) {
        for (mwSize s = 0; s < N; s++) {
            mwSize k = in_k[idx];
            mwSize id2 = (k ^ s);
            mwSize qd2 = id2 * N + s;
            DSD[idx * N + s] = Kt_pr[qd2];
        }
    }


    // 6. FHT transform and coefficient calculation
    double *a = (double *)mxMalloc(bi * N * sizeof(double));
    uint32_t *DD = (uint32_t *)mxMalloc(bi * N * sizeof(uint32_t));
    
    // Fill DD
    for (mwSize idx = 0; idx < bi; idx++) {
        mwSize k = in_k[idx];
        for (mwSize s = 0; s < N; s++) {
            DD[idx * N + s] = Kn4[s] + 2 * Kn4[k];
        }
    }
    
    // Perform FHT and generate a
    for (mwSize idx = 0; idx < bi; idx++) {
        double *row = &DSD[idx * N];
        fht(row, N);
        
        for (mwSize s = 0; s < N; s++) {
            int exponent = (int)NS3[idx * N + s];
            int sign = 1 - 2 * ((exponent >> 1) & 1);
            a[idx * N + s] = row[s] * sign / N;
        }
    }

    // 7. Filter non-zero elements and output
    if (type == 2) {
        // type=off: output empty matrices
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        plhs[1] = mxCreateNumericMatrix(0, 0, mxUINT32_CLASS, mxREAL);
    } else {
        // Count valid non-zero elements
        mwSize valid_count = 0;
        for (mwSize i = 0; i < bi * N; i++) {
            if (fabs(a[i]) >= 1e-13) valid_count++;
        }

        plhs[0] = mxCreateDoubleMatrix(1, valid_count, mxREAL);
        double *a_out = mxGetPr(plhs[0]);
        
        mxArray *DD_out;
        if (type == 0) {
            DD_out = mxCreateNumericMatrix(n, valid_count, mxUINT8_CLASS, mxREAL);
        } else {
            DD_out = mxCreateNumericMatrix(1, valid_count, mxUINT32_CLASS, mxREAL);
        }
        plhs[1] = DD_out;
        
        mwSize ptr = 0;
        for (mwSize i = 0; i < bi * N; i++) {
            if (fabs(a[i]) >= 1e-13) {
                a_out[ptr] = a[i];
                
                if (type == 0) {
                    uint8_t *dd_ptr = (uint8_t *)mxGetData(DD_out);
                    t2f(&DD[i], 1, n, &dd_ptr[ptr * n]);
                } else {
                    uint32_t *dd_ptr = (uint32_t *)mxGetData(DD_out);
                    dd_ptr[ptr] = DD[i];
                }
                ptr++;
            }
        }
    }

    // Free allocated memory
    mxFree(NS3);
    mxFree(Kn4);
    mxFree(in_k);
    mxFree(DSD);
    mxFree(a);
    mxFree(DD);
}