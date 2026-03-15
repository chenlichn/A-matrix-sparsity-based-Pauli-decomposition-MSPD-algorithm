#include "mex.h"
#include "matrix.h"
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <intrin.h>
#include <stdio.h>

// Convert mwSize integer j from binary vector to quaternary number
mwSize binary_to_quaternary_mwsize(mwSize j) {
    if (j == 0) {
        return 0;
    }

    unsigned long highest_bit;
#ifdef _WIN64
    unsigned long leading_zeros = __lzcnt64(j);
    highest_bit = 63 - leading_zeros;
#else
    unsigned long leading_zeros = __lzcnt(j);
    highest_bit = 31 - leading_zeros;
#endif

    mwSize result = 0;
    for (unsigned long k = highest_bit; k != (unsigned long)-1; --k) {
        mwSize bit = (j >> k) & 1;
        result = (result << 2) + bit;
    }

    return result;
}

// Compute dot product of binary vectors of two integers
int binary_vector_dot_mwsize(mwSize i, mwSize j) {
    mwSize common_ones = i & j;

#ifdef _WIN64
    return _mm_popcnt_u64((uint64_t)common_ones);
#else
    return _mm_popcnt_u32((uint32_t)common_ones);
#endif
}

// Check if n is a power of 2
static int is_pow2(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Fast Hadamard Transform (in-place, real implementation; extendable for complex FHT)
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

static void get_bandwidth(const mxArray *mat, int *Lb, int *Ub) {
    mwSize m = mxGetM(mat);
    mwSize n = mxGetN(mat);
    *Lb = m;
    *Ub = 0;

    // Process dense matrix only, get data pointer directly
    const double *pr = mxGetPr(mat);

    // Check last column, first row element (i=0, j=n-1)
    if (n > 0) {
        mwIndex last_col_first_row_idx = (n-1)*m + 0;
        if (fabs(pr[last_col_first_row_idx]) > 1e-13) {
            *Ub = m - 1;
        }
    }

    // Check first column, last row element (i=m-1, j=0)
    if (m > 0) {
        mwIndex last_row_first_col_idx = 0*m + (m-1);
        if (fabs(pr[last_row_first_col_idx]) > 1e-13) {
            *Lb = n - 1;
        }
    }

    // Return early if both bandwidths reach maximum
    if (*Ub == m - 1 && *Lb == n - 1) {
        return;
    }

    // Traverse all elements to compute bandwidth
    for (mwIndex j = 0; j < n; j++) {
        if (*Ub == m - 1 && j == n - 1) {
            continue;
        }

        for (mwIndex i = 0; i < m; i++) {
            if (*Lb == n - 1 && i == m - 1) {
                continue;
            }

            if (fabs(pr[j*m + i]) > 1e-13) {
                int di = i - j;
                if (di < *Lb) {
                    *Lb = di;
                }
                if (di > *Ub) {
                    *Ub = di;
                }
            }
        }
    }
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
    // Input arguments check
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("sppaolidec_w1:inputCount", "3 input arguments required: h, Kt, type");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("sppaolidec_w1:outputCount", "2 output arguments required: a and DD");
    }

    // Parse first input: h (numeric or 'off' string)
    int input_bandwidth = -1;
    if (mxIsChar(prhs[0])) {
        char *bandwidth_str = mxArrayToString(prhs[0]);
        if (strcmp(bandwidth_str, "off") == 0) {
            input_bandwidth = -1;
            mexPrintf("Bandwidth set to 'off', will compute internally\n");
        } else {
            mexErrMsgIdAndTxt("sppaolidec_w1:invalidBandwidthStr", "h only supports string 'off'");
        }
        mxFree(bandwidth_str);
    } else if (mxIsNumeric(prhs[0])) {
        input_bandwidth = (int)mxGetScalar(prhs[0]);
        if (input_bandwidth <= 0) {
            mexErrMsgIdAndTxt("sppaolidec_w1:invalidBandwidthNum", "h must be a positive integer");
        }
        mexPrintf("Using input bandwidth value: %d\n", input_bandwidth);
    } else {
        mexErrMsgIdAndTxt("sppaolidec_w1:invalidBandwidthType", "h must be numeric or string 'off'");
    }

    // Parse second input: Kt matrix
    const mxArray *Kt = prhs[1];
    if (mxIsComplex(Kt)) {
        mexErrMsgIdAndTxt("sppaolidec_w1:complexInput", "Only real matrix input is supported");
    }
    mwSize N = mxGetM(Kt);
    mwSize N_cols = mxGetN(Kt);
    if (N_cols != N) {
        mexErrMsgIdAndTxt("sppaolidec_w1:notSquare", "Input matrix must be square");
    }

    // Parse third input: type (0,1,2; 2=off)
    int type = (int)mxGetScalar(prhs[2]);
    if (type < 0 || type > 2) {
        mexErrMsgIdAndTxt("sppaolidec_w1:invalidType", "type must be 0, 1, or 2 (2=off)");
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
        mexErrMsgIdAndTxt("sppaolidec_w1:notPowerOf2", "Matrix order must be power of 2");
    }

    // 1. Compute bandwidth
    clock_t t_start, t_end;
    int Lb, Ub;
    int h;
    if (input_bandwidth > 0) {
        h = input_bandwidth;
    } else {
        t_start = clock();
        get_bandwidth(Kt, &Lb, &Ub);
        t_end = clock();
        h = (abs(Lb) > abs(Ub)) ? abs(Lb) : abs(Ub);
        mexPrintf("[Timer] Bandwidth calculation time: %.8f ms, h=%d\n", get_elapsed_ms(t_start, t_end), h);
    }

    // 2. Compute Kn4
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
    double *a_re = (double *)mxMalloc(bi * N * sizeof(double));
    double *a_im = (double *)mxMalloc(bi * N * sizeof(double));
    uint32_t *DD = (uint32_t *)mxMalloc(bi * N * sizeof(uint32_t));

    // Fill DD
    for (mwSize idx = 0; idx < bi; idx++) {
        mwSize k = in_k[idx];
        for (mwSize s = 0; s < N; s++) {
            DD[idx * N + s] = Kn4[s] + 2 * Kn4[k];
        }
    }

    // Perform FHT and generate coefficients
    for (mwSize idx = 0; idx < bi; idx++) {
        double *row = &DSD[idx * N];
        fht(row, N);

        for (mwSize s = 0; s < N; s++) {
            int exponent = (int)NS3[idx * N + s];
            int sign = 1 - 2 * ((exponent >> 1) & 1);
            int y = (exponent & 1);
            a_re[idx * N + s] = row[s] * (sign * (1-y)) / N;
            a_im[idx * N + s] = row[s] * (sign * y) / N;
        }
    }

    // 7. Filter non-zero elements and output
    if (type == 2) {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
        plhs[1] = mxCreateNumericMatrix(0, 0, mxUINT32_CLASS, mxREAL);
    } else {
        mwSize valid_count = 0;
        for (mwSize i = 0; i < bi * N; i++) {
            if (a_re[i]*a_re[i] + a_im[i]*a_im[i] >= 1e-26) {
                valid_count++;
            }
        }

        plhs[0] = mxCreateDoubleMatrix(1, valid_count, mxCOMPLEX);
        double *a_out_re = mxGetPr(plhs[0]);
        double *a_out_im = mxGetPi(plhs[0]);

        mxArray *DD_out;
        if (type == 0) {
            DD_out = mxCreateNumericMatrix(n, valid_count, mxUINT8_CLASS, mxREAL);
        } else {
            DD_out = mxCreateNumericMatrix(1, valid_count, mxUINT32_CLASS, mxREAL);
        }
        plhs[1] = DD_out;

        mwSize ptr = 0;
        for (mwSize i = 0; i < bi * N; i++) {
            if (a_re[i]*a_re[i] + a_im[i]*a_im[i] >= 1e-26) {
                a_out_re[ptr] = a_re[i];
                a_out_im[ptr] = a_im[i];

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
    mxFree(a_re);
    mxFree(a_im);
    mxFree(DD);
}