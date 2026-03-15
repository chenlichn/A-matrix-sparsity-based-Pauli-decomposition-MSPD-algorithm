#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// ########################### MSVC Compatibility Fix ###########################
#ifdef _MSC_VER
#include <intrin.h>
#define __builtin_popcount(x) _mm_popcnt_u32((unsigned int)x)
#define __builtin_popcountll(x) _mm_popcnt_u64((unsigned __int64)x)
#define __builtin_clz(x) _lzcnt_u32((unsigned int)x)
#define __builtin_clzll(x) _lzcnt_u64((unsigned __int64)x)
#endif
// #################################################################

// 1. Define MATLAB compatible types
typedef size_t mwSize;
typedef size_t mwIndex;

// 2. Declare all functions
static PyObject* sppaolidec_w1(PyObject* self, PyObject* args);
static int is_pow2(int n);
static double get_elapsed_ms(clock_t start, clock_t end);
static void get_bandwidth(const double *pr, mwSize m, mwSize n, int *Lb, int *Ub);
mwSize binary_to_quaternary_mwsize(mwSize j);
int binary_vector_dot_mwsize(mwSize i, mwSize j);
static void fht(double *x, int n);
static void t2f(const uint32_t *D, mwSize m, mwSize n, uint8_t *M);

// Method list
static PyMethodDef SppaoliSMethods[] = {
    {"sppaolidec_w1", sppaolidec_w1, METH_VARARGS, "Pauli Decomposition (Symmetric Matrix)"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef sppasmodule = {
    PyModuleDef_HEAD_INIT,
    "paoli_s",
    "Symmetric Matrix Pauli Decomposition Module",
    -1,
    SppaoliSMethods
};

// Check if a number is a power of 2
static int is_pow2(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Calculate elapsed time in milliseconds
static double get_elapsed_ms(clock_t start, clock_t end) {
    return ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
}

// Convert binary vector of mwSize integer j to specified quaternary number
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

// Calculate dot product of binary vectors of two integers
int binary_vector_dot_mwsize(mwSize i, mwSize j) {
    mwSize common_ones = i & j;
#ifdef _WIN64
    return _mm_popcnt_u64((uint64_t)common_ones);
#else
    return _mm_popcnt_u32((uint32_t)common_ones);
#endif
}

// Fast Hadamard Transform
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

// Calculate bandwidth (optimized for symmetric matrices)
static void get_bandwidth(const double *pr, mwSize m, mwSize n, int *Lb, int *Ub) {
    *Ub = 0;  // Initialize upper bandwidth to 0

    for (mwIndex j = 0; j < n; j++) {
        // Check if there is a non-zero element in the last row of current column
        if (fabs(pr[j*m + (m-1)]) > 1e-13) {
            *Ub = m - 1 - j;
            *Lb = *Ub;  // For symmetric matrix, lower bandwidth equals upper bandwidth
            return;
        }

        // Check all non-zero elements in current column, update upper bandwidth
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

// Core decomposition function (Symmetric Matrix)
static PyObject* sppaolidec_w1(PyObject* self, PyObject* args) {
    PyObject *bandwidth_obj, *Kt_obj;
    int type_;

    if (!PyArg_ParseTuple(args, "OOi", &bandwidth_obj, &Kt_obj, &type_)) {
        return NULL;
    }

    // Parse bandwidth parameter
    int input_bandwidth = -1;
    if (PyUnicode_Check(bandwidth_obj)) {
        const char *bandwidth_str = PyUnicode_AsUTF8(bandwidth_obj);
        if (strcmp(bandwidth_str, "off") == 0) {
            input_bandwidth = -1;
            printf("Bandwidth parameter is 'off', bandwidth will be calculated internally\n");
        } else {
            PyErr_SetString(PyExc_ValueError, "h only supports the string 'off'");
            return NULL;
        }
    } else if (PyNumber_Check(bandwidth_obj)) {
        input_bandwidth = (int)PyLong_AsLong(bandwidth_obj);
        if (input_bandwidth <= 0) {
            PyErr_SetString(PyExc_ValueError, "h must be an integer greater than 0");
            return NULL;
        }
        printf("Using input bandwidth value: %d\n", input_bandwidth);
    } else {
        PyErr_SetString(PyExc_TypeError, "h must be a number or the string 'off'");
        return NULL;
    }

    // Parse Kt matrix
    PyArrayObject *Kt_arr = (PyArrayObject*)Kt_obj;
    if (!PyArray_Check(Kt_arr)) {
        PyErr_SetString(PyExc_TypeError, "Kt must be a numpy array");
        return NULL;
    }
    if (PyArray_TYPE(Kt_arr) != NPY_DOUBLE || PyArray_ISCOMPLEX(Kt_arr)) {
        PyErr_SetString(PyExc_TypeError, "Input matrix must be a double-precision real matrix");
        return NULL;
    }

    mwSize N = PyArray_DIM(Kt_arr, 0);
    mwSize N_cols = PyArray_DIM(Kt_arr, 1);
    if (N != N_cols) {
        PyErr_SetString(PyExc_ValueError, "Input matrix must be square");
        return NULL;
    }

    // Parse type parameter
    if (type_ < 0 || type_ > 2) {
        PyErr_SetString(PyExc_ValueError, "type must be 0, 1, or 2 (2=off)");
        return NULL;
    }
    if (type_ == 2) {
        printf("skip non-zero filtering and output nothing\n");
    }

    // Check if matrix order is a power of 2
    int n = 0;
    mwSize temp = N;
    while (temp > 1) {
        temp >>= 1;
        n++;
    }
    if ((1 << n) != N) {
        PyErr_SetString(PyExc_ValueError, "Matrix order must be a power of 2");
        return NULL;
    }

    double *Kt_pr = (double*)PyArray_DATA(Kt_arr);
    if (!Kt_pr) {
        PyErr_SetString(PyExc_MemoryError, "Failed to access matrix data");
        return NULL;
    }

    // 1. Calculate bandwidth
    clock_t t_start, t_end;
    int Lb, Ub;
    int h;
    if (input_bandwidth > 0) {
        h = input_bandwidth;
    } else {
        t_start = clock();
        get_bandwidth(Kt_pr, N, N, &Lb, &Ub);
        t_end = clock();
        h = (abs(Lb) > abs(Ub)) ? abs(Lb) : abs(Ub);
        printf("[Timer] Bandwidth calculation time: %.8f ms, result h=%d\n", get_elapsed_ms(t_start, t_end), h);
    }

    // 2. Compute Kn4 array
    mwSize *Kn4 = (mwSize*)malloc(N * sizeof(mwSize));
    if (!Kn4) {
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed for Kn4 array");
        return NULL;
    }
    for (mwSize s = 0; s < N; s++) {
        Kn4[s] = binary_to_quaternary_mwsize(s);
    }


    // 3. Generate in_k array
    mwSize bi = 0;
    mwSize *in_k = (mwSize*)malloc(N * sizeof(mwSize));
    if (!in_k) {
        free(Kn4);
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed for in_k array");
        return NULL;
    }

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
    in_k = (mwSize*)realloc(in_k, bi * sizeof(mwSize));
    if (!in_k) {
        free(Kn4);
        PyErr_SetString(PyExc_MemoryError, "Reallocation failed for in_k array");
        return NULL;
    }


    // 4. Compute NS3 matrix
    mwSize *NS3 = (mwSize*)malloc(bi * N * sizeof(mwSize));
    if (!NS3) {
        free(Kn4);
        free(in_k);
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed for NS3 matrix");
        return NULL;
    }

    for (mwSize idx = 0; idx < bi; idx++) {
        mwSize i = in_k[idx];
        for (mwSize j = 0; j < N; j++) {
            NS3[idx * N + j] = binary_vector_dot_mwsize(i, j);
        }
    }


    // 5. Compute DSD matrix
    double *DSD = (double*)malloc(N * bi * sizeof(double));
    if (!DSD) {
        free(NS3);
        free(Kn4);
        free(in_k);
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed for DSD matrix");
        return NULL;
    }

    for (mwSize idx = 0; idx < bi; idx++) {
        for (mwSize s = 0; s < N; s++) {
            mwSize k = in_k[idx];
            mwSize id2 = (k ^ s);
            mwSize qd2 = id2 * N + s;
            DSD[idx * N + s] = Kt_pr[qd2];
        }
    }


    // 6. FHT transformation and coefficient calculation
    double *a = (double*)malloc(bi * N * sizeof(double));
    uint32_t *DD = (uint32_t*)malloc(bi * N * sizeof(uint32_t));
    if (!a || !DD) {
        free(NS3);
        free(Kn4);
        free(in_k);
        free(DSD);
        free(a);
        free(DD);
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed for a/DD arrays");
        return NULL;
    }

    // Fill DD
    for (mwSize idx = 0; idx < bi; idx++) {
        mwSize k = in_k[idx];
        for (mwSize s = 0; s < N; s++) {
            DD[idx * N + s] = Kn4[s] + 2 * Kn4[k];
        }
    }

    // Compute FHT and generate a
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
    PyObject *a_out = NULL;
    PyObject *DD_out = NULL;

    if (type_ == 2) {
        // type=off: output empty matrix
        npy_intp empty_dim = 0;
        a_out = PyArray_SimpleNew(1, &empty_dim, NPY_DOUBLE);
        DD_out = PyArray_SimpleNew(1, &empty_dim, NPY_UINT32);
    } else {
        // Filter non-zero elements
        mwSize valid_count = 0;
        for (mwSize i = 0; i < bi * N; i++) {
            if (fabs(a[i]) >= 1e-13) valid_count++;
        }

        // Create output arrays
        npy_intp valid_count_npy = (npy_intp)valid_count;
        a_out = PyArray_SimpleNew(1, &valid_count_npy, NPY_DOUBLE);
        if (type_ == 0) {
            npy_intp dd_dims[2] = {(npy_intp)n, valid_count_npy};
            DD_out = PyArray_New(&PyArray_Type, 2, dd_dims, NPY_UINT8,
                               NULL, NULL, 0, NPY_FORTRANORDER, NULL);
        } else {
            DD_out = PyArray_SimpleNew(1, &valid_count_npy, NPY_UINT32);
        }

        if (!a_out || !DD_out) {
            Py_XDECREF(a_out);
            Py_XDECREF(DD_out);
            free(NS3);
            free(Kn4);
            free(in_k);
            free(DSD);
            free(a);
            free(DD);
            PyErr_SetString(PyExc_MemoryError, "Failed to create output arrays");
            return NULL;
        }

        // Fill output data
        double *a_ptr = (double*)PyArray_DATA((PyArrayObject*)a_out);
        mwSize ptr = 0;

        if (type_ == 0) {
            uint8_t *dd_ptr = (uint8_t*)PyArray_DATA((PyArrayObject*)DD_out);
            for (mwSize i = 0; i < bi * N; i++) {
                if (fabs(a[i]) >= 1e-13) {
                    a_ptr[ptr] = a[i];
                    t2f(&DD[i], 1, n, &dd_ptr[ptr * n]);
                    ptr++;
                }
            }
        } else {
            uint32_t *dd_ptr = (uint32_t*)PyArray_DATA((PyArrayObject*)DD_out);
            for (mwSize i = 0; i < bi * N; i++) {
                if (fabs(a[i]) >= 1e-13) {
                    a_ptr[ptr] = a[i];
                    dd_ptr[ptr] = DD[i];
                    ptr++;
                }
            }
        }
    }

    // Release memory
    free(NS3);
    free(Kn4);
    free(in_k);
    free(DSD);
    free(a);
    free(DD);

    // Construct return result
    PyObject *result = PyTuple_New(2);
    PyTuple_SetItem(result, 0, a_out);
    PyTuple_SetItem(result, 1, DD_out);
    return result;
}

// Module initialization
PyMODINIT_FUNC PyInit_paoli_s(void) {
    import_array();
    return PyModule_Create(&sppasmodule);
}