#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// ########################### MSVC Compatibility Fix ###########################
// Replace GCC built-in functions for Windows MSVC
#ifdef _MSC_VER
#include <intrin.h>  // Header for MSVC built-in functions
#define __builtin_popcount(x) _mm_popcnt_u32((unsigned int)x)
#define __builtin_popcountll(x) _mm_popcnt_u64((unsigned __int64)x)
#define __builtin_clz(x) _lzcnt_u32((unsigned int)x)
#define __builtin_clzll(x) _lzcnt_u64((unsigned __int64)x)
#endif
// #################################################################

// 1. Define mwSize first
typedef size_t mwSize;

// 2. Predefine double_complex
typedef struct {
    double real;
    double imag;
} double_complex;
#define MAKE_COMPLEX(r, i) ((double_complex){r, i})

// 3. Declare all functions
static PyObject* sppaolidec_w1(PyObject* self, PyObject* args);
static int is_pow2(int n);
static double complex_abs(double_complex c);
int binary_vector_dot_mwsize(mwSize i, mwSize j);
mwSize binary_to_quaternary_mwsize(mwSize j);
static void fht(double *x, int n);
static void t2f(const uint32_t *D, mwSize m, mwSize n, uint8_t *M);
static double get_elapsed_ms(clock_t start, clock_t end);
static void get_bandwidth(const double *pr, mwSize m, mwSize n, int *Lb, int *Ub);

// Method list
static PyMethodDef SppaoliMethods[] = {
    {"sppaolidec_w1", sppaolidec_w1, METH_VARARGS, "Pauli Decomposition (Asymmetric Matrix)"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef sppamodule = {
    PyModuleDef_HEAD_INIT,
    "paoli_is",
    "Asymmetric Matrix Pauli Decomposition Module",
    -1,
    SppaoliMethods
};

// Check if a number is a power of 2
static int is_pow2(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Calculate complex magnitude
static double complex_abs(double_complex c) {
    return sqrt(c.real * c.real + c.imag * c.imag);
}

// Binary vector dot product (compatible with Linux GCC / Windows MSVC)
int binary_vector_dot_mwsize(mwSize i, mwSize j) {
    mwSize common_ones = i & j;
    if (sizeof(mwSize) == 8) {
        return __builtin_popcountll((uint64_t)common_ones);  // 64-bit
    } else {
        return __builtin_popcount((uint32_t)common_ones);    // 32-bit
    }
}

// Binary to quaternary conversion (compatible with MSVC)
mwSize binary_to_quaternary_mwsize(mwSize j) {
    if (j == 0) {
        return 0;
    }
    unsigned long highest_bit;
    if (sizeof(mwSize) == 8) {  // 64-bit system
        highest_bit = 63 - __builtin_clzll((uint64_t)j);
    } else {  // 32-bit system
        highest_bit = 31 - __builtin_clz((uint32_t)j);
    }
    mwSize result = 0;
    for (unsigned long k = highest_bit; k != (unsigned long)-1; --k) {
        mwSize bit = (j >> k) & 1;
        result = (result << 2) + bit;
    }
    return result;
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

// Decimal to quaternary vector conversion
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

// Calculate elapsed time in milliseconds
static double get_elapsed_ms(clock_t start, clock_t end) {
    return ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
}

// Calculate matrix bandwidth
static void get_bandwidth(const double *pr, mwSize m, mwSize n, int *Lb, int *Ub) {
    *Lb = (int)m;
    *Ub = 0;

    if (n > 0) {
        mwSize last_col_first_row_idx = (n-1)*m + 0;
        if (fabs(pr[last_col_first_row_idx]) > 1e-13) {
            *Ub = (int)m - 1;
        }
    }

    if (m > 0) {
        mwSize last_row_first_col_idx = 0*m + (m-1);
        if (fabs(pr[last_row_first_col_idx]) > 1e-13) {
            *Lb = (int)n - 1;
        }
    }

    if (*Ub == (int)m - 1 && *Lb == (int)n - 1) {
        return;
    }

    for (mwSize j = 0; j < n; j++) {
        if (*Ub == (int)m - 1 && j == n - 1) {
            continue;
        }
        for (mwSize i = 0; i < m; i++) {
            if (*Lb == (int)n - 1 && i == m - 1) {
                continue;
            }
            if (fabs(pr[j * m + i]) > 1e-13) {
                int di = (int)i - (int)j;
                if (di < *Lb) *Lb = di;
                if (di > *Ub) *Ub = di;
            }
        }
    }
}

// Core decomposition function
static PyObject* sppaolidec_w1(PyObject* self, PyObject* args) {
    PyObject *bandwidth_obj, *Kt_obj;
    int type_;

    if (!PyArg_ParseTuple(args, "OOi", &bandwidth_obj, &Kt_obj, &type_)) {
        return NULL;
    }

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

    if (type_ < 0 || type_ > 2) {
        PyErr_SetString(PyExc_ValueError, "type must be 0, 1, or 2 (2=off)");
        return NULL;
    }
    if (type_ == 2) {
        printf("skip non-zero filtering and output nothing\n");
    }

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

    mwSize *in_k = (mwSize*)malloc(N * sizeof(mwSize));
    if (!in_k) {
        free(Kn4);
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed for in_k array");
        return NULL;
    }


    // 3. Generate in_k array
    mwSize bi = 0;
    for (int i_offset = -1; i_offset <= n - 1; i_offset++) {
        mwSize pow2i = (i_offset >= 0) ? (1 << i_offset) : 1;
        mwSize pow2j = (i_offset >= 0) ? (1 << i_offset) : 0;
        mwSize start_j = (pow2i > h) ? (pow2j - h + 1) : 1;
        mwSize end_j = pow2i;
        if (end_j < start_j) continue;
        for (mwSize j = start_j; j <= end_j; j++) {
            mwSize k = j + pow2j - 1;
            if (k >= N) continue;
            in_k[bi] = k;
            bi++;
            if (bi >= N) goto end_in_k_loop;
        }
        if (bi >= N) goto end_in_k_loop;
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
        mwSize k = in_k[idx];
        for (mwSize s = 0; s < N; s++) {
            mwSize id2 = k ^ s;
            DSD[idx * N + s] = Kt_pr[s * N + id2];
        }
    }
    // 6. FHT transformation and coefficient calculation
    double *a_re = (double*)malloc(bi * N * sizeof(double));
    double *a_im = (double*)malloc(bi * N * sizeof(double));
    uint32_t *DD_buf = (uint32_t*)malloc(bi * N * sizeof(uint32_t));
    if (!a_re || !a_im || !DD_buf) {
        free(NS3);
        free(Kn4);
        free(in_k);
        free(DSD);
        free(a_re);
        free(a_im);
        free(DD_buf);
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed for a/DD buffers");
        return NULL;
    }

    for (mwSize idx = 0; idx < bi; idx++) {
        mwSize k = in_k[idx];
        for (mwSize s = 0; s < N; s++) {
            DD_buf[idx * N + s] = (uint32_t)(Kn4[s] + 2 * Kn4[k]);
        }
    }

    for (mwSize idx = 0; idx < bi; idx++) {
        fht(&DSD[idx * N], (int)N);
        for (mwSize s = 0; s < N; s++) {
            int exponent = (int)NS3[idx * N + s];
            int sign = 1 - 2 * ((exponent >> 1) & 1);
            int y = (exponent & 1);
            a_re[idx * N + s] = DSD[idx * N + s] * (sign * (1 - y)) / N;
            a_im[idx * N + s] = DSD[idx * N + s] * (sign * y) / N;
        }
    }
    // 7. Filter non-zero elements and output
    PyObject *a_out = NULL;
    PyObject *DD_out = NULL;

    if (type_ == 2) {
        npy_intp empty_dim = 0;
        a_out = PyArray_SimpleNew(1, &empty_dim, NPY_COMPLEX128);
        DD_out = PyArray_SimpleNew(1, &empty_dim, NPY_UINT32);

    } else {
        mwSize valid_count = 0;
        for (mwSize i = 0; i < bi * N; i++) {
            if (a_re[i]*a_re[i] + a_im[i]*a_im[i] >= 1e-26) {
                valid_count++;
            }
        }

        npy_intp valid_count_npy = (npy_intp)valid_count;
        a_out = PyArray_SimpleNew(1, &valid_count_npy, NPY_COMPLEX128);
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
            free(a_re);
            free(a_im);
            free(DD_buf);
            PyErr_SetString(PyExc_MemoryError, "Failed to create output arrays");
            return NULL;
        }

        double *a_ptr = (double*)PyArray_DATA((PyArrayObject*)a_out);
        mwSize ptr = 0;
        if (type_ == 0) {
            uint8_t *dd_ptr = (uint8_t*)PyArray_DATA((PyArrayObject*)DD_out);
            for (mwSize i = 0; i < bi * N; i++) {
                if (a_re[i]*a_re[i] + a_im[i]*a_im[i] >= 1e-26) {
                    a_ptr[ptr * 2] = a_re[i];
                    a_ptr[ptr * 2 + 1] = a_im[i];
                    t2f(&DD_buf[i], 1, n, &dd_ptr[ptr * n]);
                    ptr++;
                }
            }
        } else {
            uint32_t *dd_ptr = (uint32_t*)PyArray_DATA((PyArrayObject*)DD_out);
            for (mwSize i = 0; i < bi * N; i++) {
                if (a_re[i]*a_re[i] + a_im[i]*a_im[i] >= 1e-26) {
                    a_ptr[ptr * 2] = a_re[i];
                    a_ptr[ptr * 2 + 1] = a_im[i];
                    dd_ptr[ptr] = DD_buf[i];
                    ptr++;
                }
            }
        }
    }


    free(NS3);
    free(Kn4);
    free(in_k);
    free(DSD);
    free(a_re);
    free(a_im);
    free(DD_buf);

    PyObject *result = PyTuple_New(2);
    PyTuple_SetItem(result, 0, a_out);
    PyTuple_SetItem(result, 1, DD_out);
    return result;
}

// Module initialization
PyMODINIT_FUNC PyInit_paoli_is(void) {
    import_array();
    return PyModule_Create(&sppamodule);
}