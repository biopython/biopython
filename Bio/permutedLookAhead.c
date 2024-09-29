#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h> // for FLT_MIN

#define NEG_INF -INFINITY

// Define a struct for holding column scores and their original indices
// This struct is used to keep track of the scores for each nucleotide (A, C, G, T) and the original index of the column
typedef struct {
    float scores[4];       // Scores for nucleotides 'A', 'C', 'G', 'T'
    int original_index;    // Original index of the column in the scoring matrix
} Column;

// Comparison function for sorting columns by their maximum score
// This function is used by qsort to sort columns in descending order based on their maximum scores
int compare_columns(const void *a, const void *b) {
    float max_a = FLT_MIN;
    float max_b = FLT_MIN;
    const Column *col_a = (const Column *)a;
    const Column *col_b = (const Column *)b;
    for (int i = 0; i < 4; i++) {
        if (col_a->scores[i] > max_a) max_a = col_a->scores[i];
        if (col_b->scores[i] > max_b) max_b = col_b->scores[i];
    }
    return (max_b > max_a) - (max_b < max_a); // Return -1 if max_b > max_a, 1 if max_b < max_a, 0 if equal (for descending order)
}

/*  
    Implements permuted look ahead method for scoring a given DNA sequence s of length l
    and a scoring matrix of width m, reporting sites above a given threshold (t).

    Receives as an args object the sequence (string), the scoring matrix (PyObject),
    the width of the matrix (m) and the score sequence threshold (t).
    
    The permuted lookahead method for PSSM scoring provides a speedup to the original
    PSSM scoring procedure by rearranging the matrix rows to prioritize higher intermediate 
    thresholds and aborting the evaluation of a sequence segment when it cannot fulfill 
    the required threshold with the remaining column scores.
*/
static PyObject* permuted_look_ahead(PyObject *self, PyObject *args) { 
    const char *s;
    PyObject *scoring_matrix_obj;
    int m;
    float t;

    // Check given args looking for string, object, float and integer ("sOfi")
    if (!PyArg_ParseTuple(args, "sOfi", &s, &scoring_matrix_obj, &t, &m)) {
        return NULL;
    }

    // Define the length of the DNA sequence
    size_t l = strlen(s);

    // Allocate dynamic memory for columns (array of dimension m)
    Column *columns = (Column *)malloc(m * sizeof(Column));
    if (!columns) return PyErr_NoMemory();

    // Allocate and fill columns from scoring matrix
    // i is the column index (range m)
    // j is the row index (range 4, as it represents "ACGT")
    for (int i = 0; i < m; i++) {
        PyObject *row = PyList_GetItem(scoring_matrix_obj, i);
        if (!PyList_Check(row)) {
            free(columns);
            return NULL;
        }
        for (int j = 0; j < 4; j++) {
            columns[i].scores[j] = (float)PyFloat_AsDouble(PyList_GetItem(row, j));
        }
        columns[i].original_index = i;
    }

    // Sort columns based on maximum score using qsort and the comparison function
    qsort(columns, m, sizeof(Column), compare_columns);

    // Allocate memory for max_scores, Z, and T_intermediate
    float *max_scores = (float *)malloc(m * sizeof(float));
    float *Z = (float *)malloc(m * sizeof(float));
    float *T_intermediate = (float *)malloc(m * sizeof(float));

    if (!max_scores || !Z || !T_intermediate) {
        free(columns);
        if (max_scores) free(max_scores);
        if (Z) free(Z);
        if (T_intermediate) free(T_intermediate);
        return PyErr_NoMemory();
    }

    // Calculate max_scores for each column
    // This is an array that holds, at each position of the scoring matrix, the maximum score that can be obtained with the columns at that position
    for (int i = 0; i < m; i++) {
        max_scores[i] = FLT_MIN;
        for (int j = 0; j < 4; j++) {
            if (columns[i].scores[j] != NEG_INF && columns[i].scores[j] > max_scores[i]) {
                max_scores[i] = columns[i].scores[j];
            }
        }
    }
    // Calculate Z(i) values
    // This array holds the cumulative maximum possible score from each position to the end of the scoring matrix
    // Starts from the second-to-last position and goes backwards, then adds the maximum score of the next position to the current Z value
    Z[m-1] = 0.0; // Last position has no remaining score
    for (int i = m - 2; i >= 0; i--) {
        Z[i] = Z[i + 1] + max_scores[i + 1];
    }

    // Calculate intermediate thresholds T_intermediate
    // This array holds the intermediate thresholds at each position
    // Each threshold is the final threshold - the cumulative max score Z(i)
    for (int i = 0; i < m; i++) {
        T_intermediate[i] = t - Z[i];
    }

    // Create a Python object list to hold the search results
    PyObject *result_list = PyList_New(0);

    // Process each segment in sequence s, from sequence start to l-m
    for (int start = 0; start <= l - m; start++) {
        float total_score = 0.0;

        // Iterate over each position in the segment
        for (int j = 0; j < m; j++) {
            char nucleotide = s[start + columns[j].original_index];
            int score_index = 0;

            // Use switch-case to handle both uppercase and lowercase nucleotides
            switch (nucleotide) {
                case 'C': 
                case 'c': 
                    score_index = 1; 
                    break;
                case 'G': 
                case 'g': 
                    score_index = 2; 
                    break;
                case 'T': 
                case 't': 
                    score_index = 3; 
                    break;
                case 'A': 
                case 'a': 
                    score_index = 0; 
                    break;
                default: break;
            }

            // Add the score for the current nucleotide
            total_score += columns[j].scores[score_index];

            // Early exit if it's impossible to reach the threshold
            if (total_score < T_intermediate[j]) {
                break;
            }
        }

        // If the total score >= threshold, add the start index to the result list
        if (total_score >= t) {
            PyObject *start_index = PyLong_FromLong(start);
            PyList_Append(result_list, start_index);
            Py_DECREF(start_index);
        }
    }

    // Free allocated memory
    free(columns);
    free(max_scores);
    free(Z);
    free(T_intermediate);

    // Return the result list
    return result_list;
}

// Define the methods that will be available in the Python module
static PyMethodDef motif_methods_pl[] = {
    {"permuted_look_ahead", permuted_look_ahead, METH_VARARGS, "Find DNA motifs using permuted lookahead scoring."},
    {NULL, NULL, 0, NULL} // Sentinel value, indicate the end of the array
};

// Define the Python module itself
static struct PyModuleDef permuted_look_ahead_module = {
    PyModuleDef_HEAD_INIT, // Macro to initialize the module definition
    "permuted_look_ahead",  // Name of the module
    NULL,                  // Module documentation (NULL)
    -1,                    // -1 as the module keeps state in global variables
    motif_methods_pl       
};

// Initialization function for the module
// This function is called when the module is imported in Python
PyMODINIT_FUNC PyInit_permuted_look_ahead(void) {
    return PyModule_Create(&permuted_look_ahead_module); // Create the module and return a PyObject representing it
}

