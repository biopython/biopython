#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h> // for FLT_MIN

#define NEG_INF -INFINITY

/*
    Implements look ahead method for scoring a given DNA sequence s of length l
    and a scoring matrix of width m, reporting sites above a given threshold (t)

    Receives as an args object the sequence (string), the scoring matrix (PyObject),
    the width of the matrix (m) and the score sequence threshold (t).
    
    The lookahead method for PSSM scoring provides a speedup to the original
    PSSM scoring procedure by aborting the evaluation of a sequence segment when
    it cannot fulfill the required threshold with the remaining column scores.
*/
static PyObject* look_ahead(PyObject *self, PyObject *args) {
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
  
    // Allocate dynamic memory for scoring matrix (row of dimension m)
    float **scoring_matrix = (float **)malloc(m * sizeof(float *));
    if (!scoring_matrix) return PyErr_NoMemory();
    
    // Allocate dynamic memory for scoring matrix (column of dimension 4)
    // i is the column index (range m)
    // j is the row index (range 4, as it represents "ACGT")
    for (int i = 0; i < m; i++) {
        scoring_matrix[i] = (float *)malloc(4 * sizeof(float));
        if (!scoring_matrix[i]) {
            // If allocation fails, free allocated memory
            for (int j = 0; j < i; j++) {
                free(scoring_matrix[j]);
            }
            free(scoring_matrix);
            return PyErr_NoMemory();
        }

        // Get row data for this particular column from Python scoring matrix object
        PyObject *row = PyList_GetItem(scoring_matrix_obj, i);
        if (!PyList_Check(row)) {
            // If we don't have a correct entry, free allocated memory
            for (int j = 0; j <= i; j++) {
                free(scoring_matrix[j]);
            }
            free(scoring_matrix);
            return NULL;
        }
        // Fill up column with the row data
        for (int j = 0; j < 4; j++) {
            scoring_matrix[i][j] = (float)PyFloat_AsDouble(PyList_GetItem(row, j));
        }
    }

    // Calculate initial max_possible_score at each position
    // This is an array that holds, at each position of the scoring
    // matrix, the maximum score that can be obtained 
    // with the columns at that position
    float* max_scores = (float*)malloc(m * sizeof(float)); 
    for (int i = 0; i < m; i++) {
        max_scores[i] = FLT_MIN;
        for (int j = 0; j < 4; j++) {
            if (scoring_matrix[i][j] != NEG_INF && scoring_matrix[i][j] > max_scores[i]) {
                max_scores[i] = scoring_matrix[i][j];
            }
        }
    }

    // Calculate Z(i) values
    // This array holds the cumulative maximum possible score
    // from each position to the end of the scoring matrix
    // Starts from the second-to-last position and go backwards,
    // then adds the maximum score of the next position to the current Z value
    float* Z = (float*)malloc(m * sizeof(float));
    Z[m-1] = 0.0; // Last position has no remaining score
    for (int i = m - 2; i >= 0; i--) {
        Z[i] = Z[i + 1] + max_scores[i + 1];
    }

    // Calculate intermediate thresholds T_i
    // This array holds the intermediate thresholds at each position
    // Each threshold is the final threshold - the cumulative max score Z(i)
    float* T_intermediate = (float*)malloc(m * sizeof(float));
    for (int i = 0; i < m; i++) {
        T_intermediate[i] = t - Z[i];
    }

    // Create a Python object list to hold the search results
    PyObject *result_list = PyList_New(0);

    // Process each segment in sequence s, from sequence start to l-m
    for (int start = 0; start <= l - m; start++) {
        float total_score = 0.0;

        // Iterate over each position in the segment
        // Check position in segment
        for (int j = 0; j < m; j++) {
            char nucleotide = s[start + j];
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
                // For 'A' and 'a', score_index remains 0
                case 'A':
                case 'a':
                    score_index = 0;
                    break;
                default:
                    break;
            }

            // Add the score for the current nucleotide
            total_score += scoring_matrix[j][score_index];

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
    for (int i = 0; i < m; i++) {
        free(scoring_matrix[i]);
    }
    free(scoring_matrix);
    free(max_scores);
    free(Z);
    free(T_intermediate);

    return result_list;
}

// Define the methods that will be available in the Python module
static PyMethodDef motif_methods_la[] = {
    // Register the look_ahead function with Python, indicate it accepts arguments (METH_VARARGS)
    {"look_ahead", look_ahead, METH_VARARGS, "Find DNA motifs using lookahead scoring."},
    {NULL, NULL, 0, NULL} // Sentinel value, indicate the end of the array
};

// Define the Python module itself
static struct PyModuleDef look_ahead_module = {
    PyModuleDef_HEAD_INIT, // Macro to initialize the module definition
    "look_ahead",          // Name of the module
    NULL,                  // Module documentation (NULL)
    -1,                    // -1 as the module keeps state in global variables
    motif_methods_la       
};

// Initialization function for the module
// This function is called when the module is imported in Python
PyMODINIT_FUNC PyInit_look_ahead(void) {
    return PyModule_Create(&look_ahead_module); // Create the module and return a PyObject representing it
}
