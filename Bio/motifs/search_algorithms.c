/*Copyright (C) 2025, Joaquim Calvera i Madaula
*
* This file is part of the Biopython distribution and governed by your
* choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
* Please see the LICENSE file that should have been included as part of this
* package.
*/
#include "search_algorithms.h"
#include "darray.h"
#include <stdlib.h>
#include <time.h>

AlgorithmType parse_algorithm(const char* name){
    if (strcmp(name, "lookahead") == 0) return LOOKAHEAD;
    if (strcmp(name, "permuted") == 0) return PERMUTED_LOOKAHEAD;
    if (strcmp(name, "superalphabet") == 0) return SUPERALPHABET;
    return UNKNOWN_ALGORITHM;
}

/*This algorithm is based on the idea that, at each step of the score computation,
* there is a maximum possible score that can still be added from the remaining positions.
* Therefore, if the current score falls below (Threshold - S),
* where S is the sum of the maximum values in the scoring matrix for the remaining positions,
* we can safely discard the current subsequence.
* 
* Preprocessing:
*   - For each motif position (i) (row), store in min_scores[i] the minimum score a 
*     subsequence must have reached up to that position to still be able to reach the 
*     threshold by the end:
*       min_scores[i] = threshold - sum of the maximum values in the matrix from position i onward.
*
* Lookahead:
*   - For each subsequence sequence:
*       - Initialize score = 0
*       - For each motif (j):
*           - Add the corresponding matrix score: score += matrix[j * 4 + col] (col is the base at sequence[j+i])
*           - If the partial score is less than min_scores[j], break early (prune)
*       - If the full subsequence score ≥ threshold, store the match.
*
* Arguments:
*   - sequence:   Genomic sequence to be analyzed.
*   - s:          Length of the sequence.
*   - matrix:     PSSM matrix.
*   - m:          Motif length.
*   - threshold:  Minimum score for a subsequence to be considered a match.
*   - scores:     DArray to store matching positions and their scores.
*   - base_lookup: int[256] array to get the value of a base.
*/

int search_lookahead(const char sequence[], Py_ssize_t s, double* matrix, Py_ssize_t m, 
                    float threshold, DArray* scores, const int* base_lookup){
   
    Py_ssize_t i, j;
    float max_score;

    float* min_scores = malloc(m * sizeof(float)); 
    if(!min_scores){
        PyErr_NoMemory();
        return -1;
    }

    //Fill min_scores with (threshold - S)
    min_scores[m-1] = threshold;
    for(i = m-2; i >= 0; i--)
    {
        max_score = -1.0e30;
        for(j = 0; j < 4; j++)
        {
            if(matrix[(i+1)*4+j] > max_score)
                max_score = matrix[(i+1)*4+j];
        }
        min_scores[i] = min_scores[i+1] - max_score;
    }

    //lookahead search
    for (i = 0; i <= s - m; i++) {
        float score = 0.0;
        Py_ssize_t col = 0;
        for (j = 0; j < m && col >= 0; j++) {
            char base = sequence[i + j];
            col = base_lookup[base];
            if (col >= 0){
                score += matrix[j * 4 + col];
                if (score < min_scores[j])
                    col = -1;
            }
        }
        if (col > -1) {
            if (insert_in_darray(scores, i, score) < 0) {
                free(min_scores);
                return -1;
            }
        }
    }
    free(min_scores);
    return 0;
}


/* Comparator function for qsort. 
*  Compares two PermEntry elements based on their priority value (descending order).
*/
int compare_permutation(const void* a, const void* b) {
    float diff;
    diff = ((PermEntry*)b)->priority - ((PermEntry*)a)->priority;
    return (diff > 0) - (diff < 0);
}


/*This function computes an optimal permutation of the PSSM matrix rows to enable the most 
* discriminatory ones to be evaluated first.
* A motif position is considered more discriminatory when the maximum attainable score is
* larger than the expected score (computed based on background frequencies).
* 
* For every motif position (row) i:
*     - maxval = maximum score across A, C, G, T
*     - expval = expected score:  sum_j(matrix[i][j] * bkg[j])
*     - priority = maxval - expval
* 
* Rows are sorted in descending order of priority. 
* The resulting order is stored in the perm array.
*/
void compute_permutation(double* matrix, Py_ssize_t m, const float bkg[4], int* perm) {

    Py_ssize_t i, j;

    PermEntry* entries = malloc(m * sizeof(PermEntry));
    if (!entries) {
        PyErr_NoMemory();
        return;
    }

    for (i = 0; i < m; i++) {
        float maxval = -1e30;
        float expval = 0.0;

        for (j = 0; j < 4; j++) {
            float val = matrix[i * 4 + j];
            if (val > maxval)
                maxval = val;
            expval += val * bkg[j];
        }
        entries[i].index = i;
        entries[i].priority = maxval - expval;
    }
    qsort(entries, m, sizeof(PermEntry), compare_permutation);
    for (i = 0; i < m; i++) 
        perm[i] = entries[i].index;

    free(entries);
}

/*Performs motif search using the permuted lookahead technique.
*
* This is an optimized version of the sequential lookahead algorithm,
* where matrix rows (positions in the motif) are evaluated in a permuted
* order to maximize the chance of early rejection of low-scoring subsequences.
* 
* Preprocessing:
*
*   - Initialize 'perm' with the best possible order to scan the motif postions (matrix rows),
*      this is done using the 'compute_permutation' function.
*   - For each motif position (i) (row), store in min_scores[i] the minimum score a 
*     subsequence must have reached up to that position to still be able to reach the 
*     threshold by the end:
*       min_scores[i] = threshold - sum of the maximum values in the matrix from position i onward.
*
* Lookahead:
*   - For each subsequence sequence:
*       - Initialize score = 0
*       - For each motif perm[j]:
*           - Add the corresponding matrix score: score += matrix[perm[j] * 4 + col]
*           - If the partial score is less than min_scores[j], break early (prune)
*       - If the full subsequence score ≥ threshold, store the match.
*
* Arguments:
*   - sequence:   Genomic sequence to be analyzed.
*   - s:          Length of the sequence.
*   - matrix:     Original PWM matrix (using the standard alphabet).
*   - m:          Motif length (based on the original alphabet).
*   - threshold:  Minimum score for a subsequence to be considered a match.
*   - scores:     DArray to store matching positions and their scores.
*   - base_lookup: int[256] array to get the value of a base.
*/
int search_permuted_lookahead(const char sequence[], Py_ssize_t s,  double* matrix, Py_ssize_t m, 
                            float threshold, DArray* scores, const int* base_lookup){

    Py_ssize_t i, j;
    float bkg[4] = {0.25,0.25,0.25,0.25};  
    
    int* perm = malloc(m * sizeof(int)); //stores the order of evaluation of the subsequence positions
    if (!perm) {
        PyErr_NoMemory();
        return -1;
    }

    compute_permutation(matrix, m, bkg, perm);

    float *min_scores = malloc(m * sizeof(float)); 
    if (!min_scores) {
        PyErr_NoMemory();
        free(perm);
        return -1;
    }

    // Fill min_scores with (threshold - S) in the permutated order.
    min_scores[m-1] = threshold;
    for (i = m - 2; i >= 0; i--) {
        float max_score = -1e30;
        for (j = 0; j < 4; j++) {
            float val = matrix[perm[i+1] * 4 + j];
            if (val > max_score) 
                max_score = val;
        }
        min_scores[i] = min_scores[i+1] - max_score;
    }

    // Search lookahead in the permutated order
    for (i = 0; i <= s - m; i++) {
        float score = 0.0;
        Py_ssize_t col = 0;
        for (j = 0; j < m && col > -1; j++) {
            char base = sequence[i + perm[j]];
            col = base_lookup[base];
            if (col > -1){
                score += matrix[perm[j] * 4 + col];
                if (score < min_scores[j])
                    col = -1;
            }
        }
        if (col > -1) {
            if (insert_in_darray(scores, i, score) < 0) {
                free(min_scores);
                free(perm);
                return -1;
            }
        }
    }

    free(min_scores);
    free(perm);
    return 0;
}


/*Precomputes the "super-alphabet" scoring matrix mhat. 
* Each block of q rows in the original PWM is combined into one.
* Args:
*     - mat:   Original pssm using the original alphabet [A,C,G,T]
*     - m:     Length of the motiv using the original alphabet 
*     - q:     Length of the symbols of the new sup-eralphabet
*     - mhat:  New pssm matrix
*/
int compute_mhat(const double* mat, Py_ssize_t m, Py_ssize_t q, double* mhat) {
    Py_ssize_t num_blocks = (m + q - 1) / q;           //length of the motif using the new alphabet ⌈m / q⌉
    Py_ssize_t sigma_q = 1;     //length of the new superalphabet 4^q
    for (Py_ssize_t i = 0; i < q; i++) {
        sigma_q *= 4;
    }
    int *digits = malloc(q*sizeof(int));
    if(!digits){
        PyErr_NoMemory();
        return -1;
    }

    for (Py_ssize_t j = 0; j < num_blocks; j++) {
        Py_ssize_t remaining = m - j * q;                          //if m / q is not divisible the last block 
        Py_ssize_t effective_q = remaining < q ? remaining : q;    //will be smaller than q.

        //for every symbol in the super-alphabet AAA-->TTT
        for (Py_ssize_t sym = 0; sym < sigma_q; sym++) {
            //Decode sym as a base-4 number → digits[0..q-1]
            Py_ssize_t tmp = sym;
            for (Py_ssize_t h = q - 1; h >= 0; h--) {
                digits[h] = tmp % 4;
                tmp /= 4;
            }

            //Compute the sum of scores.
            //If the last block is smaller than q, we use only the valid positions.
            float sum = 0.0f;
            for (Py_ssize_t h = 0; h < effective_q; h++) {
                Py_ssize_t row = j * q + h;
                Py_ssize_t base = digits[h];
                sum += mat[row * 4 + base];  // mat[row][base]
            }

            mhat[j * sigma_q + sym] = sum;
        }
    }
    free(digits);
}

/* This algorithm is an improvement on the brute force approach by operating on a "super-alphabet"
* instead of the original nucleotide alphabet. Each q-tuple from the original alphabet is a new symbol of the superalphabet.
* The original matrix is preprocessed into an equivalent scoring matrix mhat for super-alphabet symbols. 
* This reduces the number of steps required during search, allowing faster matching.
*   
*   To decode alphabet symbols in the sequence:
*   Given a DNA subsequence (qtuple) length "len", we calculate its corresponding index
*   in the sorted array (AAA ---> TTT) of all possible sequences of that length using
*   base-4 encoding:
*       If q = 3 and 
*       seq = "CAT": index = 1*4^2 + 0*4^1 + 3*4^0 = 16 + 0 + 3 = 19
*       seq = "CA" (padded to "CAA"): index = 1*4^2 + 0*4^1 + 0*4^0 = 16
*       
*
* Arguments:
*   - sequence:   Genomic sequence to be analyzed.
*   - s:          Length of the sequence.
*   - matrix:     Original PWM matrix (using the standard alphabet).
*   - m:          Motif length (based on the original alphabet).
*   - threshold:  Minimum score for a subsequence to be considered a match.
*   - scores:     DArray to store matching positions and their scores.
*   - q:          Length of the q-grams (super-alphabet symbols).
*   - base_lookup: int[256] array to get the value of a base.
*/
int search_superalphabet(const char sequence[], Py_ssize_t s, double* matrix, Py_ssize_t m, 
                        float threshold, DArray* scores, Py_ssize_t q, const int* base_lookup){

    int result = -1;
    Py_ssize_t i, j;
    Py_ssize_t sa_size = 1;  //number of q-tuples in the super-alphabet (4^q) 
    for (i = 0; i < q; i++) {
        sa_size *= 4;
    }
    Py_ssize_t sa_motif_length = (m + q -1)/q;  //⌈m / q⌉ length of the motif using the new super-alphabet

    //Init the super-alphabet scoring matrix.
    double* sa_matrix = malloc(sa_motif_length * sa_size * sizeof(double));
    if (!sa_matrix) {
        PyErr_NoMemory();
        return -1;
    }

    //Compute the super-alphabet scoring matrix.
    result = compute_mhat(matrix, m, q, sa_matrix);
    if(result < 0){
        free(sa_matrix);
        return -1;
    }

    for (i = 0; i <= s - m; i++) {
        float score = 0.0;
        int valid = 1;

        for (j = 0; j < sa_motif_length && valid; j++) {
            Py_ssize_t offset = i + j * q;
            Py_ssize_t remaining = m - j * q;
            Py_ssize_t len = remaining < q ? remaining : q;

            //qtuple to index
            Py_ssize_t idx = 0;
            for (Py_ssize_t x = 0; x < q; x++) {
                Py_ssize_t b;
                if (x < len) {
                    b = base_lookup[sequence[offset+x]];
                    if (b == -1){  
                        idx = -1; //if the letter is not valid -> skip subsequence.
                        break;
                    }
                } else {
                    b = 0;   //if the sequence is shorter than q we pad with "A"
                }
                idx = idx * 4 + b;
            }
            if (idx < 0) 
                valid = 0;
            else
                score += sa_matrix[j * sa_size + idx];
        }
        //if the sequen
        if (valid && score >= threshold) {
            if (insert_in_darray(scores, i, score) < 0) {
                free(sa_matrix);
                return -1;
            }
        }
    }

    free(sa_matrix);
    return 0;
}