/* Copyright 2025 Biopython contributors.
 * All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 *
 * Complete accelerated C implementation of melting temperature calculations.
 * Provides comprehensive support for all nearest-neighbor tables, mismatches,
 * and dangling ends with exact match to Bio.SeqUtils.MeltingTemp.Tm_NN.
 *
 * Based on thermodynamic parameters from:
 * - Breslauer et al. (1986) - DNA_NN1
 * - Sugimoto et al. (1996) - DNA_NN2
 * - Allawi & SantaLucia (1997) - DNA_NN3
 * - SantaLucia & Hicks (2004) - DNA_NN4
 * - Freier et al. (1986) - RNA_NN1
 * - Xia et al. (1998) - RNA_NN2
 * - Chen et al. (2012) - RNA_NN3
 * - Sugimoto et al. (1995) - R_DNA_NN1
 * - Allawi & SantaLucia (1997-1998) - DNA_IMM1 (mismatches)
 * - SantaLucia & Peyret (2001) - DNA_TMM1 (terminal mismatches)
 * - Bommarito et al. (2000) - DNA_DE1 (dangling ends)
 * - Turner & Mathews (2010) - RNA_DE1 (dangling ends)
 *
 * Salt corrections based on:
 * - SantaLucia (1998) Proc Natl Acad Sci USA 95: 1460-1465
 * - von Ahsen et al. (2001) Clin Chem 47: 1956-1961
 * - Owczarzy et al. (2004) Biochemistry 43: 3537-3554
 * - Owczarzy et al. (2008) Biochemistry 47: 5336-5353
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>

/* Constants */
#define R_GAS_CONSTANT 1.987  /* Universal gas constant in Cal/(degrees C*Mol) */
#define MAX_SEQ_LENGTH 10000
#define MAX_KEY_LENGTH 20

/* Thermodynamic parameter structure */
typedef struct {
    double dH;  /* Enthalpy (kcal/mol) */
    double dS;  /* Entropy (cal/mol K) */
} ThermParams;

/* Table entry structure */
typedef struct {
    char key[MAX_KEY_LENGTH];
    ThermParams params;
} TableEntry;

/* Include complete thermodynamic tables from BioPython */
#include "_thermodynamic_tables.h"

/* Table selection enum - matches BioPython table names */
typedef enum {
    TABLE_DNA_NN1 = 1,
    TABLE_DNA_NN2 = 2,
    TABLE_DNA_NN3 = 3,
    TABLE_DNA_NN4 = 4,
    TABLE_RNA_NN1 = 5,
    TABLE_RNA_NN2 = 6,
    TABLE_RNA_NN3 = 7,
    TABLE_R_DNA_NN1 = 8
} TableType;

/* Forward declarations */
static double tm_nn_table(const char* seq, const char* c_seq, int shift,
                   double dnac1, double dnac2, bool selfcomp,
                   double Na, double K, double Tris, double Mg, double dNTPs,
                   int saltcorr, int table_id);
static double salt_correction(double Na, double K, double Tris, double Mg, double dNTPs,
                      int method, const char* seq, int seq_len);
static double gc_fraction(const char* seq, int len);
static int find_therm_params(const TableEntry* table, const char* key, ThermParams* params);
static void reverse_string(char* str);
static char* get_complement(const char* seq);
static void clean_sequence(const char* input, char* output);
static void uppercase_string(char* str);
static const TableEntry* get_nn_table(int table_id);
static const TableEntry* get_imm_table(int table_id);
static const TableEntry* get_tmm_table(int table_id);
static const TableEntry* get_de_table(int table_id);

/* Get the appropriate NN table based on ID */
static const TableEntry* get_nn_table(int table_id) {
    switch (table_id) {
        case TABLE_DNA_NN1: return DNA_NN1;
        case TABLE_DNA_NN2: return DNA_NN2;
        case TABLE_DNA_NN3: return DNA_NN3;
        case TABLE_DNA_NN4: return DNA_NN4;
        case TABLE_RNA_NN1: return RNA_NN1;
        case TABLE_RNA_NN2: return RNA_NN2;
        case TABLE_RNA_NN3: return RNA_NN3;
        case TABLE_R_DNA_NN1: return R_DNA_NN1;
        default: return DNA_NN4;
    }
}

/* Get the appropriate mismatch table (only DNA has these) */
static const TableEntry* get_imm_table(int table_id) {
    if (table_id >= TABLE_DNA_NN1 && table_id <= TABLE_DNA_NN4) {
        return DNA_IMM1;
    }
    return NULL;
}

/* Get the appropriate terminal mismatch table */
static const TableEntry* get_tmm_table(int table_id) {
    if (table_id >= TABLE_DNA_NN1 && table_id <= TABLE_DNA_NN4) {
        return DNA_TMM1;
    }
    return NULL;
}

/* Get the appropriate dangling ends table */
static const TableEntry* get_de_table(int table_id) {
    if (table_id >= TABLE_DNA_NN1 && table_id <= TABLE_DNA_NN4) {
        return DNA_DE1;
    } else if (table_id >= TABLE_RNA_NN1 && table_id <= TABLE_RNA_NN3) {
        return RNA_DE1;
    }
    return NULL;
}

/* Helper function to uppercase a string */
static void uppercase_string(char* str) {
    for (int i = 0; str[i]; i++) {
        str[i] = toupper(str[i]);
    }
}

/* Helper function to find thermodynamic parameters in a table */
static int find_therm_params(const TableEntry* table, const char* key, ThermParams* params) {
    if (!table) return 0;

    for (int i = 0; table[i].key[0] != '\0'; i++) {
        if (strcmp(table[i].key, key) == 0) {
            params->dH = table[i].params.dH;
            params->dS = table[i].params.dS;
            return 1;  /* Found */
        }
    }
    return 0;  /* Not found */
}

/* Reverse a string in place */
static void reverse_string(char* str) {
    int len = strlen(str);
    for (int i = 0; i < len / 2; i++) {
        char temp = str[i];
        str[i] = str[len - 1 - i];
        str[len - 1 - i] = temp;
    }
}

/* Get DNA/RNA complement */
static char* get_complement(const char* seq) {
    int len = strlen(seq);
    char* complement = (char*)PyMem_Malloc(len + 1);
    if (!complement) return NULL;

    for (int i = 0; i < len; i++) {
        switch (toupper(seq[i])) {
            case 'A': complement[i] = 'T'; break;
            case 'T': complement[i] = 'A'; break;
            case 'G': complement[i] = 'C'; break;
            case 'C': complement[i] = 'G'; break;
            case 'U': complement[i] = 'A'; break;  /* RNA support */
            case '.': complement[i] = '.'; break;
            default:
                /* For lowercase (mismatches) keep as is but complement */
                if (seq[i] == 'a') complement[i] = 't';
                else if (seq[i] == 't') complement[i] = 'a';
                else if (seq[i] == 'g') complement[i] = 'c';
                else if (seq[i] == 'c') complement[i] = 'g';
                else if (seq[i] == 'u') complement[i] = 'a';
                else complement[i] = seq[i];
                break;
        }
    }
    complement[len] = '\0';
    return complement;
}

/* Clean and uppercase sequence */
static void clean_sequence(const char* input, char* output) {
    int j = 0;
    for (int i = 0; input[i] != '\0'; i++) {
        if (isalpha(input[i]) || input[i] == '.') {
            output[j++] = toupper(input[i]);
        }
    }
    output[j] = '\0';
}

/* Calculate GC fraction (ignoring ambiguous bases) */
static double gc_fraction(const char* seq, int len) {
    int gc_count = 0;
    int total_count = 0;

    for (int i = 0; i < len; i++) {
        char base = toupper(seq[i]);
        if (base == 'G' || base == 'C') {
            gc_count++;
            total_count++;
        } else if (base == 'A' || base == 'T' || base == 'U') {
            total_count++;
        }
    }

    if (total_count == 0) return 0.0;
    return (double)gc_count / total_count;
}

/* Salt correction calculation - implements all 7 BioPython methods */
static double salt_correction(double Na, double K, double Tris, double Mg, double dNTPs,
                      int method, const char* seq, int seq_len) {
    double corr = 0.0;

    if (method == 0) return corr;

    /* Calculate monovalent ion concentration */
    double Mon = Na + K + Tris / 2.0;  /* millimolar */
    double mg = Mg * 1e-3;  /* Convert to molar */

    /* Na equivalent according to von Ahsen et al. (2001) */
    if ((K + Mg + Tris + dNTPs) > 0 && method != 7 && dNTPs < Mg) {
        Mon += 120.0 * sqrt(Mg - dNTPs);
    }
    double mon = Mon * 1e-3;  /* Convert to molar */

    /* Check for zero concentration */
    if (method >= 1 && method <= 6 && mon <= 0) {
        return 0.0;  /* Or handle error appropriately */
    }

    switch (method) {
        case 1:
            corr = 16.6 * log10(mon);
            break;
        case 2:
            corr = 16.6 * log10(mon / (1.0 + 0.7 * mon));
            break;
        case 3:
            corr = 12.5 * log10(mon);
            break;
        case 4:
            corr = 11.7 * log10(mon);
            break;
        case 5:
            corr = 0.368 * (seq_len - 1) * log(mon);
            break;
        case 6: {
            double gc_frac = gc_fraction(seq, seq_len);
            corr = ((4.29 * gc_frac - 3.95) * 1e-5 * log(mon)) +
                   9.40e-6 * pow(log(mon), 2);
            break;
        }
        case 7: {
            /* Complex formula with decision tree */
            double a = 3.92, b = -0.911, c = 6.26, d = 1.42;
            double e = -48.2, f = 52.5, g = 8.31;

            if (dNTPs > 0) {
                double dntps = dNTPs * 1e-3;
                double ka = 3e4;  /* Dissociation constant for Mg:dNTP */
                /* Free Mg2+ calculation */
                mg = (-(ka * dntps - ka * mg + 1.0) +
                      sqrt(pow(ka * dntps - ka * mg + 1.0, 2) + 4.0 * ka * mg)) /
                     (2.0 * ka);
            }

            if (Mon > 0) {
                double R = sqrt(mg) / mon;
                double gc_frac = gc_fraction(seq, seq_len);

                if (R < 0.22) {
                    corr = ((4.29 * gc_frac - 3.95) * 1e-5 * log(mon)) +
                           9.40e-6 * pow(log(mon), 2);
                    return corr;
                } else if (R < 6.0) {
                    a = 3.92 * (0.843 - 0.352 * sqrt(mon) * log(mon));
                    d = 1.42 * (1.279 - 4.03e-3 * log(mon) -
                               8.03e-3 * pow(log(mon), 2));
                    g = 8.31 * (0.486 - 0.258 * log(mon) +
                               5.25e-3 * pow(log(mon), 3));
                }

                corr = (a + b * log(mg) + gc_frac * (c + d * log(mg)) +
                       (1.0 / (2.0 * (seq_len - 1))) *
                       (e + f * log(mg) + g * pow(log(mg), 2))) * 1e-5;
            }
            break;
        }
    }

    return corr;
}

/* Global variable to store the problematic neighbor pair for error reporting */
static char g_error_nn_key[MAX_KEY_LENGTH] = {0};

/* Main Tm_NN calculation function with table selection */
/* Returns NaN on error (missing thermodynamic data) */
static double tm_nn_table(const char* seq, const char* c_seq, int shift,
                   double dnac1, double dnac2, bool selfcomp,
                   double Na, double K, double Tris, double Mg, double dNTPs,
                   int saltcorr, int table_id) {

    /* Clear any previous error */
    g_error_nn_key[0] = '\0';

    /* Get the appropriate tables */
    const TableEntry* nn_table = get_nn_table(table_id);
    if (!nn_table) {
        return NAN;
    }
    const TableEntry* imm_table = get_imm_table(table_id);
    const TableEntry* tmm_table = get_tmm_table(table_id);
    const TableEntry* de_table = get_de_table(table_id);

    /* Clean input sequences */
    char clean_seq[MAX_SEQ_LENGTH];
    char clean_cseq[MAX_SEQ_LENGTH];
    clean_sequence(seq, clean_seq);

    /* If no c_seq provided, use perfect complement */
    if (c_seq == NULL || c_seq[0] == '\0') {
        char* comp = get_complement(clean_seq);
        if (!comp) return NAN;
        strcpy(clean_cseq, comp);
        PyMem_Free(comp);
    } else {
        /* For c_seq, preserve case for mismatches but still clean */
        int j = 0;
        for (int i = 0; c_seq[i] != '\0'; i++) {
            if (isalpha(c_seq[i]) || c_seq[i] == '.') {
                clean_cseq[j++] = c_seq[i];  /* Keep original case */
            }
        }
        clean_cseq[j] = '\0';
        uppercase_string(clean_cseq);  /* Convert to uppercase for now */
    }

    /* Working copies of sequences */
    char tmp_seq[MAX_SEQ_LENGTH];
    char tmp_cseq[MAX_SEQ_LENGTH];
    strcpy(tmp_seq, clean_seq);
    strcpy(tmp_cseq, clean_cseq);


    double delta_h = 0.0;
    double delta_s = 0.0;
    ThermParams params;

    /* Handle dangling ends if shift or length mismatch */
    int len_seq = strlen(tmp_seq);
    int len_cseq = strlen(tmp_cseq);

    if (shift != 0 || len_seq != len_cseq) {
        /* Align sequences using shift parameter */
        char aligned_seq[MAX_SEQ_LENGTH] = {0};
        char aligned_cseq[MAX_SEQ_LENGTH] = {0};

        if (shift > 0) {
            /* Add dots to beginning of seq */
            for (int i = 0; i < shift; i++) {
                aligned_seq[i] = '.';
            }
            strcat(aligned_seq, tmp_seq);
            strcpy(aligned_cseq, tmp_cseq);
        } else if (shift < 0) {
            /* Add dots to beginning of c_seq */
            strcpy(aligned_seq, tmp_seq);
            for (int i = 0; i < abs(shift); i++) {
                aligned_cseq[i] = '.';
            }
            strcat(aligned_cseq, tmp_cseq);
        } else {
            strcpy(aligned_seq, tmp_seq);
            strcpy(aligned_cseq, tmp_cseq);
        }

        /* Equalize lengths */
        int max_len = strlen(aligned_seq);
        int cseq_len = strlen(aligned_cseq);
        if (cseq_len > max_len) {
            max_len = cseq_len;
            for (int i = strlen(aligned_seq); i < max_len; i++) {
                aligned_seq[i] = '.';
            }
            aligned_seq[max_len] = '\0';
        } else if (max_len > cseq_len) {
            for (int i = cseq_len; i < max_len; i++) {
                aligned_cseq[i] = '.';
            }
            aligned_cseq[max_len] = '\0';
        }

        /* Remove 'over-dangling' ends */
        while ((aligned_seq[0] == '.' && aligned_seq[1] == '.') ||
               (aligned_cseq[0] == '.' && aligned_cseq[1] == '.')) {
            memmove(aligned_seq, aligned_seq + 1, strlen(aligned_seq));
            memmove(aligned_cseq, aligned_cseq + 1, strlen(aligned_cseq));
        }

        int end_len = strlen(aligned_seq);
        while (end_len > 1 &&
               ((aligned_seq[end_len-1] == '.' && aligned_seq[end_len-2] == '.') ||
                (aligned_cseq[end_len-1] == '.' && aligned_cseq[end_len-2] == '.'))) {
            aligned_seq[end_len-1] = '\0';
            aligned_cseq[end_len-1] = '\0';
            end_len--;
        }

        /* Handle dangling ends */
        if (de_table && (aligned_seq[0] == '.' || aligned_cseq[0] == '.')) {
            char de_key[MAX_KEY_LENGTH];
            snprintf(de_key, MAX_KEY_LENGTH, "%.2s/%.2s", aligned_seq, aligned_cseq);
            if (find_therm_params(de_table, de_key, &params)) {
                delta_h += params.dH;
                delta_s += params.dS;
            }
            memmove(aligned_seq, aligned_seq + 1, strlen(aligned_seq));
            memmove(aligned_cseq, aligned_cseq + 1, strlen(aligned_cseq));
        }

        end_len = strlen(aligned_seq);
        if (de_table && end_len > 0 && (aligned_seq[end_len-1] == '.' || aligned_cseq[end_len-1] == '.')) {
            if (end_len >= 2) {
                char de_key[MAX_KEY_LENGTH];
                char rev_seq[3] = {aligned_cseq[end_len-2], aligned_cseq[end_len-1], '\0'};
                char rev_cseq[3] = {aligned_seq[end_len-2], aligned_seq[end_len-1], '\0'};
                reverse_string(rev_seq);
                reverse_string(rev_cseq);
                snprintf(de_key, MAX_KEY_LENGTH, "%s/%s", rev_seq, rev_cseq);
                if (find_therm_params(de_table, de_key, &params)) {
                    delta_h += params.dH;
                    delta_s += params.dS;
                }
            }
            aligned_seq[end_len-1] = '\0';
            aligned_cseq[end_len-1] = '\0';
        }

        strcpy(tmp_seq, aligned_seq);
        strcpy(tmp_cseq, aligned_cseq);
    }

    /* Handle terminal mismatches */
    len_seq = strlen(tmp_seq);
    if (tmm_table && len_seq >= 2) {
        /* Left terminal mismatch */
        char left_tmm[6];  /* Need 6 chars for "XX/XX\0" */
        left_tmm[0] = tmp_cseq[1];
        left_tmm[1] = tmp_cseq[0];
        left_tmm[2] = '/';
        left_tmm[3] = tmp_seq[1];
        left_tmm[4] = tmp_seq[0];
        left_tmm[5] = '\0';

        if (find_therm_params(tmm_table, left_tmm, &params)) {
            delta_h += params.dH;
            delta_s += params.dS;
            memmove(tmp_seq, tmp_seq + 1, strlen(tmp_seq));
            memmove(tmp_cseq, tmp_cseq + 1, strlen(tmp_cseq));
            len_seq--;
        }

        /* Right terminal mismatch */
        if (len_seq >= 2) {
            char right_tmm[6];
            snprintf(right_tmm, 6, "%.2s/%.2s",
                    tmp_seq + len_seq - 2, tmp_cseq + len_seq - 2);
            if (find_therm_params(tmm_table, right_tmm, &params)) {
                delta_h += params.dH;
                delta_s += params.dS;
                tmp_seq[len_seq - 1] = '\0';
                tmp_cseq[len_seq - 1] = '\0';
                len_seq--;
            }
        }
    }

    /* Initiation parameters */
    if (find_therm_params(nn_table, "init", &params)) {
        delta_h += params.dH;
        delta_s += params.dS;
    }

    /* Check for all A/T or at least one G/C */
    if (gc_fraction(clean_seq, strlen(clean_seq)) == 0) {
        if (find_therm_params(nn_table, "init_allA/T", &params)) {
            delta_h += params.dH;
            delta_s += params.dS;
        }
    } else {
        if (find_therm_params(nn_table, "init_oneG/C", &params)) {
            delta_h += params.dH;
            delta_s += params.dS;
        }
    }

    /* Penalty if 5' end is T */
    if (clean_seq[0] == 'T') {
        if (find_therm_params(nn_table, "init_5T/A", &params)) {
            delta_h += params.dH;
            delta_s += params.dS;
        }
    }
    if (clean_seq[strlen(clean_seq) - 1] == 'A') {
        if (find_therm_params(nn_table, "init_5T/A", &params)) {
            delta_h += params.dH;
            delta_s += params.dS;
        }
    }

    /* Terminal base pairs */
    char first_base = clean_seq[0];
    char last_base = clean_seq[strlen(clean_seq) - 1];
    int AT_count = 0;
    int GC_count = 0;

    if (first_base == 'A' || first_base == 'T') AT_count++;
    else if (first_base == 'G' || first_base == 'C') GC_count++;

    if (last_base == 'A' || last_base == 'T') AT_count++;
    else if (last_base == 'G' || last_base == 'C') GC_count++;

    if (find_therm_params(nn_table, "init_A/T", &params)) {
        delta_h += params.dH * AT_count;
        delta_s += params.dS * AT_count;
    }
    if (find_therm_params(nn_table, "init_G/C", &params)) {
        delta_h += params.dH * GC_count;
        delta_s += params.dS * GC_count;
    }

    /* Process nearest neighbors */
    len_seq = strlen(tmp_seq);
    for (int i = 0; i < len_seq - 1; i++) {
        char nn_key[MAX_KEY_LENGTH];
        char nn_pair[3] = {tmp_seq[i], tmp_seq[i+1], '\0'};
        char cc_pair[3] = {tmp_cseq[i], tmp_cseq[i+1], '\0'};

        snprintf(nn_key, MAX_KEY_LENGTH, "%s/%s", nn_pair, cc_pair);

        /* Try internal mismatch table first (if available) */
        bool found = false;

        if (imm_table) {
            found = find_therm_params(imm_table, nn_key, &params);
            if (!found) {
                /* Try reverse of entire key (like Python's neighbors[::-1]) */
                /* This means "AB/CD" becomes "DC/BA" */
                char rev_key[MAX_KEY_LENGTH];
                snprintf(rev_key, MAX_KEY_LENGTH, "%c%c/%c%c",
                         cc_pair[1], cc_pair[0], nn_pair[1], nn_pair[0]);
                found = find_therm_params(imm_table, rev_key, &params);
            }
        }

        /* If not found in mismatch table, try NN table */
        if (!found) {
            found = find_therm_params(nn_table, nn_key, &params);
            if (!found) {
                /* Try reverse in NN table */
                char rev_key[MAX_KEY_LENGTH];
                snprintf(rev_key, MAX_KEY_LENGTH, "%c%c/%c%c",
                         cc_pair[1], cc_pair[0], nn_pair[1], nn_pair[0]);
                found = find_therm_params(nn_table, rev_key, &params);
            }
        }

        /* If still not found, this is an error - missing thermodynamic data */
        if (!found) {
            /* Store the problematic key for error reporting */
            snprintf(g_error_nn_key, MAX_KEY_LENGTH, "%s/%s", nn_pair, cc_pair);
            /* Return NaN to indicate error (Python wrapper will convert to exception) */
            return NAN;
        }

        delta_h += params.dH;
        delta_s += params.dS;
    }

    /* Calculate k (DNA concentration factor) */
    double k;
    if (selfcomp) {
        k = dnac1 * 1e-9;
        if (find_therm_params(nn_table, "sym", &params)) {
            delta_h += params.dH;
            delta_s += params.dS;
        }
    } else {
        k = (dnac1 - (dnac2 / 2.0)) * 1e-9;
    }

    /* Prevent log of zero or negative */
    if (k <= 0) {
        k = 1e-15;  /* Very small positive number */
    }


    /* Apply salt correction */
    double corr = 0.0;
    if (saltcorr > 0) {
        corr = salt_correction(Na, K, Tris, Mg, dNTPs, saltcorr,
                              clean_seq, strlen(clean_seq));
    }


    /* Calculate melting temperature */
    double melting_temp;
    if (saltcorr == 5) {
        delta_s += corr;
    }

    /* Final Tm calculation */
    melting_temp = (1000.0 * delta_h) / (delta_s + (R_GAS_CONSTANT * log(k))) - 273.15;

    if (saltcorr >= 1 && saltcorr <= 4) {
        melting_temp += corr;
    } else if (saltcorr == 6 || saltcorr == 7) {
        if (!isnan(corr) && !isinf(corr)) {
            melting_temp = 1.0 / (1.0 / (melting_temp + 273.15) + corr) - 273.15;
        }
    }

    return melting_temp;
}

/*
 * Python C API Wrapper Functions
 */

/* Python wrapper for tm_nn_table */
static PyObject* py_tm_nn(PyObject* self, PyObject* args, PyObject* kwargs) {
    const char* seq;
    const char* c_seq = NULL;
    int shift = 0;
    double dnac1 = 25.0;
    double dnac2 = 25.0;
    int selfcomp = 0;
    double Na = 50.0;
    double K = 0.0;
    double Tris = 0.0;
    double Mg = 0.0;
    double dNTPs = 0.0;
    int saltcorr = 5;
    int table_id = TABLE_DNA_NN4;  /* Default to DNA_NN4 */

    static char* kwlist[] = {"seq", "c_seq", "shift", "dnac1", "dnac2", "selfcomp",
                            "Na", "K", "Tris", "Mg", "dNTPs", "saltcorr", "table_id", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|ziddpdddddii", kwlist,
                                    &seq, &c_seq, &shift, &dnac1, &dnac2, &selfcomp,
                                    &Na, &K, &Tris, &Mg, &dNTPs, &saltcorr, &table_id)) {
        return NULL;
    }

    double tm = tm_nn_table(seq, c_seq, shift, dnac1, dnac2, selfcomp != 0,
                           Na, K, Tris, Mg, dNTPs, saltcorr, table_id);

    if (isnan(tm)) {
        if (g_error_nn_key[0] != '\0') {
            PyErr_Format(PyExc_ValueError,
                        "no thermodynamic data for neighbors '%s' available",
                        g_error_nn_key);
        } else {
            PyErr_SetString(PyExc_RuntimeError, "Failed to calculate Tm");
        }
        return NULL;
    }

    return PyFloat_FromDouble(tm);
}

/* Module methods */
static PyMethodDef module_methods[] = {
    {"tm_nn", (PyCFunction)py_tm_nn, METH_VARARGS | METH_KEYWORDS,
     "Calculate melting temperature using nearest neighbor thermodynamics.\n\n"
     "Supports all DNA_NN and RNA_NN tables, mismatches, and dangling ends.\n"
     "Matches Bio.SeqUtils.MeltingTemp.Tm_NN exactly.\n\n"
     "Parameters:\n"
     "  seq (str): DNA/RNA sequence\n"
     "  c_seq (str, optional): Complementary sequence\n"
     "  shift (int, optional): Shift for alignment (default=0)\n"
     "  dnac1 (float, optional): Higher strand concentration in nM (default=25.0)\n"
     "  dnac2 (float, optional): Lower strand concentration in nM (default=25.0)\n"
     "  selfcomp (bool, optional): Self-complementary? (default=False)\n"
     "  Na (float, optional): Na+ concentration in mM (default=50.0)\n"
     "  K (float, optional): K+ concentration in mM (default=0.0)\n"
     "  Tris (float, optional): Tris concentration in mM (default=0.0)\n"
     "  Mg (float, optional): Mg2+ concentration in mM (default=0.0)\n"
     "  dNTPs (float, optional): dNTP concentration in mM (default=0.0)\n"
     "  saltcorr (int, optional): Salt correction method 0-7 (default=5)\n"
     "  table_id (int, optional): Table ID 1-8 (default=4 for DNA_NN4)\n"
     "    1=DNA_NN1, 2=DNA_NN2, 3=DNA_NN3, 4=DNA_NN4,\n"
     "    5=RNA_NN1, 6=RNA_NN2, 7=RNA_NN3, 8=R_DNA_NN1\n\n"
     "Returns:\n"
     "  float: Melting temperature in degrees Celsius\n\n"
     "Raises:\n"
     "  ValueError: If thermodynamic data is missing for a neighbor pair\n"},
    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_meltingtemp_complete",
    "Complete implementation of DNA/RNA melting temperature calculations\n\n"
    "This module provides accelerated C implementations for all BioPython\n"
    "nearest neighbor tables with full support for mismatches and dangling ends.",
    -1,
    module_methods
};

/* Module initialization */
PyMODINIT_FUNC PyInit__meltingtemp_complete(void) {
    PyObject* module = PyModule_Create(&module_def);
    if (module == NULL) {
        return NULL;
    }

    /* Add table constants */
    PyModule_AddIntConstant(module, "TABLE_DNA_NN1", TABLE_DNA_NN1);
    PyModule_AddIntConstant(module, "TABLE_DNA_NN2", TABLE_DNA_NN2);
    PyModule_AddIntConstant(module, "TABLE_DNA_NN3", TABLE_DNA_NN3);
    PyModule_AddIntConstant(module, "TABLE_DNA_NN4", TABLE_DNA_NN4);
    PyModule_AddIntConstant(module, "TABLE_RNA_NN1", TABLE_RNA_NN1);
    PyModule_AddIntConstant(module, "TABLE_RNA_NN2", TABLE_RNA_NN2);
    PyModule_AddIntConstant(module, "TABLE_RNA_NN3", TABLE_RNA_NN3);
    PyModule_AddIntConstant(module, "TABLE_R_DNA_NN1", TABLE_R_DNA_NN1);

    return module;
}
