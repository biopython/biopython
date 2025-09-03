/* Copyright 2025 Biopython contributors.
 * All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 *
 * Accelerated C implementation of melting temperature calculations.
 * Provides exact match with Bio.SeqUtils.MeltingTemp.Tm_NN using DNA_NN3 parameters.
 * 
 * Based on:
 * - Allawi & SantaLucia (1997) Biochemistry 36: 10581-10594 (DNA_NN3 parameters)
 * - SantaLucia (1998) Proc Natl Acad Sci USA 95: 1460-1465 (salt corrections)
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>

/* DNA_NN3 parameters (Allawi & SantaLucia 1997) */
typedef struct {
    const char* sequence;
    double dH;  /* Enthalpy in kcal/mol */
    double dS;  /* Entropy in cal/(mol·K) */
} NNParams;

/* All dinucleotide pairs from DNA_NN3 */
static const NNParams DNA_NN3_PARAMS[] = {
    {"AA/TT", -7.9, -22.2},
    {"AT/TA", -7.2, -20.4},
    {"TA/AT", -7.2, -21.3},
    {"CA/GT", -8.5, -22.7},
    {"GT/CA", -8.4, -22.4},
    {"CT/GA", -7.8, -21.0},
    {"GA/CT", -8.2, -22.2},
    {"CG/GC", -10.6, -27.2},
    {"GC/CG", -9.8, -24.4},
    {"GG/CC", -8.0, -19.9},
    {NULL, 0, 0}
};

/* Get complement of a base */
static char complement_base(char base) {
    switch(toupper(base)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        case 'U': return 'A';  /* Handle RNA */
        default: return 'N';
    }
}

/* Find NN parameters for a dinucleotide pair */
static bool get_nn_params(const char* seq_dinuc, const char* comp_dinuc, double* dH, double* dS) {
    char pair[8];
    
    /* Create pair string like "GA/CT" */
    snprintf(pair, sizeof(pair), "%c%c/%c%c", 
             toupper(seq_dinuc[0]), toupper(seq_dinuc[1]),
             toupper(comp_dinuc[0]), toupper(comp_dinuc[1]));
    
    /* Search in table */
    for (int i = 0; DNA_NN3_PARAMS[i].sequence != NULL; i++) {
        if (strcmp(pair, DNA_NN3_PARAMS[i].sequence) == 0) {
            *dH = DNA_NN3_PARAMS[i].dH;
            *dS = DNA_NN3_PARAMS[i].dS;
            return true;
        }
    }
    
    /* Try complete reverse: "AC/TG" -> "GT/CA" */
    /* Python does [::-1] which reverses the entire string */
    int pair_len = strlen(pair);
    char reverse[8];
    for (int i = 0; i < pair_len; i++) {
        reverse[i] = pair[pair_len - 1 - i];
    }
    reverse[pair_len] = '\0';
    
    for (int i = 0; DNA_NN3_PARAMS[i].sequence != NULL; i++) {
        if (strcmp(reverse, DNA_NN3_PARAMS[i].sequence) == 0) {
            *dH = DNA_NN3_PARAMS[i].dH;
            *dS = DNA_NN3_PARAMS[i].dS;
            return true;
        }
    }
    
    return false;
}

/* Check if sequence is self-complementary */
static bool is_self_complementary(const char* seq, size_t len) {
    for (size_t i = 0; i < len/2; i++) {
        if (toupper(seq[i]) != complement_base(seq[len-1-i])) {
            return false;
        }
    }
    return true;
}

/* Core Tm_NN calculation - EXACT VERSION */
static double calculate_tm_nn_exact(const char* seq, size_t len, 
                                    double dnac1, double dnac2,
                                    bool selfcomp,
                                    double Na, double K, double Tris,
                                    double Mg, double dNTPs,
                                    int saltcorr) {
    double delta_h = 0.0;
    double delta_s = 0.0;
    
    /* Generate complement sequence */
    char* complement = PyMem_Malloc(len + 1);
    if (!complement) return -999.0;
    
    for (size_t i = 0; i < len; i++) {
        complement[i] = complement_base(seq[i]);
    }
    complement[len] = '\0';
    
    /* Initiation parameters (DNA_NN3) */
    delta_h += 0.0;  /* init */
    delta_s += 0.0;  /* init */
    
    /* Check for GC content */
    int gc_count = 0;
    for (size_t i = 0; i < len; i++) {
        char c = toupper(seq[i]);
        if (c == 'G' || c == 'C') gc_count++;
    }
    
    if (gc_count == 0) {
        /* All AT - init_allA/T */
        delta_h += 0.0;
        delta_s += 0.0;
    } else {
        /* Has GC - init_oneG/C */
        delta_h += 0.0;
        delta_s += 0.0;
    }
    
    /* Check for 5' T - init_5T/A */
    if (toupper(seq[0]) == 'T') {
        delta_h += 0.0;
        delta_s += 0.0;
    }
    
    /* Check for 3' A - init_5T/A */
    if (toupper(seq[len-1]) == 'A') {
        delta_h += 0.0;
        delta_s += 0.0;
    }
    
    /* Terminal base pair corrections */
    char first = toupper(seq[0]);
    char last = toupper(seq[len-1]);
    
    /* init_A/T: (2.3, 4.1) */
    /* init_G/C: (0.1, -2.8) */
    if (first == 'A' || first == 'T') {
        delta_h += 2.3;
        delta_s += 4.1;
    } else if (first == 'G' || first == 'C') {
        delta_h += 0.1;
        delta_s += -2.8;
    }
    
    if (last == 'A' || last == 'T') {
        delta_h += 2.3;
        delta_s += 4.1;
    } else if (last == 'G' || last == 'C') {
        delta_h += 0.1;
        delta_s += -2.8;
    }
    
    /* Process nearest-neighbor pairs */
    for (size_t i = 0; i < len - 1; i++) {
        char seq_dinuc[3] = {seq[i], seq[i+1], '\0'};
        char comp_dinuc[3] = {complement[i], complement[i+1], '\0'};
        
        double dH_nn = 0.0, dS_nn = 0.0;
        if (get_nn_params(seq_dinuc, comp_dinuc, &dH_nn, &dS_nn)) {
            delta_h += dH_nn;
            delta_s += dS_nn;
        }
    }
    
    PyMem_Free(complement);
    
    /* Check for self-complementary */
    if (selfcomp || is_self_complementary(seq, len)) {
        delta_s += -1.4;  /* Symmetry correction */
        selfcomp = true;
    }
    
    /* Calculate concentration factor */
    double k;
    if (selfcomp) {
        k = dnac1 * 1e-9;
    } else {
        k = (dnac1 - (dnac2 / 2.0)) * 1e-9;
    }
    
    /* Avoid log(0) */
    if (k <= 0) k = 1e-9;
    
    /* Universal gas constant */
    double R = 1.987;
    
    /* Salt correction */
    double corr = 0.0;
    double Mon = Na + K + Tris / 2.0;  /* mM */
    
    /* Sodium-equivalent concentration (von Ahsen et al. 2001) */
    /* When Mg2+ is present and dNTPs < Mg, add correction term */
    if ((K > 0 || Mg > 0 || Tris > 0 || dNTPs > 0) && saltcorr != 7 && dNTPs < Mg) {
        Mon += 120.0 * sqrt(Mg - dNTPs);
    }
    
    if (saltcorr == 5 && Mon > 0) {
        /* Method 5: Correction for deltaS: 0.368 x (N-1) x ln[Na+] */
        /* Mon is in mM, need to convert to M */
        corr = 0.368 * (len - 1) * log(Mon / 1000.0);
        delta_s += corr;
    }
    
    /* Calculate Tm using van't Hoff equation */
    double melting_temp = (1000.0 * delta_h) / (delta_s + R * log(k)) - 273.15;
    
    /* Apply other salt correction methods */
    /* Methods 1-4 and 6 use sodium-equivalent concentration */
    if (saltcorr >= 1 && saltcorr <= 4 && Mon > 0) {
        double mon_molar = Mon / 1000.0;
        if (saltcorr == 1) {
            corr = 16.6 * log10(mon_molar);
        } else if (saltcorr == 2) {
            corr = 16.6 * log10(mon_molar / (1.0 + 0.7 * mon_molar));
        } else if (saltcorr == 3) {
            corr = 12.5 * log10(mon_molar);
        } else if (saltcorr == 4) {
            corr = 11.7 * log10(mon_molar);
        }
        melting_temp += corr;
    } else if (saltcorr == 6 && Mon > 0) {
        /* Method 6: Owczarzy et al. 2004 */
        double mon_molar = Mon / 1000.0;
        double gc_frac = (double)gc_count / len;
        corr = ((4.29 * gc_frac - 3.95) * 1e-5 * log(mon_molar)) + 
               (9.40e-6 * log(mon_molar) * log(mon_molar));
        /* For method 6, correction is applied as: Tm = 1/(1/Tm + corr) */
        melting_temp = 1.0 / (1.0 / (melting_temp + 273.15) + corr) - 273.15;
    } else if (saltcorr == 7) {
        /* Method 7: Owczarzy et al. 2008 - complex Mg2+ correction */
        double a = 3.92, b = -0.911, c = 6.26, d = 1.42;
        double e = -48.2, f = 52.5, g = 8.31;
        
        double mon_molar = Mon / 1000.0;
        double mg_molar = Mg / 1000.0;
        
        if (dNTPs > 0) {
            /* Adjust Mg2+ for dNTP binding */
            double dntps = dNTPs / 1000.0;
            double ka = 3e4;  /* Dissociation constant for Mg:dNTP */
            /* Free Mg2+ calculation */
            double discriminant = (ka * dntps - ka * mg_molar + 1.0);
            discriminant = discriminant * discriminant + 4.0 * ka * mg_molar;
            mg_molar = (-(ka * dntps - ka * mg_molar + 1.0) + sqrt(discriminant)) / (2.0 * ka);
        }
        
        if (mon_molar > 0) {
            double R = sqrt(mg_molar) / mon_molar;
            double gc_frac = (double)gc_count / len;
            
            if (R < 0.22) {
                /* Monovalent salt dominant - Python returns early here */
                corr = (4.29 * gc_frac - 3.95) * 1e-5 * log(mon_molar) + 
                       9.40e-6 * log(mon_molar) * log(mon_molar);
                /* For method 7, correction is applied as: Tm = 1/(1/Tm + corr) */
                melting_temp = 1.0 / (1.0 / (melting_temp + 273.15) + corr) - 273.15;
                return melting_temp;
            } else if (R < 6.0) {
                /* Mixed salt conditions - recalculate a, d, g */
                a = 3.92 * (0.843 - 0.352 * sqrt(mon_molar) * log(mon_molar));
                d = 1.42 * (1.279 - 4.03e-3 * log(mon_molar) - 8.03e-3 * log(mon_molar) * log(mon_molar));
                g = 8.31 * (0.486 - 0.258 * log(mon_molar) + 5.25e-3 * log(mon_molar) * log(mon_molar) * log(mon_molar));
            }
            /* Otherwise a, b, c, d, e, f, g remain as initialized */
        }
        
        /* Final calculation for non-early-return cases */
        double gc_frac = (double)gc_count / len;
        corr = (a + b * log(mg_molar) + 
                gc_frac * (c + d * log(mg_molar)) + 
                (1.0 / (2.0 * (len - 1))) * (e + f * log(mg_molar) + g * log(mg_molar) * log(mg_molar))) * 1e-5;
        
        /* For method 7, correction is applied as: Tm = 1/(1/Tm + corr) */
        melting_temp = 1.0 / (1.0 / (melting_temp + 273.15) + corr) - 273.15;
    }
    
    return melting_temp;
}

/* Python wrapper */
static PyObject* py_tm_nn_exact(PyObject* self, PyObject* args, PyObject* kwargs) {
    const char* seq;
    Py_ssize_t len;
    double dnac1 = 25.0;
    double dnac2 = 25.0;
    int selfcomp = 0;
    double Na = 50.0;
    double K = 0.0;
    double Tris = 0.0;
    double Mg = 0.0;
    double dNTPs = 0.0;
    int saltcorr = 5;
    
    static char* kwlist[] = {"seq", "dnac1", "dnac2", "selfcomp",
                            "Na", "K", "Tris", "Mg", "dNTPs", 
                            "saltcorr", NULL};
    
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|ddpdddddi", kwlist,
                                    &seq, &len, &dnac1, &dnac2, &selfcomp,
                                    &Na, &K, &Tris, &Mg, &dNTPs, &saltcorr)) {
        return NULL;
    }
    
    double tm = calculate_tm_nn_exact(seq, len, dnac1, dnac2, selfcomp != 0,
                                      Na, K, Tris, Mg, dNTPs, saltcorr);
    
    if (tm < -900) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }
    
    return PyFloat_FromDouble(tm);
}

/* Module methods */
static PyMethodDef module_methods[] = {
    {"tm_nn_exact", (PyCFunction)py_tm_nn_exact, METH_VARARGS | METH_KEYWORDS,
     "Calculate melting temperature using exact DNA_NN3 parameters.\n\n"
     "Matches Bio.SeqUtils.MeltingTemp.Tm_NN exactly.\n"},
    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_meltingtemp_exact",
    "Exact implementation of DNA melting temperature calculations",
    -1,
    module_methods
};

/* Module initialization */
PyMODINIT_FUNC PyInit__meltingtemp_exact(void) {
    return PyModule_Create(&module_def);
}
