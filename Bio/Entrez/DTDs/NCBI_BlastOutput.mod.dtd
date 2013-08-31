<!-- ============================================
     ::DATATOOL:: Generated from "blastxml.asn"
     ::DATATOOL:: by application DATATOOL version 2.0.0
     ::DATATOOL:: on 08/02/2010 23:05:14
     ============================================ -->

<!-- ============================================ -->
<!-- This section is mapped from module "NCBI-BlastOutput"
================================================= -->

<!--$Id: blastxml.asn 100080 2007-03-12 16:05:35Z kazimird $ -->


<!ELEMENT BlastOutput (
        BlastOutput_program, 
        BlastOutput_version, 
        BlastOutput_reference, 
        BlastOutput_db, 
        BlastOutput_query-ID, 
        BlastOutput_query-def, 
        BlastOutput_query-len, 
        BlastOutput_query-seq?, 
        BlastOutput_param, 
        BlastOutput_iterations, 
        BlastOutput_mbstat?)>

<!-- BLAST program: blastp, tblastx etc. -->
<!ELEMENT BlastOutput_program (#PCDATA)>

<!-- Program version  -->
<!ELEMENT BlastOutput_version (#PCDATA)>

<!-- Steven, David, Tom and others -->
<!ELEMENT BlastOutput_reference (#PCDATA)>

<!-- BLAST Database name -->
<!ELEMENT BlastOutput_db (#PCDATA)>

<!-- SeqId of query -->
<!ELEMENT BlastOutput_query-ID (#PCDATA)>

<!-- Definition line of query -->
<!ELEMENT BlastOutput_query-def (#PCDATA)>

<!-- length of query sequence -->
<!ELEMENT BlastOutput_query-len (%INTEGER;)>

<!-- query sequence itself -->
<!ELEMENT BlastOutput_query-seq (#PCDATA)>

<!-- search parameters -->
<!ELEMENT BlastOutput_param (Parameters)>

<!ELEMENT BlastOutput_iterations (Iteration*)>

<!-- Mega BLAST search statistics -->
<!ELEMENT BlastOutput_mbstat (Statistics)>


<!ELEMENT Iteration (
        Iteration_iter-num, 
        Iteration_query-ID?, 
        Iteration_query-def?, 
        Iteration_query-len?, 
        Iteration_hits?, 
        Iteration_stat?, 
        Iteration_message?)>

<!-- iteration number -->
<!ELEMENT Iteration_iter-num (%INTEGER;)>

<!-- SeqId of query -->
<!ELEMENT Iteration_query-ID (#PCDATA)>

<!-- Definition line of query -->
<!ELEMENT Iteration_query-def (#PCDATA)>

<!-- length of query sequence -->
<!ELEMENT Iteration_query-len (%INTEGER;)>

<!-- Hits one for every db sequence -->
<!ELEMENT Iteration_hits (Hit*)>

<!-- search statistics             -->
<!ELEMENT Iteration_stat (Statistics)>

<!-- Some (error?) information -->
<!ELEMENT Iteration_message (#PCDATA)>


<!ELEMENT Parameters (
        Parameters_matrix?, 
        Parameters_expect, 
        Parameters_include?, 
        Parameters_sc-match?, 
        Parameters_sc-mismatch?, 
        Parameters_gap-open, 
        Parameters_gap-extend, 
        Parameters_filter?, 
        Parameters_pattern?, 
        Parameters_entrez-query?)>

<!-- Matrix used (-M) -->
<!ELEMENT Parameters_matrix (#PCDATA)>

<!-- Expectation threshold (-e) -->
<!ELEMENT Parameters_expect (%REAL;)>

<!-- Inclusion threshold (-h) -->
<!ELEMENT Parameters_include (%REAL;)>

<!-- match score for NT (-r) -->
<!ELEMENT Parameters_sc-match (%INTEGER;)>

<!-- mismatch score for NT (-q) -->
<!ELEMENT Parameters_sc-mismatch (%INTEGER;)>

<!-- Gap opening cost (-G) -->
<!ELEMENT Parameters_gap-open (%INTEGER;)>

<!-- Gap extension cost (-E) -->
<!ELEMENT Parameters_gap-extend (%INTEGER;)>

<!-- Filtering options (-F) -->
<!ELEMENT Parameters_filter (#PCDATA)>

<!-- PHI-BLAST pattern -->
<!ELEMENT Parameters_pattern (#PCDATA)>

<!-- Limit of request to Entrez query -->
<!ELEMENT Parameters_entrez-query (#PCDATA)>


<!ELEMENT Statistics (
        Statistics_db-num, 
        Statistics_db-len, 
        Statistics_hsp-len, 
        Statistics_eff-space, 
        Statistics_kappa, 
        Statistics_lambda, 
        Statistics_entropy)>

<!-- Number of sequences in BLAST db -->
<!ELEMENT Statistics_db-num (%INTEGER;)>

<!-- Length of BLAST db -->
<!ELEMENT Statistics_db-len (%INTEGER;)>

<!-- Effective HSP length -->
<!ELEMENT Statistics_hsp-len (%INTEGER;)>

<!-- Effective search space -->
<!ELEMENT Statistics_eff-space (%REAL;)>

<!-- Karlin-Altschul parameter K -->
<!ELEMENT Statistics_kappa (%REAL;)>

<!-- Karlin-Altschul parameter Lambda -->
<!ELEMENT Statistics_lambda (%REAL;)>

<!-- Karlin-Altschul parameter H -->
<!ELEMENT Statistics_entropy (%REAL;)>


<!ELEMENT Hit (
        Hit_num, 
        Hit_id, 
        Hit_def, 
        Hit_accession, 
        Hit_len, 
        Hit_hsps?)>

<!-- hit number -->
<!ELEMENT Hit_num (%INTEGER;)>

<!-- SeqId of subject -->
<!ELEMENT Hit_id (#PCDATA)>

<!-- definition line of subject -->
<!ELEMENT Hit_def (#PCDATA)>

<!-- accession -->
<!ELEMENT Hit_accession (#PCDATA)>

<!-- length of subject -->
<!ELEMENT Hit_len (%INTEGER;)>

<!-- all HSP regions for the given subject -->
<!ELEMENT Hit_hsps (Hsp*)>


<!ELEMENT Hsp (
        Hsp_num, 
        Hsp_bit-score, 
        Hsp_score, 
        Hsp_evalue, 
        Hsp_query-from, 
        Hsp_query-to, 
        Hsp_hit-from, 
        Hsp_hit-to, 
        Hsp_pattern-from?, 
        Hsp_pattern-to?, 
        Hsp_query-frame?, 
        Hsp_hit-frame?, 
        Hsp_identity?, 
        Hsp_positive?, 
        Hsp_gaps?, 
        Hsp_align-len?, 
        Hsp_density?, 
        Hsp_qseq, 
        Hsp_hseq, 
        Hsp_midline?)>

<!-- HSP number -->
<!ELEMENT Hsp_num (%INTEGER;)>

<!-- score (in bits) of HSP -->
<!ELEMENT Hsp_bit-score (%REAL;)>

<!-- score of HSP -->
<!ELEMENT Hsp_score (%REAL;)>

<!-- e-value of HSP -->
<!ELEMENT Hsp_evalue (%REAL;)>

<!-- start of HSP in query -->
<!ELEMENT Hsp_query-from (%INTEGER;)>

<!-- end of HSP -->
<!ELEMENT Hsp_query-to (%INTEGER;)>

<!-- start of HSP in subject -->
<!ELEMENT Hsp_hit-from (%INTEGER;)>

<!-- end of HSP in subject -->
<!ELEMENT Hsp_hit-to (%INTEGER;)>

<!-- start of PHI-BLAST pattern -->
<!ELEMENT Hsp_pattern-from (%INTEGER;)>

<!-- end of PHI-BLAST pattern -->
<!ELEMENT Hsp_pattern-to (%INTEGER;)>

<!-- translation frame of query -->
<!ELEMENT Hsp_query-frame (%INTEGER;)>

<!-- translation frame of subject -->
<!ELEMENT Hsp_hit-frame (%INTEGER;)>

<!-- number of identities in HSP -->
<!ELEMENT Hsp_identity (%INTEGER;)>

<!-- number of positives in HSP -->
<!ELEMENT Hsp_positive (%INTEGER;)>

<!-- number of gaps in HSP -->
<!ELEMENT Hsp_gaps (%INTEGER;)>

<!-- length of the alignment used -->
<!ELEMENT Hsp_align-len (%INTEGER;)>

<!-- score density -->
<!ELEMENT Hsp_density (%INTEGER;)>

<!-- alignment string for the query (with gaps) -->
<!ELEMENT Hsp_qseq (#PCDATA)>

<!-- alignment string for subject (with gaps) -->
<!ELEMENT Hsp_hseq (#PCDATA)>

<!-- formating middle line -->
<!ELEMENT Hsp_midline (#PCDATA)>

