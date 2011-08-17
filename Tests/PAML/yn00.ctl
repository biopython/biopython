      seqfile = PAML/alignment.phylip * sequence data file name
      outfile = temp.out           * main result file
      verbose = 1  * 1: detailed output (list sequences), 0: concise output

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

    weighting = 0  * weighting pathways between codons (0/1)?
   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)?
       ndata = 1


* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
