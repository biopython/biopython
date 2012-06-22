Files ex1.* were originally from the pysam project. File ex1.fa
contains two sequences cut from the human genome build36. They
were extracted with command:

  samtools faidx human_b36.fa 2:2043966-2045540 20:67967-69550

Sequence names were changed manually for simplicity. File ex1.sam.gz
contains MAQ alignments extracted with:

  (samtools view NA18507_maq.bam 2:2044001-2045500;
   samtools view NA18507_maq.bam 20:68001-69500)

and processed with `samtools fixmate' to make it self-consistent as a
standalone alignment, and then converted to ex1.bam using `samtools
import'.

Files ex1_header.* were derived from ex1.* by adding a SAM header
(@HD and @SQ lines). The contents of File ex1_refresh.bam are the
same as ex1_header.bam but it was generated with a newer version
of samtools, and uses a different BGZF block strategy. The older
file ex1_header.bam has all its blocks the same size (64KB), the
newer file ex1_refresh.bam gives the header its own block and also
avoids splitting reads between blocks.
