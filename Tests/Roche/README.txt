Created with Roche tools from E3MFGYR02.sff, extracting 10 random reads:

* E3MFGYR02_random_10_reads.sff - has an .mft1.00 style index and manifest
* E3MFGYR02_random_10_reads.fasta
* E3MFGYR02_random_10_reads.qual
* E3MFGYR02_random_10_reads_no_trim.fasta
* E3MFGYR02_random_10_reads_no_trim.qual
* E3MFGYR02_no_manifest.sff - has an .srt1.00 style index and no manifest

Modified with a dummy index block in various positions, SFF software
should ignore an unknown index block gracefully. The Roche tools are
happy to read these files:

* E3MFGYR02_alt_index_at_start.sff
* E3MFGYR02_alt_index_in_middle.sff
* E3MFGYR02_alt_index_at_end.sff

Modified to move the ".mft1.00" index block. The Roche tools do not
like these files, and Roche told me in an email exchange that their
(undocumented) index blocks should only be at the end of the file.
However, I still found these useful for testing:

* E3MFGYR02_index_at_start.sff
* E3MFGYR02_index_in_middle.sff

Additional SFF files modified from real data to rename the reads. The
Roche tools are happy to read these files. They make excellent tests
for parsing the index block since the read names are different lengths:

* greek.sff  - has an .srt1.00 style index
* paired.sff - has an .mft1.00 style index and manifest

