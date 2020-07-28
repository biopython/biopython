<h1 id="chapter-xxx-codon-alignment">Chapter XXX Codon Alignment</h1>
<p>This chapter is about Codon Alignments, which is a special case of nucleotide alignment in which the trinucleotides correspond directly to amino acids in the translated protein product. Codon Alignment carries information that can be used for many evolutionary analysis.</p>
<p>This chapter has been divided into four parts to explain the codon alignment support in Biopython. First, a general introduction about the basic classes in <code>Bio.codonalign</code> will be given. Then, a typical procedure of how to obtain a codon alignment within Biopython is then discussed. Next, some simple applications of codon alignment, such as dN/dS ratio estimation and neutrality test and so forth will be covered. Finally, IO support of codon alignment will help user to conduct analysis that cannot be done within Biopython.</p>
<h2 id="x.1-codonseq-class">X.1 <code>CodonSeq</code> Class</h2>
<p><code>Bio.codonalign.CodonSeq</code> object is the base object in Codon Alignment. It is similar to <code>Bio.Seq</code> but with some extra attributes. To obtain a simple <code>CodonSeq</code> object, you just need to give a <code>str</code> object of nucleotide sequence whose length is a multiple of 3 (This can be violated if you have <code>rf_table</code> argument). For example:</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio.codonalign import CodonSeq
&gt;&gt;&gt; codon_seq = CodonSeq(&quot;AAATTTCCCGGG&quot;)
&gt;&gt;&gt; codon_seq
CodonSeq(&#39;AAATTTCCCGGG&#39;)</code></pre>
<p>An error will raise up if the input sequence is not a multiple of 3.</p>
<pre class="pycon"><code>&gt;&gt;&gt; codon_seq = CodonSeq(&quot;AAATTTCCCGG&quot;)
Traceback (most recent call last):
...
ValueError: Sequence length is not a multiple of three (i.e. a whole number of codons)</code></pre>
<p>The slice of <code>CodonSeq</code> is exactly the same with <code>Seq</code> and it will always return a <code>Seq</code> object if you sliced a <code>CodonSeq</code>. For example:</p>
<pre class="pycon"><code>&gt;&gt;&gt; codon_seq1 = CodonSeq(&quot;AAA---CCCGGG&quot;)
&gt;&gt;&gt; codon_seq1
CodonSeq(&#39;AAA---CCCGGG&#39;)
&gt;&gt;&gt; codon_seq1[:6]
Seq(&#39;AAA---&#39;)
&gt;&gt;&gt; codon_seq1[1:5]
Seq(&#39;AA--&#39;)</code></pre>
<p><code>CodonSeq</code> objects have a <code>rf_table</code> attribute that dictates how the <code>CodonSeq</code> will be translated by indicating the starting position of each codon in the sequence). This is useful if your sequence is known to have frameshift events or pseudogene that has insertion or deletion. You might notice that in the previous example, you haven’t specify the <code>rf_table</code> when initiate a <code>CodonSeq</code> object. In fact, <code>CodonSeq</code> object will automatically assign a <code>rf_table</code> to the <code>CodonSeq</code> if you don’t say anything about it.</p>
<pre class="pycon"><code>&gt;&gt;&gt; codon_seq1 = CodonSeq(&quot;AAACCCGGG&quot;)
&gt;&gt;&gt; codon_seq1
CodonSeq(&#39;AAACCCGGG&#39;)
&gt;&gt;&gt; codon_seq1.rf_table
[0, 3, 6]
&gt;&gt;&gt; codon_seq1.translate()
&#39;KPG&#39;
&gt;&gt;&gt; codon_seq2 = CodonSeq(&quot;AAACCCGG&quot;, rf_table=[0, 3, 5])
&gt;&gt;&gt; codon_seq2.rf_table
[0, 3, 5]
&gt;&gt;&gt; codon_seq2.translate()
&#39;KPR&#39;</code></pre>
<p>In the example, we didn’t assign <code>rf_table</code> to <code>codon_seq1</code>. By default, <code>CodonSeq</code> will automatically generate a <code>rf_table</code> to the coding sequence assuming no frameshift events. In this case, it is <code>[0, 3, 6]</code>, which means the first codon in the sequence starts at position 0, the second codon in the sequence starts at position 3, and the third codon in the sequence starts at position 6. In <code>codon_seq2</code>, we only have 8 nucleotides in the sequence, but with <code>rf_table</code> option specified. In this case, the third codon starts at the 5th position of the sequence rather than the 6th. And the <code>translate()</code> function will use the <code>rf_table</code> to get the translated amino acid sequence.</p>
<p>Another thing to keep in mind is that <code>rf_table</code> will only be applied to ungapped nucleotide sequence. This makes <code>rf_table</code> to be interchangeable between <code>CodonSeq</code> with the same sequence but different gaps inserted. For example,</p>
<pre class="pycon"><code>&gt;&gt;&gt; codon_seq1 = CodonSeq(&quot;AAACCC---GGG&quot;)
&gt;&gt;&gt; codon_seq1.rf_table
[0, 3, 6]
&gt;&gt;&gt; codon_seq1.translate()
&#39;KPG&#39;
&gt;&gt;&gt; codon_seq1.full_translate()
&#39;KP-G&#39;</code></pre>
<p>We can see that the <code>rf_table</code> of <code>codon_seq1</code> is still <code>[0, 3, 6]</code>, even though we have gaps added. The <code>translate()</code> function will skip the gaps and return the ungapped amino acid sequence. If gapped protein sequence is what you need, <code>full_translate()</code> comes to help.</p>
<p>It is also easy to convert <code>Seq</code> object to <code>CodonSeq</code> object, but it is the user’s responsibility to ensure all the necessary information is correct for a <code>CodonSeq</code> (mainly <code>rf_table</code>).</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio.Seq import Seq
&gt;&gt;&gt; codon_seq = CodonSeq()
&gt;&gt;&gt; seq = Seq(&#39;AAAAAA&#39;)
&gt;&gt;&gt; codon_seq.from_seq(seq)
CodonSeq(&#39;AAAAAA&#39;)
&gt;&gt;&gt; seq = Seq(&#39;AAAAA&#39;)
&gt;&gt;&gt; codon_seq.from_seq(seq)
Traceback (most recent call last):
...
ValueError: Sequence length is not a multiple of three (i.e. a whole number of codons)
&gt;&gt;&gt; codon_seq.from_seq(seq, rf_table=(0, 2))
CodonSeq(&#39;AAAAA&#39;)</code></pre>
<h2 id="x.2-codonalignment-class">X.2 <code>CodonAlignment</code> Class</h2>
<p>The <code>CodonAlignment</code> class is another new class in <code>Codon.Align</code>. Its aim is to store codon alignment data and apply various analysis upon it. Similar to <code>MultipleSeqAlignment</code>, you can use numpy style slice to a <code>CodonAlignment</code>. However, once you sliced, the returned result will always be a <code>MultipleSeqAlignment</code> object.</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio.codonalign import CodonSeq, CodonAlignment
&gt;&gt;&gt; from Bio.SeqRecord import SeqRecord
&gt;&gt;&gt; a = SeqRecord(CodonSeq(&quot;AAAACGTCG&quot;), id=&quot;Alpha&quot;)
&gt;&gt;&gt; b = SeqRecord(CodonSeq(&quot;AAA---TCG&quot;), id=&quot;Beta&quot;)
&gt;&gt;&gt; c = SeqRecord(CodonSeq(&quot;AAAAGGTGG&quot;), id=&quot;Gamma&quot;)
&gt;&gt;&gt; codon_aln = CodonAlignment([a, b, c])
&gt;&gt;&gt; print(codon_aln)
CodonAlignment with 3 rows and 9 columns (3 codons)
AAAACGTCG Alpha
AAA---TCG Beta
AAAAGGTGG Gamma
&gt;&gt;&gt; codon_aln[0]
SeqRecord(seq=CodonSeq(&#39;AAAACGTCG&#39;), id=&#39;Alpha&#39;, name=&#39;&lt;unknown name&gt;&#39;, description=&#39;&lt;unknown description&gt;&#39;, dbxrefs=[])
&gt;&gt;&gt; print(codon_aln[:, 3])
A-A
&gt;&gt;&gt; print(codon_aln[1:, 3:10])
Alignment with 2 rows and 6 columns
---TCG Beta
AGGTGG Gamma</code></pre>
<p>You can write out <code>CodonAlignment</code> object just as what you do with <code>MultipleSeqAlignment</code>.</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio import AlignIO
&gt;&gt;&gt; AlignIO.write(codon_aln, &#39;example.aln&#39;, &#39;clustal&#39;)
1</code></pre>
<p>An alignment file called <code>example.aln</code> can then be found in your current working directory. You can write <code>CodonAlignment</code> out in any MSA format that Biopython supports.</p>
<p>Currently, you are not able to read MSA data as a <code>CodonAlignment</code> object directly (because of dealing with <code>rf_table</code> issue for each sequence). However, you can read the alignment data in as a <code>MultipleSeqAlignment</code> object and convert them into <code>CodonAlignment</code> object using <code>from_msa()</code> class method. For example,</p>
<pre class="pycon"><code>&gt;&gt;&gt; aln = AlignIO.read(&#39;example.aln&#39;, &#39;clustal&#39;)
&gt;&gt;&gt; codon_aln = CodonAlignment()
&gt;&gt;&gt; print(codon_aln.from_msa(aln))
CodonAlignment with 3 rows and 9 columns (3 codons)
AAAACGTCG Alpha
AAA---TCG Beta
AAAAGGTGG Gamma</code></pre>
<p>Note, the <code>from_msa()</code> method assume there is no frameshift events occurs in your alignment. Its behavior is not guaranteed if your sequence contains frameshift events!!</p>
<p>There is a couple of methods that can be applied to <code>CodonAlignment</code> class for evolutionary analysis. We will cover them more in X.4.</p>
<h2 id="x.3-build-a-codon-alignment">X.3 Build a Codon Alignment</h2>
<p>Building a codon alignment is the first step of many evolutionary anaysis. But how to do that? <code>Bio.codonalign</code> provides you an easy function <code>build()</code> to achieve all. The data you need to prepare in advance is a protein alignment and a set of DNA sequences that can be translated into the protein sequences in the alignment.</p>
<p><code>codonalign.build</code> method requires two mandatory arguments. The first one should be a protein <code>MultipleSeqAlignment</code> object and the second one is a list of nucleotide <code>SeqRecord</code> object. By default, <code>codonalign.build</code> assumes the order of the alignment and nucleotide sequences are in the same. For example:</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio import codonalign
&gt;&gt;&gt; from Bio.Align import MultipleSeqAlignment
&gt;&gt;&gt; from Bio.SeqRecord import SeqRecord
&gt;&gt;&gt; from Bio.Seq import Seq
&gt;&gt;&gt; nucl1 = SeqRecord(Seq(&#39;AAATTTCCCGGG&#39;), id=&#39;nucl1&#39;)
&gt;&gt;&gt; nucl2 = SeqRecord(Seq(&#39;AAATTACCCGCG&#39;), id=&#39;nucl2&#39;)
&gt;&gt;&gt; nucl3 = SeqRecord(Seq(&#39;ATATTACCCGGG&#39;), id=&#39;nucl3&#39;)
&gt;&gt;&gt; prot1 = SeqRecord(nucl1.seq.translate(), id=&#39;prot1&#39;)
&gt;&gt;&gt; prot2 = SeqRecord(nucl2.seq.translate(), id=&#39;prot2&#39;)
&gt;&gt;&gt; prot3 = SeqRecord(nucl3.seq.translate(), id=&#39;prot3&#39;)
&gt;&gt;&gt; aln = MultipleSeqAlignment([prot1, prot2, prot3])
&gt;&gt;&gt; codon_aln = codonalign.build(aln, [nucl1, nucl2, nucl3])
&gt;&gt;&gt; print(codon_aln)
CodonAlignment with 3 rows and 12 columns (4 codons)
AAATTTCCCGGG nucl1
AAATTACCCGCG nucl2
ATATTACCCGGG nucl3</code></pre>
<p>In the above example, <code>codonalign.build</code> will try to match <code>nucl1</code> with <code>prot1</code>, <code>nucl2</code> with <code>prot2</code> and <code>nucl3</code> with <code>prot3</code>, i.e., assuming the order of records in <code>aln</code> and <code>[nucl1, nucl2, nucl3]</code> is the same.</p>
<p><code>codonalign.build</code> method is also able to handle key match. In this case, records with same id are paired. For example:</p>
<pre class="pycon"><code>&gt;&gt;&gt; nucl1 = SeqRecord(Seq(&#39;AAATTTCCCGGG&#39;), id=&#39;nucl1&#39;)
&gt;&gt;&gt; nucl2 = SeqRecord(Seq(&#39;AAATTACCCGCG&#39;), id=&#39;nucl2&#39;)
&gt;&gt;&gt; nucl3 = SeqRecord(Seq(&#39;ATATTACCCGGG&#39;), id=&#39;nucl3&#39;)
&gt;&gt;&gt; prot1 = SeqRecord(nucl1.seq.translate(), id=&#39;prot1&#39;)
&gt;&gt;&gt; prot2 = SeqRecord(nucl2.seq.translate(), id=&#39;prot2&#39;)
&gt;&gt;&gt; prot3 = SeqRecord(nucl3.seq.translate(), id=&#39;prot3&#39;)
&gt;&gt;&gt; aln = MultipleSeqAlignment([prot1, prot2, prot3])
&gt;&gt;&gt; nucl = {&#39;prot1&#39;: nucl1, &#39;prot2&#39;: nucl2, &#39;prot3&#39;: nucl3}
&gt;&gt;&gt; codon_aln = codonalign.build(aln, nucl)
&gt;&gt;&gt; print(codon_aln)
CodonAlignment with 3 rows and 12 columns (4 codons)
AAATTTCCCGGG nucl1
AAATTACCCGCG nucl2
ATATTACCCGGG nucl3</code></pre>
<p>This option is useful if you read nucleotide sequences using <code>SeqIO.index</code> method, in which case the nucleotide dict with be generated automatically.</p>
<p>Sometimes, you are neither not able to ensure the same order or the same id. <code>codonalign.build</code> method provides you an manual approach to tell the program nucleotide sequence and protein sequence correspondance by generating a <code>corr_dict</code>. <code>corr_dict</code> should be a dictionary that uses protein record id as key and nucleotide record id as item. Let’s look at an example:</p>
<pre class="pycon"><code>&gt;&gt;&gt; nucl1 = SeqRecord(Seq(&#39;AAATTTCCCGGG&#39;), id=&#39;nucl1&#39;)
&gt;&gt;&gt; nucl2 = SeqRecord(Seq(&#39;AAATTACCCGCG&#39;), id=&#39;nucl2&#39;)
&gt;&gt;&gt; nucl3 = SeqRecord(Seq(&#39;ATATTACCCGGG&#39;), id=&#39;nucl3&#39;)
&gt;&gt;&gt; prot1 = SeqRecord(nucl1.seq.translate(), id=&#39;prot1&#39;)
&gt;&gt;&gt; prot2 = SeqRecord(nucl2.seq.translate(), id=&#39;prot2&#39;)
&gt;&gt;&gt; prot3 = SeqRecord(nucl3.seq.translate(), id=&#39;prot3&#39;)
&gt;&gt;&gt; aln = MultipleSeqAlignment([prot1, prot2, prot3])
&gt;&gt;&gt; corr_dict = {&#39;prot1&#39;: &#39;nucl1&#39;, &#39;prot2&#39;: &#39;nucl2&#39;, &#39;prot3&#39;: &#39;nucl3&#39;}
&gt;&gt;&gt; codon_aln = codonalign.build(aln, [nucl3, nucl1, nucl2], corr_dict=corr_dict)
&gt;&gt;&gt; print(codon_aln)
CodonAlignment with 3 rows and 12 columns (4 codons)
AAATTTCCCGGG nucl1
AAATTACCCGCG nucl2
ATATTACCCGGG nucl3</code></pre>
<p>We can see, even though the second argument of <code>codonalign.build</code> is not in the same order with <code>aln</code> in the above example, the <code>corr_dict</code> tells the program to pair protein records and nucleotide records. And we are still able to obtain the correct <code>codonalignment</code> object.</p>
<p>The underlying algorithm of <code>codonalign.build</code> method is very similar to <code>pal2nal</code> (a very famous perl script to build codon alignment). <code>codonalign.build</code> will first translate protein sequences into a long degenerate regular expression and tries to find a match in its corresponding nucleotide sequence. When translation fails, it divides protein sequence into several small anchors and tries to match each anchor to the nucleotide sequence to figure out where the mismatch and frameshift events lie. Other options available for <code>codonalign.build</code> includes <code>anchor_len</code> (default 10) and <code>max_score</code> (maximum tolerance of unexpected events, default 10). You may want to refer the Biopython build-in help to get more information about these options.</p>
<p>Now let’s look at a real example of building codon alignment. Here we will use epidermal growth factor (EGFR) gene to demonstrate how to obtain codon alignment. To reduce your effort, we have already collected EGFR sequences for Homo sapiens, Bos taurus, Rattus norvegicus, Sus scrofa and Drosophila melanogaster. The three files used in this example (<code>egfr_nucl.fa</code> with the nucleotide sequences of EGFR, <code>egfr_pro.aln</code> with the EGFR protein sequence alignment in <code>clustal</code> format, and <code>egfr_id</code> with the id correspondance between protein records and nucleotide records) is available from the ‘Tests/codonalign‘ directory in the Biopython distribution. You can then try the following code (make sure the files are in your current python working directory):</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio import SeqIO, AlignIO
&gt;&gt;&gt; nucl = SeqIO.parse(&#39;egfr_nucl.fa&#39;, &#39;fasta&#39;)
&gt;&gt;&gt; prot = AlignIO.read(&#39;egfr_pro.aln&#39;, &#39;clustal&#39;)
&gt;&gt;&gt; id_corr = {i.split()[0]: i.split()[1] for i in open(&#39;egfr_id&#39;).readlines()}
&gt;&gt;&gt; aln = codonalign.build(prot, nucl, corr_dict=id_corr)
/biopython/Bio/codonalign/__init__.py:568: UserWarning: gi|47522840|ref|NP_999172.1|(L 449) does not correspond to gi|47522839|ref|NM_214007.1|(ATG)
  % (pro.id, aa, aa_num, nucl.id, this_codon))
&gt;&gt;&gt; print(aln)
CodonAlignment with 6 rows and 4446 columns (1482 codons)
ATGATGATTATCAGCATGTGGATGAGCATATCGCGAGGATTGTGGGACAGCAGCTCC...GTG gi|24657088|ref|NM_057410.3|
---------------------ATGCTGCTGCGACGGCGCAACGGCCCCTGCCCCTTC...GTG gi|24657104|ref|NM_057411.3|
------------------------------ATGAAAAAGCACGAG------------...GCC gi|302179500|gb|HM749883.1|
------------------------------ATGCGACGCTCCTGGGCGGGCGGCGCC...GCA gi|47522839|ref|NM_214007.1|
------------------------------ATGCGACCCTCCGGGACGGCCGGGGCA...GCA gi|41327737|ref|NM_005228.3|
------------------------------ATGCGACCCTCAGGGACTGCGAGAACC...GCA gi|6478867|gb|M37394.2|RATEGFR</code></pre>
<p>We can see, while building the codon alignment a mismatch event is found. And this is shown as a UserWarning.</p>
<h2 id="x.4-codon-alignment-application">X.4 Codon Alignment Application</h2>
<p>The most important application of codon alignment is to estimate nonsynonymous substitutions per site (dN) and synonymous substitutions per site (dS). <code>codonalign</code> currently support three counting based methods (NG86, LWL85, YN00) and maximum likelihood method to estimate dN and dS. The function to conduct dN, dS estimation is called <code>cal_dn_ds</code>. When you obtained a codon alignment, it is quite easy to calculate dN and dS. For example (assuming you have EGFR codon alignmnet in the python working space):</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio.codonalign.codonseq import cal_dn_ds
&gt;&gt;&gt; print(aln)
CodonAlignment with 6 rows and 4446 columns (1482 codons)
ATGATGATTATCAGCATGTGGATGAGCATATCGCGAGGATTGTGGGACAGCAGCTCC...GTG gi|24657088|ref|NM_057410.3|
---------------------ATGCTGCTGCGACGGCGCAACGGCCCCTGCCCCTTC...GTG gi|24657104|ref|NM_057411.3|
------------------------------ATGAAAAAGCACGAG------------...GCC gi|302179500|gb|HM749883.1|
------------------------------ATGCGACGCTCCTGGGCGGGCGGCGCC...GCA gi|47522839|ref|NM_214007.1|
------------------------------ATGCGACCCTCCGGGACGGCCGGGGCA...GCA gi|41327737|ref|NM_005228.3|
------------------------------ATGCGACCCTCAGGGACTGCGAGAACC...GCA gi|6478867|gb|M37394.2|RATEGFR
&gt;&gt;&gt; dN, dS = cal_dn_ds(aln[0], aln[1], method=&#39;NG86&#39;)
&gt;&gt;&gt; print(dN, dS)
0.0209078305058 0.0178371876389
&gt;&gt;&gt; dN, dS = cal_dn_ds(aln[0], aln[1], method=&#39;LWL85&#39;)
&gt;&gt;&gt; print(dN, dS)
0.0203061425453 0.0163935691992
&gt;&gt;&gt; dN, dS = cal_dn_ds(aln[0], aln[1], method=&#39;YN00&#39;)
&gt;&gt;&gt; print(dN, dS)
0.0198195580321 0.0221560648799
&gt;&gt;&gt; dN, dS = cal_dn_ds(aln[0], aln[1], method=&#39;ML&#39;)
&gt;&gt;&gt; print(dN, dS)
0.0193877676103 0.0217247139962</code></pre>
<p>If you are using maximum likelihood methdo to estimate dN and dS, you are also able to specify equilibrium codon frequency to <code>cfreq</code> argument. Available options include <code>F1x4</code>, <code>F3x4</code> and <code>F61</code>.</p>
<p>It is also possible to get dN and dS matrix or a tree from a <code>CodonAlignment</code> object.</p>
<pre class="pycon"><code>&gt;&gt;&gt; dn_matrix, ds_matrix = aln.get_dn_ds_matrix()
&gt;&gt;&gt; print(dn_matrix)
gi|24657088|ref|NM_057410.3|    0
gi|24657104|ref|NM_057411.3|    0.0209078305058 0
gi|302179500|gb|HM749883.1|     0.611523924924  0.61022032668   0
gi|47522839|ref|NM_214007.1|    0.614035083563  0.60401686212   0.0411803504059 0
gi|41327737|ref|NM_005228.3|    0.61415325314   0.60182631356   0.0670105144563 0.0614703609541 0
gi|6478867|gb|M37394.2|RATEGFR  0.61870883409   0.606868724887  0.0738690303483 0.0735789092792 0.0517984707257 0
gi|24657088|ref|NM_057410.3|    gi|24657104|ref|NM_057411.3|    gi|302179500|gb|HM749883.1| gi|47522839|ref|NM_214007.1|    gi|41327737|ref|NM_005228.3|    gi|6478867|gb|M37394.2|RATEGFR
&gt;&gt;&gt; dn_tree, ds_tree = aln.get_dn_ds_tree()
&gt;&gt;&gt; print(dn_tree)
Tree(rooted=True)
    Clade(branch_length=0, name=&#39;Inner5&#39;)
        Clade(branch_length=0.279185347322, name=&#39;Inner4&#39;)
            Clade(branch_length=0.00859186651689, name=&#39;Inner3&#39;)
                Clade(branch_length=0.0258992353629, name=&#39;gi|6478867|gb|M37394.2|RATEGFR&#39;)
                Clade(branch_length=0.0258992353629, name=&#39;gi|41327737|ref|NM_005228.3|&#39;)
            Clade(branch_length=0.0139009266768, name=&#39;Inner2&#39;)
                Clade(branch_length=0.020590175203, name=&#39;gi|47522839|ref|NM_214007.1|&#39;)
                Clade(branch_length=0.020590175203, name=&#39;gi|302179500|gb|HM749883.1|&#39;)
        Clade(branch_length=0.294630667432, name=&#39;Inner1&#39;)
            Clade(branch_length=0.0104539152529, name=&#39;gi|24657104|ref|NM_057411.3|&#39;)
            Clade(branch_length=0.0104539152529, name=&#39;gi|24657088|ref|NM_057410.3|&#39;)</code></pre>
<p>Another application of codon alignment that <code>codonalign</code> supports is Mcdonald-Kreitman test. This test compares the within species synonymous substitutions and nonsynonymous substitutions and between species synonymous substitutions and nonsynonymous substitutions to see if they are from the same evolutionary process. The test requires gene sequences sampled from different individuals of the same species. In the following example, we will use Adh gene from fluit fly to demonstrate how to conduct the test. The data includes 11 individuals from D. melanogaster, 4 individuals from D. simulans and 12 individuals from D. yakuba. The data is available in the ‘Tests/codonalign‘ directory in the Biopython distribution. A function called <code>mktest</code> will be used for the test.</p>
<pre class="pycon"><code>&gt;&gt;&gt; from Bio import SeqIO, AlignIO
&gt;&gt;&gt; from Bio.codonalign import build
&gt;&gt;&gt; from Bio.codonalign.codonalignment import mktest
&gt;&gt;&gt; pro_aln = AlignIO.read(&#39;adh.aln&#39;, &#39;clustal&#39;)
&gt;&gt;&gt; p = SeqIO.index(&#39;drosophilla.fasta&#39;, &#39;fasta&#39;)
&gt;&gt;&gt; codon_aln = build(pro_aln, p)
&gt;&gt;&gt; print(codon_aln)
CodonAlignment with 27 rows and 768 columns (256 codons)
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9217|emb|X57365.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9219|emb|X57366.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9221|emb|X57367.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9223|emb|X57368.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9225|emb|X57369.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9227|emb|X57370.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9229|emb|X57371.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9231|emb|X57372.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9233|emb|X57373.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9235|emb|X57374.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9237|emb|X57375.1|
ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|9239|emb|X57376.1|
ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATT...ATC gi|9097|emb|X57361.1|
ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATT...ATC gi|9099|emb|X57362.1|
ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATT...ATC gi|9101|emb|X57363.1|
ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATC...ATC gi|9103|emb|X57364.1|
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|156879|gb|M17837.1|DROADHCK
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|156877|gb|M17836.1|DROADHCJ
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|156875|gb|M17835.1|DROADHCI
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|156873|gb|M17834.1|DROADHCH
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|156871|gb|M17833.1|DROADHCG
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATT...ATC gi|156863|gb|M19547.1|DROADHCC
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|156869|gb|M17832.1|DROADHCF
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATT...ATC gi|156867|gb|M17831.1|DROADHCE
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATT...ATC gi|156865|gb|M17830.1|DROADHCD
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATT...ATC gi|156861|gb|M17828.1|DROADHCB
ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATT...ATC gi|156859|gb|M17827.1|DROADHCA

&gt;&gt;&gt; print(mktest([codon_aln[1:12], codon_aln[12:16], codon_aln[16:]]))
0.00206457257254</code></pre>
<p>In the above example, <code>codon_aln[1:12]</code> belongs to D. melanogaster, <code>codon_aln[12:16]</code> belongs to D. simulans and <code>codon_aln[16:]</code> belongs to D. yakuba. <code>mktest</code> will return the p-value of the test. We can see in this case, 0.00206 &lt;&lt; 0.01, therefore, the gene is under strong negative selection according to MK test.</p>
<h2 id="x.4-future-development">X.4 Future Development</h2>
<p>Because of the limited time frame for Google Summer of Code project, some of the functions in <code>codonalign</code> is not tested comprehensively. In the following days, I will continue perfect the code and several new features will be added. I am always welcome to hear your suggestions and feature request. You are also highly encouraged to contribute to the existing code. Please do not hesitable to email me (zruan1991 at gmail dot com) when you have novel ideas that can make the code better.</p>
