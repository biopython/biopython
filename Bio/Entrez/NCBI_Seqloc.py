# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results specified by NCBI's DTD file
# NCBI_Entrezgene.mod.dtd (04/10/2008 16:04:22).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = None

booleans = (
)

integers = (
    "Seq-point_point",		# (%INTEGER;)
    "Seq-interval_from",	# (%INTEGER;)
    "PDB-seq-id_chain",		# (%INTEGER;)
    "Giimport-id_id",		# (%INTEGER;)
    "Textseq-id_version",	# (%INTEGER;)
    "Patent-seq-id_seqid",	# (%INTEGER;)
    "Seq-id_gi",		# (%INTEGER;)
    "Seq-id_gibbsq"	,	# (%INTEGER;)
    "Packed-seqpnt_points_E",	# (%INTEGER;)
    "Seq-id_gibbmt",		# (%INTEGER;)
    "Seq-interval_to",		# (%INTEGER;)
)

strings = (
    "Giimport-id_db",		# (#PCDATA)
    "Giimport-id_release",	# (#PCDATA)
    "Na-strand",		# %ENUM;
				# ATTLIST value (unknown |
				#                plus |
				#                minus |
				#                both |
				#                both-rev |
				#                other
				#               ) #REQUIRED
    "PDB-mol-id",		# (#PCDATA)
    "Seq-loc_null",		# EMPTY
    "Textseq-id_name",		# (#PCDATA)
    "Textseq-id_accession",	# (#PCDATA)
    "Textseq-id_release",	# (#PCDATA)
)

lists = (
    "Packed-seqint",		# (Seq-interval*)
    "Packed-seqpnt_points",	# (Packed-seqpnt_points_E*)
    "Seq-loc-equiv",		# (Seq-loc*)
    "Seq-loc-mix",		# (Seq-loc*)
)

dictionaries = (
    "Giimport-id",	# (Giimport-id_id, 
			#  Giimport-id_db?, 
			#  Giimport-id_release?)
    "Packed-seqpnt",		# (Packed-seqpnt_strand?, 
				#  Packed-seqpnt_id, 
 	                       #  Packed-seqpnt_fuzz?, 
				#  Packed-seqpnt_points)
    "Packed-seqpnt_strand",	# (Na-strand)
    "Packed-seqpnt_id",		# (Seq-id)
    "Packed-seqpnt_fuzz",	# (Int-fuzz)
    "Patent-seq-id",	# (Patent-seq-id_seqid, Patent-seq-id_cit)
    "Patent-seq-id_cit",	# (Id-pat)
    "PDB-seq-id",	# (PDB-seq-id_mol, PDB-seq-id_chain?, PDB-seq-id_rel?)
    "PDB-seq-id_mol",	# (PDB-mol-id)
    "PDB-seq-id_rel",	# (Date)
    "Seq-id",		# (Seq-id_local |
			#  Seq-id_gibbsq |
			#  Seq-id_gibbmt |
			#  Seq-id_giim |
			#  Seq-id_genbank |
			#  Seq-id_embl |
			#  Seq-id_pir |
			#  Seq-id_swissprot |
			#  Seq-id_patent |
			#  Seq-id_other |
			#  Seq-id_general |
			#  Seq-id_gi |
			#  Seq-id_ddbj |
			#  Seq-id_prf |
			#  Seq-id_pdb |
			#  Seq-id_tpg |
			#  Seq-id_tpe |
			#  Seq-id_tpd |
			#  Seq-id_gpipe)
    "Seq-id_local",	# (Object-id)
    "Seq-id_giim",	# (Giimport-id)
    "Seq-id_genbank",	# (Textseq-id)
    "Seq-id_embl",	# (Textseq-id>
    "Seq-id_pir",	# (Textseq-id)
    "Seq-id_swissprot",	# (Textseq-id)
    "Seq-id_patent",	# (Patent-seq-id)
    "Seq-id_other",	# (Textseq-id)
    "Seq-id_general",	# (Dbtag)
    "Seq-id_ddbj",	# (Textseq-id)
    "Seq-id_prf",	# (Textseq-id)
    "Seq-id_pdb",	# (PDB-seq-id)
    "Seq-id_tpg",	# (Textseq-id)
    "Seq-id_tpe",	# (Textseq-id)
    "Seq-id_tpd",	# (Textseq-id)
    "Seq-id_gpipe",	# (Textseq-id)
    "Seq-interval",	# (Seq-interval_from,
			#  Seq-interval_to,
			#  Seq-interval_strand?,
			#  Seq-interval_id,
			#  Seq-interval_fuzz-from?,
			#  Seq-interval_fuzz-to?)
    "Seq-interval_fuzz-from",	# (Int-fuzz)
    "Seq-interval_fuzz-to",	# (Int-fuzz)
    "Seq-interval_id",		# (Seq-id)
    "Seq-interval_strand",	# (Na-strand)
    "Seq-loc",		# (Seq-loc_null | 
			#  Seq-loc_empty | 
			#  Seq-loc_whole | 
			#  Seq-loc_int | 
			#  Seq-loc_packed-int | 
			#  Seq-loc_pnt | 
			#  Seq-loc_packed-pnt | 
			#  Seq-loc_mix | 
			#  Seq-loc_equiv | 
			#  Seq-loc_bond | 
			#  Seq-loc_feat)
    "Seq-loc_bond",	# (Seq-bond)
    "Seq-loc_empty",	# (Seq-id)
    "Seq-loc_equiv",	# (Seq-loc-equiv)
    "Seq-loc_feat",	# (Feat-id)
    "Seq-loc_int",	# (Seq-interval)
    "Seq-loc_mix",	# (Seq-loc-mix)
    "Seq-loc_packed-int",	# (Packed-seqint)
    "Seq-loc_packed-pnt",	# (Packed-seqpnt)
    "Seq-loc_pnt",	# (Seq-point)
    "Seq-loc_whole",	# (Seq-id)
    "Seq-point",	# (Seq-point_point, 
			#  Seq-point_strand?, 
			#  Seq-point_id, 
			#  Seq-point_fuzz?)
    "Seq-point_strand",	# (Na-strand)>
    "Seq-point_id",	# (Seq-id)
    "Seq-point_fuzz",	# (Int-fuzz)
    "Seq-bond",		# (Seq-bond_a, Seq-bond_b?)
    "Seq-bond_a",	# (Seq-point)
    "Seq-bond_b",	# (Seq-point)
    "Textseq-id",	# (Textseq-id_name?, 
			#  Textseq-id_accession?, 
			#  Textseq-id_release?, 
			#  Textseq-id_version?)
)

structures = {}

items = ()

def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path = self.path[:-1]

