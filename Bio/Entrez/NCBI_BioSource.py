# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results specified by NCBI's DTD file
# NCBI_BioSource.mod.dtd (04/10/2008 16:04:22).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = None

booleans = (
)

integers = (
    "BioSource_origin",		# (%INTEGER;)
				# ATTLIST value (unknown |
				#                natural |
				#                natmut |
				#                mut |
				#                artificial |
				#                synthetic |
				#                other) #IMPLIED
    "BioSource_genome",		# (%INTEGER;)
    "SubSource_subtype",	# (%INTEGER;)
				# ATTLIST value (chromosome |
        			#                map |
        			#                clone |
        			#                subclone |
        			#                haplotype |
        			#                genotype |
        			#                sex |
        			#                cell-line |
        			#                cell-type |
        			#                tissue-type |
        			#                clone-lib |
        			#                dev-stage |
        			#                frequency |
        			#                germline |
        			#                rearranged |
        			#                lab-host |
        			#                pop-variant |
        			#                tissue-lib |
        			#                plasmid-name |
        			#                transposon-name |
        			#                insertion-seq-name |
        			#                plastid-name |
        			#                country |
        			#                segment |
        			#                endogenous-virus-name |
        			#                transgenic |
        			#                environmental-sample |
        			#                isolation-source |
        			#                lat-lon |
        			#                collection-date |
        			#                collected-by |
        			#                identified-by |
        			#                fwd-primer-seq |
        			#                rev-primer-seq |
        			#                fwd-primer-name |
        			#                rev-primer-name |
        			#                metagenomic |
        			#                other) #IMPLIED
)

strings = (
    "BioSource_is-focus",	# EMPTY
    "SubSource_name",		# (#PCDATA)
    "SubSource_attrib",		# (#PCDATA)
)

lists = (
    "BioSource_subtype",	# (SubSource*)
)

dictionaries = (
    "BioSource",	# (BioSource_genome?, 
			#  BioSource_origin?, 
			#  BioSource_org, 
			#  BioSource_subtype?, 
			#  BioSource_is-focus?)
			# ATTLIST value (unknown |
			#                genomic |
			#                chloroplast |
			#                chromoplast |
			#                kinetoplast |
			#                mitochondrion |
			#                plastid |
			#                macronuclear |
			#                extrachrom |
			#                plasmid |
			#                transposon |
			#                insertion-seq |
			#                cyanelle |
			#                proviral |
			#                virion |
			#                nucleomorph |
			#                apicoplast |
			#                leucoplast |
			#                proplastid |
			#                endogenous-virus |
			#                hydrogenosome |
			#                chromosome |
			#                chromatophore) #IMPLIED
    "BioSource_org",	# (Org-ref)
    "SubSource",	# (SubSource_subtype, SubSource_name, SubSource_attrib?)
)

structures = {}

items = ()

def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path = self.path[:-1]
