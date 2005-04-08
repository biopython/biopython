# Copyright Yair Benita Y.Benita@pharm.uu.nl
# Biopython (http://biopython.org) license applies

"""Calculate isoelectric points of polypeptides using methods of Bjellqvist.

pK values and the methos are taken from:

* Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F., Sanchez,
J.-Ch., Frutiger, S. & Hochstrasser, D.F.
The focusing positions of polypeptides in immobilized pH gradients can be predicted
from their amino acid sequences. Electrophoresis 1993, 14, 1023-1031. 

* Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
Reference points for comparisons of two-dimensional maps of proteins from
different human cell types defined in a pH scale where isoelectric points correlate
with polypeptide compositions. Electrophoresis 1994, 15, 529-539.

I designed the algorithm according to a note by David L. Tabb, available at:
http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf

"""

pKinternal = {'D':4.05, 'E':4.45, 'H':5.98, 'C':9.0, 'K':10.0, 'Y':10.0, 'R':12.0}
pKcterminal= {'D':4.55, 'E':4.75}
pKnterminal = {'A':7.59, 'M':7.0, 'S':6.93, 'P':8.36, 'T':6.82, 'V':7.44, 'E':7.7}

neg = ['D', 'E', 'C', 'Y', 'Cterm']
pos = ['R', 'K', 'H', 'Nterm']

# access this module through ProtParam.ProteinAnalysis class.
# first make a ProteinAnalysis object and then call its isoelectric_point method.
class IsoelectricPoint:
	def __init__(self, ProteinSequence, AminoAcidsContent):
		self.sequence = ProteinSequence
		self.amino_acids_content = AminoAcidsContent
	
	
	#This function calculates the total charge of the protein at a given pH.		
	def _chargeR(self, pH, Cterm=0, Nterm=0):
		NegScore = 0.0
		PosScore = 0.0
		# make a local dictionary with the relevant amino acids
		ChargedAminoAcids={}
		for i in pKinternal.keys():
			ChargedAminoAcids[i] = self.amino_acids_content[i]
		
		for i in ChargedAminoAcids.keys():
			#first calculate the score per amino acid
			score = ChargedAminoAcids[i] * self._cr_ratio(pKinternal[i], pH, i)
			# determine if it goes to positive or negative
			if i in pos:
				PosScore += score
			elif i in neg:
				NegScore += score
				
		# now add the score for C and N terminals.
		if Nterm == 0:	
			PosScore += self._cr_ratio(7.5, pH ,'Nterm')
		else:
			PosScore += self._cr_ratio(Nterm, pH ,'Nterm')
		
		if Cterm == 0:
			NegScore += self._cr_ratio(3.55, pH ,'Cterm')
		else:
			NegScore += self._cr_ratio(Cterm, pH ,'Cterm')
			
		return PosScore-NegScore
	
	# This is the action function, it tries different pH until the charge of the protein is 0 (or close).
	def pi(self):		
		if self.sequence[0] in pKnterminal.keys():
			Nterm = pKnterminal[self.sequence[0]]
		else:
			Nterm = 0
			
		if self.sequence[-1] in pKcterminal.keys():
			Cterm = pKcterminal[self.sequence[-1]]
		else:
			Cterm = 0
			
		pH = 7.0
		CutBy = 3.5
		Charge = self._chargeR(pH, Cterm, Nterm)

		while CutBy > 0.0001: # this loop is usually good enough to get the job done.
			pH1=pH+CutBy
			pH2=pH-CutBy
			Charge1 = self._chargeR(pH1, Cterm, Nterm)
			Charge2 = self._chargeR(pH2, Cterm, Nterm)

			if abs(Charge1) < abs(Charge2):
				Charge = Charge1
				pH = pH1
			else:
				Charge = Charge2
				pH = pH2
			CutBy = CutBy/2.0
			
		
		# some times the charge is still higher or lower than zero. in such case we fix it in a step-wise way.
		if abs(Charge) > 0.1:
			while abs(Charge) > 0.1:
				if Charge > 0:
					pH += 0.001
				else:
					pH -= 0.001
				Charge = self._chargeR(pH, Cterm, Nterm)
				
		return pH
	
	# This function calculates the charge on one amino acid at a specific pH.	
	def _cr_ratio(self, pk, ph, i):
		if i in pos:
			CR=10**(pk-ph)
		elif i in neg:
			CR=10**(ph-pk)
		else:
			raise KeyError ("amino acid is neither positive nor negative. check dictionary.")
		return (CR/(CR+1.0))
	

