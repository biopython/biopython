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

positive_pKs = { 'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H':5.98 }
negative_pKs = { 'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C':9.0, 'Y':10.0 }
pKcterminal= {'D':4.55, 'E':4.75}
pKnterminal = {'A':7.59, 'M':7.0, 'S':6.93, 'P':8.36, 'T':6.82, 'V':7.44, 'E':7.7}
charged_aas = ('K', 'R', 'H', 'D', 'E', 'C', 'Y')

# access this module through ProtParam.ProteinAnalysis class.
# first make a ProteinAnalysis object and then call its isoelectric_point method.
class IsoelectricPoint(object):
    def __init__(self, ProteinSequence, AminoAcidsContent):
        self.sequence = ProteinSequence        
        self.charged_aas_content = self._select_charged(AminoAcidsContent)

    # This function creates a dictionary with the contents of each charged aa, 
    # plus Cterm and Nterm.
    def _select_charged(self, AminoAcidsContent):
        charged = {}    
        for aa in charged_aas:
            charged[aa] = float(AminoAcidsContent[aa])
        charged['Nterm'] = 1.0
        charged['Cterm'] = 1.0
        return charged

    #This function calculates the total charge of the protein at a given pH.
    def _chargeR(self, pH, pos_pKs, neg_pKs):
        PositiveCharge = 0.0
        for aa, pK in pos_pKs.iteritems():         
             CR = 10**(pK-pH)
             partial_charge = CR/(CR+1.0)
             PositiveCharge += self.charged_aas_content[aa] * partial_charge 

        NegativeCharge = 0.0
        for aa, pK in neg_pKs.iteritems():         
             CR = 10**(pH-pK)
             partial_charge = CR/(CR+1.0)
             NegativeCharge += self.charged_aas_content[aa] * partial_charge 

        return PositiveCharge - NegativeCharge       

    # This is the action function, it tries different pH until the charge of the protein is 0 (or close).
    def pi(self):        
        pos_pKs = dict(positive_pKs)
        neg_pKs = dict(negative_pKs)
        nterm = self.sequence[0]
        cterm = self.sequence[-1]    
        if nterm in pKnterminal.keys():
            pos_pKs['Nterm'] = pKnterminal[nterm]
        if cterm in pKcterminal.keys():
            neg_pKs['Cterm'] = pKcterminal[cterm]

        # Bracket between pH1 and pH2
        pH = 7.0
        Charge = self._chargeR(pH, pos_pKs, neg_pKs)
        if Charge > 0.0:
            pH1 = pH
            Charge1 = Charge
            while Charge1 > 0.0:
                pH = pH1 + 1.0
                Charge = self._chargeR(pH, pos_pKs, neg_pKs)
                if Charge > 0.0:
                    pH1 = pH
                    Charge1 = Charge
                else:
                    pH2 = pH
                    Charge2 = Charge
                    break
        else:
            pH2 = pH
            Charge2 = Charge
            while Charge2 < 0.0:
                pH = pH2 - 1.0
                Charge = self._chargeR(pH, pos_pKs, neg_pKs)
                if Charge < 0.0:
                    pH2 = pH
                    Charge2 = Charge
                else:
                    pH1 = pH
                    Charge1 = Charge
                    break

        # Bisection
        while pH2 - pH1 > 0.0001 and Charge!=0.0:
            pH = (pH1 + pH2) / 2.0
            Charge = self._chargeR(pH, pos_pKs, neg_pKs)
            if Charge > 0.0:
                pH1 = pH
                Charge1 = Charge
            else:
                pH2 = pH
                Charge2 = Charge

        return pH
