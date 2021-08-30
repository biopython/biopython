#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Restriction Digest Enzymes.

Examples
--------
    >>> from Bio.Seq import Seq
    >>> from Bio.Restriction import *
    >>> pBs_mcs = 'GGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTG'
    >>> pBs_mcs += 'CAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGCCACCGCGGTGGAGCTC'
    >>> seq = Seq(pBs_mcs)  # Multiple-cloning site of pBluescript SK(-)
    >>> a = Analysis(AllEnzymes, seq)
    >>> a.print_that()              # no argument -> print all the results
    AbaSI      :  10, 12, 13, 16, 17, 18, 19, 20, 22, 23, 24, 25, 25, 26, 27...
    BmeDI      :  2, 7, 8, 8, 9, 9, 13, 14, 15, 16, 17, 18, 19, 19, 21, 21...
    YkrI       :  10, 12, 13, 16, 16, 17, 19, 20, 21, 22, 23, 24, 25, 25, 26...

    BmeDI      :  1, 2, 7, 8, 8, 9, 9, 13, 14, 15, 16, 17, 18, 19...
    AccII      :  98.
    AciI       :  86, 90, 96, 98...

    Enzymes which do not cut the sequence.

    AspLEI    BstHHI    CfoI      CviAII    FaeI      FaiI      FatI      GlaI
    HhaI      Hin1II    Hin6I     HinP1I    HpyCH4IV  HpySE526I Hsp92II   HspAI
    MaeII     MseI      NlaIII    SaqAI     TaiI      Tru1I     Tru9I...
    <BLANKLINE>
    >>> b = a.blunt()  # Analysis with blunt enzmyes
    >>> a.print_that(b)  # Print results for blunt cutters
    AccII      :  98.
    AfaI       :  4.
    AluBI      :  40, 106.
    AluI       :  40, 106.
    Bsh1236I   :  98.
    BshFI      :  10, 89.
    BsnI       :  10, 89.
    BspANI     :  10, 89...

    Enzymes which do not cut the sequence.

    FaiI      GlaI      CdiI      MlyI      SchI      SspD5I    AanI...
    <BLANKLINE>

"""  # noqa: W291, W293

from Bio.Restriction.Restriction import *  # noqa (legacy module arrangement)


#
#   OK can't put the following code in Bio.Restriction.__init__ unless
#   I put everything from Restriction in here.
#   or at least the RestrictionBatch class.
#
#   The reason for that is if I do that, I break the __contains__ method of
#   the RestrictionBatch in Restriction, which expect to find the name of
#   the enzymes in the locals() dictionary when evaluating string to see if
#   it is an enzyme.
#
#   This calls for some explanations I guess:
#       When testing for the presence of a Restriction enzyme in a
#       RestrictionBatch, the user can use:
#
#           1) a class of type 'RestrictionType'
#           2) a string of the name of the enzyme (its repr)
#               i.e:
#                   >>> from Bio.Restriction import RestrictionBatch, EcoRI
#                   >>> MyBatch = RestrictionBatch(EcoRI)
#                   >>> EcoRI in MyBatch        # the class EcoRI.
#                   True
#                   >>> 'EcoRI' in MyBatch      # a string representation
#                   True
#
#   OK, that's how it is supposed to work. And I find it quite useful.
#
#   Now if I leave the code here I got:
#                   >>> from Bio.Restriction import RestrictionBatch, EcoRI
#                   >>> MyBatch = RestrictionBatch(EcoRI)
#                   >>> EcoRI in MyBatch # the class EcoRI.
#                   True
#                   >>> 'EcoRI' in MyBatch   # a string.
#                   False

#   There is 5 ways to change that:
#       1) abandon the evaluation of string representation.
#       2) leave the code like that and hack something in RestrictionBatch.
#       3) Move back the code in Bio.Restriction.Restriction
#       4) Move RestrictionBatch here.
#       5) Remove Restriction.Restriction and move all the code in here
#
#   1) no fun in that.
#   2) there is a simpler way to do it.
#   3) I prefer to keep all the code together.
#   4) and 5) both are OK. Only a matter of preference.
#
#   So the following code has been moved back to Bio.Restriction.Restriction
#   For the user the results is transparent:
#   from Bio.Restriction import * works as before.
#

# ##
# ##   The restriction enzyme classes are created dynamically when the module
# ##   is imported. Here is the magic which allow the creation of the
# ##   restriction-enzyme classes.
# ##
# ##   The reason for the two dictionaries in Restriction_Dictionary
# ##   one for the types (which will be called pseudo-type as they really
# ##   correspond to the values that instances of RestrictionType can take)
# ##   and one for the enzymes is efficiency as the bases are evaluated
# ##   once per pseudo-type.
# ##
# ##   However Restriction is still a very inefficient module at import. But
# ##   remember that around 660 classes (which is more or less the size of
# ##   Rebase) have to be created dynamically. However, this processing take
# ##   place only once.
# ##   This inefficiency is however largely compensated by the use of metaclass
# ##   which provide a very efficient layout for the class themselves mostly
# ##   alleviating the need of if/else loops in the class methods.
# ##
# ##   It is essential to run Restriction with doc string optimisation (-OO
# ##   switch) as the doc string of 660 classes take a lot of processing.
# ##
# # CommOnly    = RestrictionBatch()    # commercial enzymes
# # NonComm     = RestrictionBatch()    # not available commercially
# # for TYPE, (bases, enzymes) in typedict.items():
# #    #
# #    #   The keys are the pseudo-types TYPE (stored as type1, type2...)
# #    #   The names are not important and are only present to differentiate
# #    #   the keys in the dict. All the pseudo-types are in fact
# #    #   RestrictionType. These names will not be used after and the pseudo-
# #    #   types are not kept in the locals() dictionary. It is therefore
# #    #   impossible to import them.
# #    #   Now, if you have look at the dictionary, you will see that not all
# #    #   the types are present as those without corresponding enzymes have
# #    #   been removed by Dictionary_Builder().
# #    #
# #    #   The values are tuples which contain
# #    #   as first element a tuple of bases (as string) and
# #    #   as second element the names of the enzymes.
# #    #
# #    #   First eval the bases.
# #    #
# #    bases = tuple(eval(x) for x in bases)
# #    #
# #    #   now create the particular value of RestrictionType for the classes
# #    #   in enzymes.
# #    #
# #    T = type.__new__(RestrictionType, 'RestrictionType', bases, {})
# #    for k in enzymes:
# #        #
# #        #   Now, we go through all the enzymes and assign them their type.
# #        #   enzymedict[k] contains the values of the attributes for this
# #        #   particular class (self.site, self.ovhg,....).
# #        #
# #        newenz = T(k, bases, enzymedict[k])
# #        #
# #        #   we add the enzymes to the corresponding batch.
# #        #
# #        #   No need to verify the enzyme is a RestrictionType -> add_nocheck
# #        #
# #        if newenz.is_comm() : CommOnly.add_nocheck(newenz)
# #        else : NonComm.add_nocheck(newenz)
# ##
# ##   AllEnzymes is a RestrictionBatch with all the enzymes from Rebase.
# ##
# # AllEnzymes = CommOnly | NonComm
# ##
# ##   Now, place the enzymes in locals so they can be imported.
# ##
# # names = [str(x) for x in AllEnzymes]
# # locals().update(dict(map(None, names, AllEnzymes)))
# ##
# ##   Limit what can be imported by from Restriction import *
# ##   Most of the classes here will never be used outside this module
# ##   (Defined,Palindromic...). It is still possible to request them
# ##   specifically
# ##
# ##   also delete the variable that are no longer needed.
# ##
# ##
# # __all__= ['Analysis', 'RestrictionBatch','AllEnzymes','CommOnly',
# #           'NonComm']+names
# # del k, x, enzymes, TYPE, bases, names
