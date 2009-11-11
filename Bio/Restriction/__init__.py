#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""
Usage:
=====

    >>> from Rana.fts import fts    # 
    >>> from Rana.Vector import *   # Just a way to get a sequence.
    >>> from Bio.Seq import Seq     # Use your prefer method here.
    >>> pbr = fts(pBR322)           #
    >>> seq = Seq(str(pbr))         #
    >>>
    >>> from Bio.Restriction import *
    >>> a = Analysis(AllEnzymes, seq, linear=False)
    >>> b = a.blunt()
    >>> a.print_that()              # no argument -> print all the results
    AasI       :  2169, 2582.
    AatII      :  4289.
    Acc16I     :  263, 1359, 1457, 3589.
    ...
        More enzymes here.
    ...
    >>> b = a.without_site()
    >>> a.print_that(b, '', '\n Enzymes which do not cut pBR322.\n\n')

     Enzymes which do not cut pBR322.

    AarI      AatI      Acc65I    AcsI      AcvI      AdeI      AflII     AgeI      
    AhlI      AleI      AloI      ApaI      ApoI      AscI      AsiAI     AsiSI     
    Asp718I   AspA2I    AsuII     AvaIII    AvrII     AxyI      BaeI      BbrPI     
    BbvCI     BclI      BcuI      BfrBI     BfrI      BglII     BlnI      BlpI      
    BmgBI     BmgI      BplI      Bpu1102I  Bpu14I    BsaXI     Bse21I    BsePI     
    BseRI     BshTI     BsiWI     Bsp119I   Bsp120I   Bsp1407I  Bsp1720I  Bsp19I    
    BspT104I  BspTI     BsrGI     BssHI     BssHII    Bst98I    BstAUI    BstBI     
    BstEII    BstPI     BstSNI    BstXI     Bsu36I    BtrI      CciNI     CelII     
    Cfr42I    Cfr9I     CpoI      Csp45I    CspAI     CspCI     CspI      DraIII    
    DrdII     Ecl136II  Eco105I   Eco147I   Eco72I    Eco81I    Eco91I    EcoICRI   
    EcoO65I   EcoRI     EcoT22I   EspI      FalI      FbaI      FseI      FunII     
    HpaI      KpnI      Ksp22I    KspAI     KspI      MabI      MfeI      MluI      
    Mph1103I  MspCI     MssI      MunI      NcoI      NotI      NsiI      NspV      
    OliI      PacI      PaeR7I    PasI      PauI      PceI      Pfl23II   PinAI     
    PmaCI     PmeI      PmlI      Ppu10I    PsiI      Psp124BI  PspAI     PspCI     
    PspEI     PspLI     PspOMI    PspXI     PsrI      RleAI     Rsr2I     RsrII     
    SacI      SacII     SanDI     SauI      SbfI      SciI      SdaI      SexAI     
    SfiI      Sfr274I   Sfr303I   SfuI      SgfI      SgrBI     SlaI      SmaI      
    SmiI      SnaBI     SpeI      SplI      SrfI      Sse232I   Sse8387I  Sse8647I  
    SseBI     SspBI     SstI      StuI      SunI      SwaI      TliI      UthSI     
    Vha464I   XapI      XbaI      XcmI      XhoI      XmaCI     XmaI      XmaJI     
    Zsp2I     

    >>>
    """

from Bio.Restriction.Restriction import *
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
#   This call for some explanations I guess:
#       When testing for the presence of a Restriction enzyme in a
#       RestrictionBatch, the user can use:
#
#           1) a class of type 'RestrictionType'
#           2) a string of the name of the enzyme (it's repr)
#               i.e:
#                   >>> from Bio.Restriction import RestrictionBatch, EcoRI
#                   >>> MyBatch = RestrictionBatch(EcoRI)
#                   >>> #!/usr/bin/env python
#                   >>> EcoRI in MyBatch        # the class EcoRI.
#                   True
#                   >>>
#                   >>> 'EcoRI' in MyBatch      # a string representation
#                   True
#
#   OK, that's how it is suppose to work. And I find it quite useful.
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
#   So the following code has been moved back to Bio.Restricion.Restriction
#   For the user the results is transparent:
#   from Bio.Restriction import * works as before.
#
        
###
###   The restriction enzyme classes are created dynamically when the module is
###   imported. Here is the magic which allow the creation of the
###   restriction-enzyme classes.
###
###   The reason for the two dictionaries in Restriction_Dictionary
###   one for the types (which will be called pseudo-type as they really
###   correspond to the values that instances of RestrictionType can take)
###   and one for the enzymes is efficiency as the bases are evaluated
###   once per pseudo-type.
###
###   However Restriction is still a very inefficient module at import. But
###   remember that around 660 classes (which is more or less the size of Rebase)
###   have to be created dynamically. However, this processing take place only
###   once.
###   This inefficiency is however largely compensated by the use of metaclass
###   which provide a very efficient layout for the class themselves mostly
###   alleviating the need of if/else loops in the class methods.
###
###   It is essential to run Restriction with doc string optimisation (-OO switch)
###   as the doc string of 660 classes take a lot of processing.
###
##CommOnly    = RestrictionBatch()    # commercial enzymes
##NonComm     = RestrictionBatch()    # not available commercially
##for TYPE, (bases, enzymes) in typedict.iteritems():
##    #
##    #   The keys are the pseudo-types TYPE (stored as type1, type2...)
##    #   The names are not important and are only present to differentiate
##    #   the keys in the dict. All the pseudo-types are in fact RestrictionType.
##    #   These names will not be used after and the pseudo-types are not
##    #   kept in the locals() dictionary. It is therefore impossible to
##    #   import them.
##    #   Now, if you have look at the dictionary, you will see that not all the
##    #   types are present as those without corresponding enzymes have been
##    #   removed by Dictionary_Builder().
##    #
##    #   The values are tuples which contain
##    #   as first element a tuple of bases (as string) and
##    #   as second element the names of the enzymes.
##    #
##    #   First eval the bases.
##    #
##    bases = tuple([eval(x) for x in bases])
##    #
##    #   now create the particular value of RestrictionType for the classes
##    #   in enzymes.
##    #
##    T = type.__new__(RestrictionType, 'RestrictionType', bases, {})
##    for k in enzymes:
##        #
##        #   Now, we go through all the enzymes and assign them their type.
##        #   enzymedict[k] contains the values of the attributes for this
##        #   particular class (self.site, self.ovhg,....).
##        #
##        newenz = T(k, bases, enzymedict[k])
##        #
##        #   we add the enzymes to the corresponding batch.
##        #
##        #   No need to verify the enzyme is a RestrictionType -> add_nocheck
##        #
##        if newenz.is_comm() : CommOnly.add_nocheck(newenz)
##        else : NonComm.add_nocheck(newenz)
###
###   AllEnzymes is a RestrictionBatch with all the enzymes from Rebase.
###
##AllEnzymes = CommOnly | NonComm
###
###   Now, place the enzymes in locals so they can be imported.
###
##names = [str(x) for x in AllEnzymes]
##locals().update(dict(map(None, names, AllEnzymes)))
###
###   Limit what can be imported by from Restriction import *
###   Most of the classes here will never be used outside this module
###   (Defined,Palindromic...). It is still possible to request them specifically
###
###   also delete the variable that are no longer needed.
###   
###
##__all__=['Analysis', 'RestrictionBatch','AllEnzymes','CommOnly','NonComm']+names
##del k, x, enzymes, TYPE, bases, names
