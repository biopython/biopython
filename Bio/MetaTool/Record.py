"""Hold MetaTool data in a straightforward format.
"""

# standard modules
import numpy.oldnumeric.matrix as Matrix

# local stuff
import Bio.MetaTool

class Metabolite( dict ):

    def __init__( self, reaction_count, metabolite_name ):
        self.reaction_count = reaction_count
        self.metabolite_name = metabolite_name

    def __str__( self ):
        out = '%d        %s\n' % ( self.reaction_count, self.metabolite_name )
        return out

class MetaboliteRole( dict ):

    def __init__( self, met, consumed, built, irreversible_vector ):
        self.met = met
        self.consumed = consumed
        self.built = built
        self.irreversible_vector = irreversible_vector

    def __str__( self ):
        out =  '%s\n' % self.met
        out = out + 'consumed in %d reactions\n' % self.consumed
        out = out + 'built in %d reactions\n' % self.built
        out = out + str( self.irreversible_vector )
        out = out + '\n'
        return out



class PathwayTransform:

    def __init__( self ):
        self.matrix = None
        self.enzymes = []
        self.reactions = []
        self.irreversible_vector = []

    def __str__( self ):
        out = ''
        out = out + '\nMatrix\n'
        out = out + str( self.matrix )
        if( len( self.enzymes ) > 0 ):
            out = out + '\n Enzymes\n'
            for enzyme in self.enzymes:
                out = out + '%s\n' % enzyme
        if( len( self.reactions ) > 0 ):
            out = out + '\n Reactions\n'
            for reaction in self.reactions:
                out = out + '%s\n' % reaction
        if( len( self.irreversible_vector ) > 0 ):
            out = out + '\n\nIrreversible\n\n'
            for scalar in self.irreversible_vector:
                out = out + '%03d    ' % int( scalar )
        out = out + '\n'
        return out



class Record:

    """Hold MetaTool output information.
    """

    def __init__( self ):
        self.external_metabolites = []
        self.internal_metabolites = []
        self.unbalanced_metabolites = []
        self.branch_metabolites = []
        self.non_branch_metabolites = []
        self.sum_is_constant_lines = []
        self.stochiometric = PathwayTransform()
        self.kernel = PathwayTransform()
        self.subsets = PathwayTransform()
        self.reduced_system = PathwayTransform()
        self.convex_basis = PathwayTransform()
        self.conservation_relations = PathwayTransform()
        self.elementary_modes = PathwayTransform()

    def __str__( self ):
        out = ''
        out = out + 'Input file name: %s\n' % self.input_file_name
        out = out + 'Number of internal metabolites: %d\n' % self.num_int_metabolites
        out = out + 'Number of reactions: %d\n' % self.num_reactions
        out = out + '\n\nExternal Metabolites\n\ncount    name\n\n'
        for metabolite in self.external_metabolites:
            out = out + str( metabolite )
        out = out + '\n\nInternal Metabolites\n\ncount    name\n\n'
        for metabolite in self.internal_metabolites:
            out = out + str( metabolite )
        out = out + '\n\nBranch Metabolites\n\n'
        for metabolite in self.branch_metabolites:
            out = out + str( metabolite )
        out = out + '\n\nNon Branch Metabolites\n\n'
        for metabolite in self.non_branch_metabolites:
            out = out + str( metabolite )
        out = out + '\n\nStochiometric\n\n'
        out = out + str( self.stochiometric )
        if( len( self.unbalanced_metabolites ) > 0 ):
            out = out + '\n\nUnbalanced Metabolites\n\n'
            for metabolite in self.unbalanced_metabolites:
                out = out + '%s\n' % metabolite
        out = out + '\n\nKernel\n\n'
        out = out + str( self.kernel )
        out = out + '\n\nSubsets\n\n'
        out = out + str( self.subsets )
        out = out + '\n\nReduced System\n\n'
        out = out + str( self.reduced_system )
        out = out + '\n\nConvex Basis\n\n'
        out = out + str( self.convex_basis )
        out = out + '\n\nConservation Relations\n\n'
        out = out + str( self.conservation_relations )
        out = out + '\n'
        for line in self.sum_is_constant_lines:
            out = out + '%s\n' % line
        out = out + '\n\nElementary Modes\n\n'
        out = out + str( self.elementary_modes )
        return out

