# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""BioPython Pathway module.

Bio.Pathway is a lightweight class library designed to support the following tasks:

 - Data interchange and preprocessing between pathway databases and analysis software.
 - Quick prototyping of pathway analysis algorithms

The basic object in the Bio.Pathway model is Interaction, which represents an arbitrary
interaction between any number of biochemical species.

Network objects are used to represent the connectivity between species in pathways
and reaction networks.

For applications where it is not necessary to explicitly represent network connectivity,
the specialized classes Reaction and System should be used in place of Interacton and
Network.

The Bio.Pathway classes, especially Interaction, are intentionally
designed to be very flexible. Their intended use are as wrappers around database
specific records, such as BIND objects. The value-added in this module is a
framework for representing collections of reactions in a way that supports
graph theoretic and numeric analysis.

Note: This module should be regarded as a prototype only. API changes are likely.
      Comments and feature requests are most welcome.

"""

from functools import reduce

from Bio.Pathway.Rep.MultiGraph import MultiGraph


class Reaction:
    """Abstraction for a biochemical transformation.

    This class represents a (potentially reversible) biochemical
    transformation of the type::

      a S1 + b S2 + ... --> c P1 + d P2 + ...

    where:
     - a, b, c, d ... are positive numeric stochiometric coefficients,
     - S1, S2, ... are substrates
     - P1, P2, ... are products

    A Reaction should be viewed as the net result of one or more individual
    reaction steps, where each step is potentially facilitated by a different
    catalyst. Support for 'Reaction algebra' will be added at some point in
    the future.

    Attributes:
     - reactants   -- dict of involved species to their stochiometric coefficients:
       reactants[S] = stochiometric constant for S
     - catalysts   -- list/tuple of tuples of catalysts required for this reaction
     - reversible  -- true iff reaction is reversible
     - data        -- reference to arbitrary additional data

    Invariants:
     - for all S in reactants: reactants[S] != 0
     - for all C in catalysts: catalysts[C] != 0

    """

    def __init__(self, reactants=None, catalysts=(), reversible=0, data=None):
        """Initialize a new Reaction object."""
        # enforce invariants on reactants:
        if reactants is None:
            self.reactants = {}
        else:
            self.reactants = reactants.copy()
            # loop over original, edit the copy
            for r, value in reactants.items():
                if value == 0:
                    del self.reactants[r]
        self.catalysts = sorted(set(catalysts))
        self.data = data
        self.reversible = reversible

    def __eq__(self, r):
        """Return true iff self is equal to r."""
        return (
            isinstance(r, Reaction)
            and self.reactants == r.reactants
            and self.catalysts == r.catalysts
            and self.data == r.data
            and self.reversible == r.reversible
        )

    def __hash__(self):
        """Return a hashcode for self."""
        t = tuple(self.species())
        return hash(t)

    def __repr__(self):
        """Return a debugging string representation of self."""
        return "Reaction(%r, %r, %r, %r)" % (
            self.reactants,
            self.catalysts,
            self.reversible,
            self.data,
        )

    def __str__(self):
        """Return a string representation of self."""
        substrates = ""
        products = ""
        all_species = sorted(self.reactants)
        for species in all_species:
            stoch = self.reactants[species]
            if stoch < 0:
                # species is a substrate:
                if substrates != "":
                    substrates = substrates + " + "
                if stoch != -1:
                    substrates = substrates + str(abs(stoch)) + " "
                substrates = substrates + str(species)
            elif stoch > 0:
                # species is a product:
                if products != "":
                    products = products + " + "
                if stoch != 1:
                    products = products + str(stoch) + " "
                products = products + str(species)
            else:
                raise AttributeError("Invalid 0 coefficient in Reaction.reactants")
        if self.reversible:
            return substrates + " <=> " + products
        else:
            return substrates + " --> " + products

    def reverse(self):
        """Return a new Reaction that is the reverse of self."""
        reactants = {}
        for r in self.reactants:
            reactants[r] = -self.reactants[r]
        return Reaction(reactants, self.catalysts, self.reversible, self.data)

    def species(self):
        """Return a list of all Species involved in self."""
        return list(self.reactants)


class System:
    """Abstraction for a collection of reactions.

    This class is used in the Bio.Pathway framework to represent an arbitrary
    collection of reactions without explicitly defined links.

    Attributes:
     - None

    """

    def __init__(self, reactions=()):
        """Initialize a new System object."""
        self.__reactions = set(reactions)

    def __repr__(self):
        """Return a debugging string representation of self."""
        return "System(" + ",".join(map(repr, self.__reactions)) + ")"

    def __str__(self):
        """Return a string representation of self."""
        return (
            "System of "
            + str(len(self.__reactions))
            + " reactions involving "
            + str(len(self.species()))
            + " species"
        )

    def add_reaction(self, reaction):
        """Add reaction to self."""
        self.__reactions.add(reaction)

    def remove_reaction(self, reaction):
        """Remove reaction from self."""
        self.__reactions.remove(reaction)

    def reactions(self):
        """Return a list of the reactions in this system.

        Note the order is arbitrary!
        """
        # TODO - Define __lt__ so that Reactions can be sorted on Python?
        return list(self.__reactions)

    def species(self):
        """Return a list of the species in this system."""
        return sorted(
            set(reduce(lambda s, x: s + x, [x.species() for x in self.reactions()], []))
        )

    def stochiometry(self):
        """Compute the stoichiometry matrix for self.

        Returns (species, reactions, stoch) where:
         - species    = ordered list of species in this system
         - reactions  = ordered list of reactions in this system
         - stoch      = 2D array where stoch[i][j] is coef of the
           jth species in the ith reaction, as defined
           by species and reactions above

        """
        # Note: This an inefficient and ugly temporary implementation.
        #       To be practical, stochiometric matrices should probably
        #       be implemented by sparse matrices, which would require
        #       NumPy dependencies.
        #
        # PS: We should implement automatic checking for NumPy here.
        species = self.species()
        reactions = self.reactions()
        stoch = [] * len(reactions)
        for i in range(len(reactions)):
            stoch[i] = 0 * len(species)
            for s in reactions[i].species():
                stoch[species.index(s)] = reactions[i].reactants[s]
        return (species, reactions, stoch)


class Interaction:
    """An arbitrary interaction between any number of species.

    This class definition is intended solely as a minimal wrapper interface that should
    be implemented and extended by more specific abstractions.

    Attributes:
     - data      -- reference to arbitrary additional data

    """

    def __init_(self, data):
        self.data = data

    def __hash__(self):
        """Return a hashcode for self."""
        return hash(self.data)

    def __repr__(self):
        """Return a debugging string representation of self."""
        return "Interaction(" + repr(self.data) + ")"

    def __str__(self):
        """Return a string representation of self."""
        return "<" + str(self.data) + ">"


class Network:
    """A set of species that are explicitly linked by interactions.

    The network is a directed multigraph with labeled edges. The nodes in the graph
    are the biochemical species involved. The edges represent an interaction between
    two species, and the edge label is a reference to the associated Interaction
    object.

    Attributes:
     - None

    """

    def __init__(self, species=()):
        """Initialize a new Network object."""
        self.__graph = MultiGraph(species)

    def __repr__(self):
        """Return a debugging string representation of this network."""
        return "<Network: __graph: " + repr(self.__graph) + ">"

    def __str__(self):
        """Return a string representation of this network."""
        return "Network of %i species and %i interactions." % (
            len(self.species()),
            len(self.interactions()),
        )

    def add_species(self, species):
        """Add species to this network."""
        self.__graph.add_node(species)

    def add_interaction(self, source, sink, interaction):
        """Add interaction to this network."""
        self.__graph.add_edge(source, sink, interaction)

    def source(self, species):
        """Return list of unique sources for species."""
        return self.__graph.parents(species)

    def source_interactions(self, species):
        """Return list of (source, interaction) pairs for species."""
        return self.__graph.parent_edges(species)

    def sink(self, species):
        """Return list of unique sinks for species."""
        return self.__graph.children(species)

    def sink_interactions(self, species):
        """Return list of (sink, interaction) pairs for species."""
        return self.__graph.child_edges(species)

    def species(self):
        """Return list of the species in this network."""
        return self.__graph.nodes()

    def interactions(self):
        """Return list of the unique interactions in this network."""
        return self.__graph.labels()
