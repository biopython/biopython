import numpy as np
from Bio.PDB.Entity import Entity


def to_networkx(neighbor_pairs):
    """Convert Biopython neighbor search results to a NetworkX Graph.

    Arguments:
    - neighbor_pairs: List of tuples (Entity, Entity) from NeighborSearch.search_all
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError(
            "NetworkX is required for this function. Install it with 'pip install networkx'."
        )

    G = nx.Graph()

    for entity1, entity2 in neighbor_pairs:
        # Biopython '-' operator calculates Euclidean distance between Entities
        distance = entity1 - entity2

        # Add edge with distance as weight
        G.add_edge(entity1, entity2, weight=float(distance))

        # Add node attributes for richer analysis
        for entity in (entity1, entity2):
            if entity not in G.nodes:
                G.nodes[entity]["level"] = entity.get_level()
                G.nodes[entity]["full_id"] = entity.get_full_id()

    return G
