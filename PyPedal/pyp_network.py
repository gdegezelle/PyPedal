#!/usr/bin/env python3

###############################################################################
# NAME: pyp_network.py
# VERSION: see PyPedal.__version__
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
# FUNCTIONS:
#   ped_to_graph()
#   find_ancestors()
#   find_descendants()
#   immediate_family()
#   count_offspring()
#   offspring_influence()
#   most_influential_offspring()
#   get_founder_descendants()
#   ---------------------------------------------------------------------------
#   get_node_degrees()
#   get_node_degree_histograms()
#   mean_geodesic()
#   graph_density()
#   dyad_census()  # not a PyPedal 4.0 API; NetworkX 3 has no equivalent
#   mean_degree_centrality()
#   mean_value()
#   get_closeness_centrality()
#   get_clustering_coefficient()
#   get_betweenness_centrality()
#   get_node_betweenness()
###############################################################################

import copy
import logging
import networkx as nx

from .pyp_errors import PyPedalUsageError

logging.basicConfig(level=logging.INFO)

def ped_to_graph(pedobj, oid=False):
    """
    Converts a PyPedal pedigree object to a NetworkX DiGraph.
    """
    di_graph = nx.DiGraph(name=pedobj.kw['pedname'])
    
    missing = str(pedobj.kw['missing_parent'])

    def _node_id(animal_id, original=False):
        if original:
            return int(pedobj.pedigree[int(animal_id) - 1].originalID)
        return int(animal_id)

    for individual in pedobj.pedigree:
        animal_node = int(individual.originalID) if oid else int(individual.animalID)
        sire_attr = missing
        dam_attr = missing
        if str(individual.sireID) != missing:
            sire_attr = str(_node_id(individual.sireID, oid))
        if str(individual.damID) != missing:
            dam_attr = str(_node_id(individual.damID, oid))
        di_graph.add_node(animal_node, sire=sire_attr, dam=dam_attr)

        if str(individual.sireID) != missing:
            di_graph.add_edge(_node_id(individual.sireID, oid), animal_node, sex='s')
        if str(individual.damID) != missing:
            di_graph.add_edge(_node_id(individual.damID, oid), animal_node, sex='d')
    return di_graph


def find_ancestors(pedgraph, anid, _ancestors=None):
    """
    Identifies the ancestors of an animal and returns them in a list.

    The unused ``_ancestors`` argument is accepted for compatibility with
    older call sites that passed an accumulator list.
    """
    try:
        found = list(nx.ancestors(pedgraph, anid))
    except nx.NetworkXError:
        found = []
    if _ancestors is None:
        return found
    if isinstance(_ancestors, dict):
        for ancestor in found:
            _ancestors.setdefault(ancestor, 1)
        return _ancestors
    _ancestors.extend(found)
    return _ancestors


def find_ancestors_g(pedgraph, anid, _ancestors=None, gens=3):
    """
    Identifies the ancestors of an animal up to a specific number of generations.

    Historical callers pass the accumulator as the third argument:
    ``find_ancestors_g(graph, animal_id, {}, gens)``. The swapped form
    ``find_ancestors_g(graph, animal_id, gens, acc)`` is also accepted.
    """
    if isinstance(_ancestors, int) and not isinstance(gens, int):
        _ancestors, gens = gens, _ancestors
    elif isinstance(_ancestors, int) and isinstance(gens, int):
        gens = _ancestors
        _ancestors = {}

    if _ancestors is None:
        _ancestors = {}
    try:
        gens = int(gens)
    except (TypeError, ValueError):
        gens = 3
    if gens <= 0:
        return _ancestors

    try:
        for parent in pedgraph.predecessors(anid):
            if isinstance(_ancestors, dict):
                if parent not in _ancestors:
                    _ancestors[parent] = gens
                    find_ancestors_g(pedgraph, parent, _ancestors, gens - 1)
            elif parent not in _ancestors:
                _ancestors.append(parent)
                find_ancestors_g(pedgraph, parent, _ancestors, gens - 1)
    except nx.NetworkXError:
        pass
    return _ancestors


def find_descendants(pedgraph, anid, _descendants=None):
    """
    Identifies the descendants of an animal and returns them in a list.

    The unused ``_descendants`` argument is accepted for compatibility with
    older call sites that passed an accumulator list.
    """
    try:
        found = list(nx.descendants(pedgraph, anid))
    except nx.NetworkXError:
        found = []
    if _descendants is None:
        return found
    if isinstance(_descendants, dict):
        for descendant in found:
            _descendants.setdefault(descendant, 1)
        return _descendants
    _descendants.extend(found)
    return _descendants


def immediate_family(pedgraph, anid):
    """
    Returns parents and offspring of an animal.
    """
    family = set(pedgraph.predecessors(anid)).union(pedgraph.successors(anid))
    return list(family)


def count_offspring(pedgraph, anid):
    """
    Returns the number of offspring of an animal.
    """
    return len(list(pedgraph.successors(anid)))


def offspring_influence(pedgraph, anid):
    """
    Returns the number of grandchildren by each child of an animal.
    """
    influence = {}
    try:
        for child in pedgraph.successors(anid):
            influence[child] = len(list(pedgraph.successors(child)))
    except nx.NetworkXError:
        pass
    return influence


def most_influential_offspring(pedgraph, anid, resolve='all'):
    """
    Returns the most influential offspring of an animal by their number of offspring.
    """
    offspring_dict = offspring_influence(pedgraph, anid)
    max_value = max(offspring_dict.values(), default=0)

    if resolve == 'all':
        return {k: v for k, v in offspring_dict.items() if v == max_value}
    elif resolve == 'first':
        return next(((k, v) for k, v in offspring_dict.items() if v == max_value), None)
    elif resolve == 'last':
        return next(((k, v) for k, v in reversed(list(offspring_dict.items())) if v == max_value), None)


def get_founder_descendants(pedgraph):
    """
    Returns a dictionary containing the descendants of each founder in the graph.
    """
    founders = {node for node in pedgraph if pedgraph.in_degree(node) == 0}
    return {founder: list(nx.descendants(pedgraph, founder)) for founder in founders}


def get_node_degrees(pg):
    """
    Returns the in, out, and total degrees for nodes in the graph.
    """
    return {
        "in": dict(pg.in_degree),
        "out": dict(pg.out_degree),
        "all": dict(pg.degree)
    }


def get_node_degree_histograms(node_degrees):
    """
    Returns histograms of the in, out, and total degrees for nodes in the graph.
    """
    return {k: {d: list(v).count(d) for d in set(v)} for k, v in node_degrees.items()}


def mean_geodesic(pg):
    """
    Calculates the mean geodesic distance between nodes in the graph.
    """
    try:
        lengths = dict(nx.all_pairs_shortest_path_length(pg))
        total_length = sum(sum(length.values()) for length in lengths.values())
        num_pairs = sum(len(length) for length in lengths.values())
        return total_length / num_pairs
    except ZeroDivisionError:
        return -999


def graph_density(pg):
    """
    Calculates the density of the graph.
    """
    return nx.density(pg)


def dyad_census(pg):
    """
    Historical 2.0.4 wrapper around NetworkX directed dyad_census.

    NetworkX 3 does not provide that function. This is not a PyPedal 4.0
    analysis API.
    """
    raise PyPedalUsageError(
        "pyp_network.dyad_census is outside the PyPedal 4.0 domain. "
        "NetworkX 3 has no directed dyad_census equivalent; the 2.0.4 "
        "local pair census is not restored."
    )


def mean_degree_centrality(pg, normalize=False):
    """
    Calculates mean in- and out-degree centralities for directed graphs.
    """
    in_centrality = sum(dict(pg.in_degree).values())
    out_centrality = sum(dict(pg.out_degree).values())
    num_nodes = pg.number_of_nodes()

    if normalize:
        total_edges = pg.size()
        in_centrality /= total_edges
        out_centrality /= total_edges

    return {"in": in_centrality / num_nodes, "out": out_centrality / num_nodes}


def mean_value(mydict):
    """
    Calculates the mean from all values in a dictionary.
    """
    try:
        return sum(mydict.values()) / len(mydict)
    except ZeroDivisionError:
        return -999


def get_closeness_centrality(pg):
    """
    Returns the closeness centrality for each node in the graph.
    """
    return nx.closeness_centrality(pg)


def get_clustering_coefficient(pg):
    """
    Returns the clustering coefficient for each node in the graph.
    """
    return nx.clustering(pg)


def get_betweenness_centrality(pg):
    """
    Returns the betweenness centrality for each node in the graph.
    """
    return nx.betweenness_centrality(pg)


def get_node_betweenness(pg):
    """
    Returns the betweenness centrality for each node in the graph.
    """
    return nx.betweenness_centrality(pg)
