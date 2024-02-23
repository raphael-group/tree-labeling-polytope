import argparse
import itertools
import math

import pandas as pd
import numpy as np
import networkx as nx
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict
from dataclasses import dataclass, replace
from typing import FrozenSet

@dataclass(frozen=True)
class Cell:
    anatomical_site: int
    identifier: int
    mutations: FrozenSet[int]

@dataclass(frozen=True)
class MutationType(Enum):
    DRIVER = 1
    PASSENGER = 2

@dataclass(frozen=True)
class EvolutionParameters:
    driver_prob: float
    driver_fitness: float
    passenger_fitness: float
    carrying_capacity: int
    mutation_rate: float
    migration_rate: float
    mean_migrations: float
    structure: str

""" 
Global variables to generate unique identifiers for cells and mutations,
and to keep track of the mutation types of each mutation.
"""
cell_id = 0
def generate_cell_id():
    global cell_id
    cell_id += 1
    return cell_id

mutation_id = 0
def generate_mutation_id():
    global mutation_id
    mutation_id += 1
    return mutation_id

mutation_type_map = {}

"""
Given a cell and a map from mutations to mutation types,
return the phenotype of the cell, which is the set of
driver mutations in the cell.
"""
def phenotype(cell):
    return frozenset(mutation for mutation in cell.mutations if mutation_type_map[mutation] == MutationType.DRIVER)

"""
Computes the probability that a cell divides and produces offspring.
"""
def compute_birth_probability(cell, base_growth_prob, driver_fitness, passenger_fitness, carrying_capacity, N_cs):
    N_c = N_cs[(cell.anatomical_site)]
    K_c = carrying_capacity # * len(phenotype(cell))
    log_p = np.log(base_growth_prob)
    for mutation in cell.mutations:
        mut_fitness = driver_fitness if mutation_type_map[mutation] == MutationType.DRIVER else passenger_fitness
        log_p += np.log((1+mut_fitness)*(1 - (N_c/K_c)))

    prob = np.exp(log_p)
    return prob

"""
Given a set of cells to migrate, a mean number of migrations,
a migration graph, and a generation of cells, the following set of functions
simulate the migration of cells to other anatomical sites, ensuring that
the migration graph topology is respected.
"""

""" Simulate the migration of cells assuming the migration graph is a polyclonal tree. """
def polyclonal_tree_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_site_fn):
    k = np.random.poisson(mean_migrations)
    migrating_cell_indices = np.random.choice(cells_to_migrate, k, replace=True)
    migrating_cell_indices = np.unique(migrating_cell_indices)
    for idx in migrating_cell_indices:
        migrating_cell = generation[idx]

        if np.random.rand() > 0.5:
            new_site = new_site_fn(site)
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        else:
            if len(list(migration_graph.successors(site))) == 0:
                continue

            new_site = np.random.choice(list(migration_graph.successors(site)))
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)

    return

""" Simulate the migration of cells assuming the migration graph is a polyclonal DAG. """
def polyclonal_dag_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_site_fn):
    k = np.random.poisson(mean_migrations)
    migrating_cell_indices = np.random.choice(cells_to_migrate, k, replace=True)
    migrating_cell_indices = np.unique(migrating_cell_indices)
    for idx in migrating_cell_indices:
        migrating_cell = generation[idx]

        if np.random.rand() > 0.5:
            new_site = new_site_fn(site)
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        else:
            potential_sites = set(migration_graph.nodes()) - set(nx.ancestors(migration_graph, site))
            if len(potential_sites) == 0:
                continue

            new_site = np.random.choice(list(potential_sites))
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)

    return

""" Simulate the migration of cells assuming the migration graph is a monoclonal tree. """
def monoclonal_tree_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_site_fn):
    idx = np.random.choice(cells_to_migrate)
    migrating_cell = generation[idx]
    new_site = new_site_fn(site)
    generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
    return

""" Simulate the migration of cells assuming the migration graph is a monoclonal DAG. """
def monoclonal_dag_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_site_fn):
    idx = np.random.choice(cells_to_migrate)
    migrating_cell = generation[idx]

    if np.random.rand() > 0.5:
        new_site = new_site_fn(site)
        generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        return
    else:
        potential_sites = set(migration_graph.nodes()) - set(nx.ancestors(migration_graph, site))
        if len(potential_sites) == 0:
            return

        new_site = np.random.choice(list(potential_sites))
        generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)

    return

""" Simulate the migration of cells assuming the migration graph is arbitrary. """
def arbitrary_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_site_fn):
    k = np.random.poisson(mean_migrations)
    migrating_cell_indices = np.random.choice(cells_to_migrate, k, replace=True)
    migrating_cell_indices = np.unique(migrating_cell_indices)
    for idx in migrating_cell_indices:
        migrating_cell = generation[idx]

        r = np.random.rand()
        if r > 0.5:
            new_site = new_site_fn(site)
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        elif r < 0.1:
            potential_sites = set(nx.ancestors(migration_graph, site))
            if len(potential_sites) == 0:
                continue

            new_site = np.random.choice(list(potential_sites))
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        else:
            potential_sites = set(migration_graph.nodes())
            if len(potential_sites) == 0:
                continue

            new_site = np.random.choice(list(potential_sites))
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)

    return

"""
Simulate the metastatic migration of cells in a generation.
Taking a generation of cells as a list, we simulate the 
migration of cells to other anatomical sites, outputting 
the a list of migrated cells of the same length and phenotypes,
but with migrated anatomical sites.
"""
def simulate_migration(T, generation, params):
    generation = generation.copy()
    
    N_cs = defaultdict(int) # map from (site, phenotype) to number of cells
    for cell in generation:
        N_cs[(cell.anatomical_site, phenotype(cell))] += 1

    anatomical_sites = set(cell.anatomical_site for cell in generation)
    all_anatomical_sites = set(cell.anatomical_site for cell in T.nodes())

    # construct migration graph
    migration_graph = nx.DiGraph()
    for site in all_anatomical_sites:
        migration_graph.add_node(site)

    for (u, v) in T.edges():
        if u.anatomical_site == v.anatomical_site:
            continue

        if migration_graph.has_edge(u.anatomical_site, v.anatomical_site):
            continue

        migration_graph.add_edge(u.anatomical_site, v.anatomical_site)

    # construct map from site to cells with that site
    sites_to_cell_idxs = defaultdict(list)
    site_to_phenotypes = defaultdict(set)
    for idx, cell in enumerate(generation):
        sites_to_cell_idxs[cell.anatomical_site].append(idx)
        if len(phenotype(cell)) > 0:
            site_to_phenotypes[cell.anatomical_site].add(phenotype(cell))

    # compute probability of migration *from* each site
    migration_probs = {}
    for site in anatomical_sites:
        migration_probs[site] = params.migration_rate * len(sites_to_cell_idxs[site])

    def new_anatomical_site_fn(parent_site):
        new_anatomical_site = max(all_anatomical_sites) + 1
        migration_graph.add_node(new_anatomical_site)
        migration_graph.add_edge(parent_site, new_anatomical_site)
        all_anatomical_sites.add(new_anatomical_site)
        return new_anatomical_site

    for site in anatomical_sites:
        if np.random.rand() > migration_probs[site]:
            continue

        cells_to_migrate = sites_to_cell_idxs[site]
        if len(cells_to_migrate) == 0:
            continue

        mean_migrations = params.mean_migrations

        if params.structure == "polyclonal_tree":
            polyclonal_tree_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_anatomical_site_fn)
        elif params.structure == "polyclonal_dag":
            polyclonal_dag_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_anatomical_site_fn)
        elif params.structure == "monoclonal_tree":
            monoclonal_tree_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_anatomical_site_fn)
        elif params.structure == "monoclonal_dag":
            monoclonal_dag_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_anatomical_site_fn)
        else:
            arbitrary_migration(mean_migrations, cells_to_migrate, site, migration_graph, generation, new_anatomical_site_fn)

    return generation

"""
Simulates division of the current generation of cells,
returning the next generation of cells and the updated
cell lineage tree.
"""
def simulate_cell_division(T, current_generation, base_growth_prob, params):
    T = T.copy()

    # count the number of cells at each anatomical site 
    N_cs = defaultdict(int) 
    for cell in current_generation:
        N_cs[(cell.anatomical_site)] += 1

    next_generation = []
    for cell in current_generation:
        p = compute_birth_probability(
            cell, base_growth_prob, params.driver_fitness, params.passenger_fitness, 
            params.carrying_capacity, N_cs
        )

        if np.random.rand() > p: # cell dies and has no offspring
            continue

        dc1_id = generate_cell_id()
        dc2_id = generate_cell_id()

        daughter_cell1 = Cell(cell.anatomical_site, dc1_id, cell.mutations)
        
        if np.random.rand() < params.mutation_rate:
            mutation = generate_mutation_id()

            daughter_cell2 = Cell(cell.anatomical_site, dc2_id, cell.mutations | frozenset([mutation]))
            if np.random.rand() < params.driver_prob:
                mutation_type_map[mutation] = MutationType.DRIVER
            else:
                mutation_type_map[mutation] = MutationType.PASSENGER
        else:
            daughter_cell2 = Cell(cell.anatomical_site, dc2_id, cell.mutations)

        next_generation.append(daughter_cell1)
        next_generation.append(daughter_cell2)

        T.add_node(daughter_cell1)
        T.add_node(daughter_cell2)

        T.add_edge(cell, daughter_cell1)
        T.add_edge(cell, daughter_cell2)

    return T, next_generation

"""
Simulate the evolution of a cancer lineage tree.
"""
def simulate_evolution(params: EvolutionParameters, n: int, num_generations: int) -> nx.DiGraph:
    T = nx.DiGraph()
    founder = Cell(0, 0, frozenset())
    T.add_node(founder)

    generations = {}
    generations[0] = [founder]
    for i in range(num_generations):
        current_generation = generations[i]
        base_growth_prob = 1/2 + (1/(2 + i))

        # display anatomical site statistics
        anatomical_site_counts = defaultdict(int)
        for cell in current_generation:
            anatomical_site_counts[cell.anatomical_site] += 1

        logger.info(f"Generation {i} has {len(current_generation)} cells.")
        logger.info(f"Generation {i} has {len(anatomical_site_counts.keys())} anatomical sites.")

        anatomical_df = pd.DataFrame(anatomical_site_counts.items(), columns=["Anatomical site", "Number of cells"])
        anatomical_df['Fraction of cells'] = anatomical_df['Number of cells'] / len(current_generation)
        # print(anatomical_df)

        T, next_generation = simulate_cell_division(
            T, current_generation, base_growth_prob, params
        )

        metastasized_next_generation = simulate_migration(
            T, next_generation, params
        )

        """ Replace cells with metastasized cells in cell lineage tree """
        for j, cell in enumerate(next_generation):
            T_children = list(T.successors(cell))
            T_parent = list(T.predecessors(cell))[0]

            T.remove_node(cell)
            T.add_node(metastasized_next_generation[j])
            T.add_edge(T_parent, metastasized_next_generation[j])
            for child in T_children:
                T.add_edge(metastasized_next_generation[j], child)

        generations[i+1] = metastasized_next_generation

    logger.info("Constructing induced subtree...")
    selected_leaves = np.random.choice(generations[args.generations], args.n, replace=False)
    marked_nodes = set(selected_leaves)
    for n in nx.dfs_postorder_nodes(T, source=founder):
        if any(marked in T.successors(n) for marked in marked_nodes):
            marked_nodes.add(n)
    
    logger.info(f"Contracting rank-2 nodes...")
    # TODO: this loops take quadratic time
    T_induced = nx.induced_subgraph(T, marked_nodes).copy()
    while True:
        has_rank_two_nodes = False
        nodes = list(T_induced.nodes)
        for n in nodes:
            if n not in T_induced.nodes: continue
            if T_induced.out_degree(n) != 1: continue
            if T_induced.in_degree(n) != 1: continue
            parent = list(T_induced.predecessors(n))[0]
            child = list(T_induced.successors(n))[0]
            if n.anatomical_site != child.anatomical_site: continue

            has_rank_two_nodes = True
            T_induced.remove_node(n)
            T_induced.add_edge(parent, child)
            break

        if not has_rank_two_nodes:
            break

    anatomical_sites = set(cell.anatomical_site for cell in T_induced.nodes())
    anatomical_sites_in_leaves = set(cell.anatomical_site for cell in selected_leaves)
    logger.info(f"Number of anatomical sites in tree: {len(anatomical_sites)}")
    logger.info(f"Number of anatomical sites in leaves: {len(anatomical_sites_in_leaves)}")

    # replace anatomical sites of cells in induced tree if they are not in the leaves
    for n in list(nx.dfs_preorder_nodes(T_induced, source=founder)):
        if n not in T_induced.nodes: continue
        if n.anatomical_site in anatomical_sites_in_leaves: continue
        if T_induced.in_degree(n) == 0: continue

        parent = list(T_induced.predecessors(n))[0]
        children = list(T_induced.successors(n))

        T.remove_node(n)
        T.add_node(Cell(parent.anatomical_site, n.identifier, n.mutations))
        T.add_edge(parent, n)
        for child in children:
            T.add_edge(n, child)

    return T_induced

"""
Stochastically pertubs T using NNI operations.
"""
def stochastic_spr(T, num_perturbations=10):
    T = T.copy()

    vertices = list(T.nodes())
    count = 0
    while count < num_perturbations:
        v = vertices[np.random.choice(len(vertices))]

        if T.in_degree(v) == 0:
            continue

        u = list(T.predecessors(v))[0]

        candidates = set(T.nodes()) - set(nx.descendants(T, v))
        candidates.remove(v)

        if len(candidates) == 0:
            continue

        w = list(candidates)[np.random.choice(len(candidates))]
        if w == u:
            continue

        T.remove_edge(u, v)
        T.add_edge(w, v)

        count += 1

    return T

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate metastatic cancer evolution along a phylogenetic tree.")
    parser.add_argument("-r", "--random-seed", help="Random seed", type=int, default=0)
    parser.add_argument("-o", "--output", help="Output prefix", default="result")
    parser.add_argument("-n", help="Number of leaves to sample", type=int, default=200)
    parser.add_argument("-e", "--errors", help="Number of errors (default: 0)", type=float, default=0)
    parser.add_argument("--generations", help="Number of generations", type=int, default=40)
    parser.add_argument("--driver-prob", help="Driver mutation probability", type=float, default=2e-7)
    parser.add_argument("--driver-fitness", help="Driver mutation fitness effect", type=float, default=0.1)
    parser.add_argument("--passenger-fitness", help="Passenger mutation fitness effect", type=float, default=0)
    parser.add_argument("--carrying-capacity", help="Carrying capacity", type=int, default=5000)
    parser.add_argument("--mutation-rate", help="Mutation rate", type=float, default=0.1)
    parser.add_argument("--migration-rate", help="Migration rate", type=float, default=1e-15)
    parser.add_argument(
        "-s", "--structure", help="Migration graph structure",
        choices=["polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag", "none"],
        default="polyclonal_dag"
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    params = EvolutionParameters(
        driver_prob=args.driver_prob,
        driver_fitness=args.driver_fitness,
        passenger_fitness=args.passenger_fitness,
        carrying_capacity=args.carrying_capacity,
        mutation_rate=args.mutation_rate,
        migration_rate=args.migration_rate,
        mean_migrations=3.0,
        structure=args.structure
    )

    np.random.seed(args.random_seed)
    T = simulate_evolution(params, args.n, args.generations)

    if args.errors > 0:
        T_perturbed = stochastic_spr(T, args.errors)
    else:
        T_perturbed = T

    migration_graph = nx.DiGraph()
    for (u, v) in T.edges():
        if u.anatomical_site == v.anatomical_site:
            continue

        migration_graph.add_edge(u.anatomical_site, v.anatomical_site)

    with open(f"{args.output}_tree_edgelist.tsv", "w") as f:
        for (i, j) in T.edges:
            f.write(f"s{i.identifier}\ts{j.identifier}\n")

    with open(f"{args.output}_perturbed_tree_edgelist.tsv", "w") as f:
        for (i, j) in T_perturbed.edges:
            f.write(f"s{i.identifier}\ts{j.identifier}\n")

    with open(f"{args.output}_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for cell in T.nodes:
            f.write(f"s{cell.identifier},{cell.anatomical_site}\n")

    with open(f"{args.output}_leaf_labeling.csv", "w") as f:
        f.write("leaf,label\n")
        for cell in T_perturbed.nodes:
            if len(list(T_perturbed.successors(cell))) == 0:
                f.write(f"s{cell.identifier},{cell.anatomical_site}\n")

    with open(f"{args.output}_migration_graph.csv", "w") as f:
        f.write("src,dst\n")
        for (i, j) in migration_graph.edges:
            f.write(f"{i},{j}\n")
