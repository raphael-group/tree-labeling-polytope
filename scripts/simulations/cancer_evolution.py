import argparse
import itertools

import pandas as pd
import numpy as np
import networkx as nx
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict
from dataclasses import dataclass
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

def phenotype(cell, mutation_type_map):
    return frozenset(mutation for mutation in cell.mutations if mutation_type_map[mutation] == MutationType.DRIVER)

def compute_birth_probability(cell, base_growth_prob, driver_fitness, passenger_fitness, carrying_capacity, mutation_type_map, N_cs):
    N_c = N_cs[(cell.anatomical_site, phenotype(cell, mutation_type_map))]
    K_c = carrying_capacity * len(phenotype(cell, mutation_type_map))
    log_p = np.log(base_growth_prob)
    for mutation in cell.mutations:
        mut_fitness = driver_fitness if mutation_type_map[mutation] == MutationType.DRIVER else passenger_fitness
        log_p += np.log((1+mut_fitness)*(1 - (N_c/K_c)))

    return np.exp(log_p)

"""
Simulate the metastatic migration of cells in a generation.
Taking a generation of cells as a list, we simulate the 
migration of cells to other anatomical sites, outputting 
the a list of migrated cells of the same length and phenotypes,
but with different anatomical sites.
"""
def simulate_migration(T, generation, mutation_type_map, migration_rate, structure, mean_migrations=3.0):
    anatomical_sites = set(cell.anatomical_site for cell in generation)
    all_anatomical_sites = set(cell.anatomical_site for cell in T.nodes())

    migration_graph = nx.DiGraph()
    for site in all_anatomical_sites:
        migration_graph.add_node(site)

    for (u, v) in T.edges():
        if u.anatomical_site == v.anatomical_site:
            continue

        if migration_graph.has_edge(u.anatomical_site, v.anatomical_site):
            continue

        migration_graph.add_edge(u.anatomical_site, v.anatomical_site)

    N_cs = defaultdict(int) # map from (site, phenotype) to number of cells
    for cell in generation:
        N_cs[(cell.anatomical_site, phenotype(cell, mutation_type_map))] += 1

    migration_probs = {}
    for site in anatomical_sites:
        phenotypes = frozenset(
            phenotype(cell, mutation_type_map) 
            for cell in generation if cell.anatomical_site == site and len(phenotype(cell, mutation_type_map)) > 0
        )

        log_migration_prob = np.log(migration_rate)
        for phen in phenotypes:
            N = N_cs[(site, phen)]
            log_migration_prob += np.log(N)
            log_migration_prob += np.log(len(phen))
        migration_probs[site] = np.exp(log_migration_prob)

    def generate_new_anatomical_site(parent_site):
        new_anatomical_site = max(all_anatomical_sites) + 1
        migration_graph.add_node(new_anatomical_site)
        migration_graph.add_edge(parent_site, new_anatomical_site)
        all_anatomical_sites.add(new_anatomical_site)
        return new_anatomical_site

    for site in anatomical_sites:
        if np.random.rand() > migration_probs[site]:
            continue

        cells_to_migrate = [(idx, cell) for idx, cell in enumerate(generation) if cell.anatomical_site == site]
        if len(cells_to_migrate) == 0:
            continue

        if structure == "polyclonal_tree":
            k = np.random.poisson(mean_migrations)
            migrating_cell_indices = np.random.choice([x for x, _ in cells_to_migrate], k, replace=True)
            migrating_cell_indices = np.unique(migrating_cell_indices)
            for idx in migrating_cell_indices:
                migrating_cell = generation[idx]

                if np.random.rand() > 0.5:
                    new_site = generate_new_anatomical_site(site)
                    generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
                else:
                    if len(list(migration_graph.successors(site))) == 0:
                        continue

                    new_site = np.random.choice(list(migration_graph.successors(site)))
                    generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        elif structure == "polyclonal_dag":
            k = np.random.poisson(mean_migrations)
            migrating_cell_indices = np.random.choice([x for x, _ in cells_to_migrate], k, replace=True)
            migrating_cell_indices = np.unique(migrating_cell_indices)
            for idx in migrating_cell_indices:
                migrating_cell = generation[idx]

                if np.random.rand() > 0.5:
                    new_site = generate_new_anatomical_site(site)
                    generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
                else:
                    if len(list(migration_graph.successors(site))) == 0:
                        continue

                    potential_sites = set(migration_graph.nodes()) - set(migration_graph.ancestors(site))
                    new_site = np.random.choice(list(potential_sites))
                    generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        elif structure == "monoclonal_tree":
            idx = np.random.choice([x for x, _ in cells_to_migrate])
            migrating_cell = generation[idx]
            new_site = generate_new_anatomical_site(site)
            generation[idx] = Cell(new_site, migrating_cell.identifier, migrating_cell.mutations)
        elif structure == "monoclonal_dag":
            raise NotImplementedError

            k = np.random.poisson(mean_migrations)
            migrating_cell_indices = np.random.choice(len(cells_to_migrate), k, replace=True)
            migrating_cell_indices = np.unique(migrating_cell_indices)

    return generation

def simulate_evolution(args):
    T = nx.DiGraph()
    founder = Cell(0, 0, frozenset())
    T.add_node(founder)

    generations = {}
    generations[0] = [founder]
    mutation_type_map = {}
    num_mutations = 0
    num_cells = 1
    for i in range(args.generations):
        current_generation = generations[i]
        anatomical_site_counts = defaultdict(int)
        for cell in current_generation:
            anatomical_site_counts[cell.anatomical_site] += 1

        logger.info(f"Generation {i} has {len(current_generation)} cells.")
        logger.info(f"Generation {i} has {len(anatomical_site_counts.keys())} anatomical sites.")

        anatomical_df = pd.DataFrame(anatomical_site_counts.items(), columns=["Anatomical site", "Number of cells"])
        anatomical_df['Fraction of cells'] = anatomical_df['Number of cells'] / len(current_generation)
        # print(anatomical_df.sort_values(by="Number of cells", ascending=False).head(10))

        N_cs = defaultdict(int) # map from (site, phenotype) to number of cells
        for cell in current_generation:
            N_cs[(cell.anatomical_site, phenotype(cell, mutation_type_map))] += 1

        next_generation = []
        for cell in current_generation:
            base_growth_prob = 1/2 + (1/(2 + i))

            p = compute_birth_probability(
                cell, base_growth_prob, args.driver_fitness, args.passenger_fitness, 
                args.carrying_capacity, mutation_type_map, N_cs
            )

            if np.random.rand() > p: # cell dies and has no offspring
                continue

            dc1_id = num_cells 
            dc2_id = num_cells + 1
            num_cells += 2

            daughter_cell1 = Cell(cell.anatomical_site, dc1_id, cell.mutations)
            
            if np.random.rand() < args.mutation_rate:
                mutation = num_mutations + 1
                num_mutations += 1

                daughter_cell2 = Cell(cell.anatomical_site, dc2_id, cell.mutations | frozenset([mutation]))
                if np.random.rand() < args.driver_prob:
                    mutation_type_map[mutation] = MutationType.DRIVER
                else:
                    mutation_type_map[mutation] = MutationType.PASSENGER
            else:
                daughter_cell2 = Cell(cell.anatomical_site, dc2_id, cell.mutations)

            ids = [node.identifier for node in T.nodes]

            next_generation.append(daughter_cell1)
            next_generation.append(daughter_cell2)

            T.add_node(daughter_cell1)
            T.add_node(daughter_cell2)

            T.add_edge(cell, daughter_cell1)
            T.add_edge(cell, daughter_cell2)

        ids = [node.identifier for node in T.nodes]

        metastasized_next_generation = simulate_migration(
            T, next_generation.copy(), mutation_type_map, args.migration_rate, args.structure
        )

        # replace cells with metastasized cells in tree
        for j, cell in enumerate(next_generation):
            T_children = list(T.successors(cell))
            T_parent = list(T.predecessors(cell))[0]

            T.remove_node(cell)
            T.add_node(metastasized_next_generation[j])
            T.add_edge(T_parent, metastasized_next_generation[j])
            for child in T_children:
                T.add_edge(metastasized_next_generation[j], child)

        ids = [node.identifier for node in T.nodes]

        generations[i+1] = metastasized_next_generation

    logger.info("Constructing induced subtree...")
    # select "induced" subtree obtained by sampling n leaves
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

    return T_induced

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate metastatic cancer evolution along a phylogenetic tree.")
    parser.add_argument("-r", "--random-seed", help="Random seed", type=int, default=0)
    parser.add_argument("-o", "--output", help="Output prefix", default="result")
    parser.add_argument("-n", help="Number of leaves to sample", type=int, default=200)
    parser.add_argument("--generations", help="Number of generations", type=int, default=10)
    parser.add_argument("--driver-prob", help="Driver mutation probability", type=float, default=1e-7)
    parser.add_argument("--driver-fitness", help="Driver mutation fitness effect", type=float, default=0.1)
    parser.add_argument("--passenger-fitness", help="Passenger mutation fitness effect", type=float, default=0)
    parser.add_argument("--carrying-capacity", help="Carrying capacity", type=int, default=5e4)
    parser.add_argument("--mutation-rate", help="Mutation rate", type=float, default=0.1)
    parser.add_argument("--migration-rate", help="Migration rate", type=float, default=1e-5)
    parser.add_argument(
        "-s", "--structure", help="Migration graph structure",
        choices=["polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag"],
        default="polyclonal_dag"
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    np.random.seed(args.random_seed)
    T = simulate_evolution(args)

    migration_graph = nx.DiGraph()
    for (u, v) in T.edges():
        if u.anatomical_site == v.anatomical_site:
            continue

        migration_graph.add_edge(u.anatomical_site, v.anatomical_site)

    with open(f"{args.output}_tree_edgelist.tsv", "w") as f:
        for (i, j) in T.edges:
            f.write(f"s{i.identifier}\ts{j.identifier}\n")

    with open(f"{args.output}_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for cell in T.nodes:
            f.write(f"s{cell.identifier},{cell.anatomical_site}\n")

    with open(f"{args.output}_leaf_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for cell in T.nodes:
            if len(list(T.successors(cell))) == 0:
                f.write(f"s{cell.identifier},{cell.anatomical_site}\n")

    with open(f"{args.output}_migration_graph.csv", "w") as f:
        f.write("src,dst\n")
        for (i, j) in migration_graph.edges:
            f.write(f"{i},{j}\n")
    
