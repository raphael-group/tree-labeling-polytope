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

def simulate_migration(T, generation, mutation_type_map, migration_rate, structure, mean_migrations=3.0):
    print(len(generation))
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

        cells_to_migrate = [cell for cell in generation if cell.anatomical_site == site]
        if len(cells_to_migrate) == 0:
            continue

        if structure == "polyclonal_tree":
            k = np.random.poisson(mean_migrations)
            migrating_cell_indices = np.random.choice(len(cells_to_migrate), k, replace=True)
            migrating_cell_indices = np.unique(migrating_cell_indices)
            for idx in migrating_cell_indices:
                if np.random.rand() > 0.5:
                    new_site = generate_new_anatomical_site(site)
                    generation[idx] = Cell(new_site, cells_to_migrate[idx].mutations)
                else:
                    if len(list(migration_graph.successors(site))) == 0:
                        continue

                    new_site = np.random.choice(list(migration_graph.successors(site)))
                    generation[idx] = Cell(new_site, cells_to_migrate[idx].mutations)
        elif structure == "polyclonal_dag":
            k = np.random.poisson(mean_migrations)
            migrating_cell_indices = np.random.choice(len(cells_to_migrate), k, replace=True)
            migrating_cell_indices = np.unique(migrating_cell_indices)
            for idx in migrating_cell_indices:
                if np.random.rand() > 0.5:
                    new_site = generate_new_anatomical_site(site)
                    generation[idx] = Cell(new_site, cells_to_migrate[idx].mutations)
                else:
                    if len(list(migration_graph.successors(site))) == 0:
                        continue

                    potential_sites = set(migration_graph.nodes()) - set(migration_graph.ancestors(site))
                    new_site = np.random.choice(list(potential_sites))
                    generation[idx] = Cell(new_site, cells_to_migrate[idx].mutations)
        elif structure == "monoclonal_tree":
            migrating_cell_idx = np.random.randint(0, len(cells_to_migrate))
            migrating_cell = cells_to_migrate[migrating_cell_idx]
            new_site = generate_new_anatomical_site(site)
            generation[migrating_cell_idx] = Cell(new_site, migrating_cell.mutations)
        elif structure == "monoclonal_dag":
            k = np.random.poisson(mean_migrations)
            pass

    return generation

def simulate_evolution(args):
    T = nx.DiGraph()
    founder = Cell(0, frozenset())
    T.add_node(founder)

    generations = {}
    generations[0] = [founder]
    mutation_type_map = {}
    num_mutations = 0
    for i in range(args.generations):
        current_generation = generations[i]
        anatomical_site_counts = defaultdict(int)
        for cell in current_generation:
            anatomical_site_counts[cell.anatomical_site] += 1

        logger.info(f"Generation {i} has {len(current_generation)} cells.")
        logger.info(f"Generation {i} has {len(anatomical_site_counts.keys())} anatomical sites.")

        anatomical_df = pd.DataFrame(anatomical_site_counts.items(), columns=["Anatomical site", "Number of cells"])
        anatomical_df['Fraction of cells'] = anatomical_df['Number of cells'] / len(current_generation)
        logger.info(anatomical_df.sort_values(by="Number of cells", ascending=False).head(10))

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

            if np.random.rand() > p: # cell dies
                continue

            daughter_cell1 = Cell(cell.anatomical_site, cell.mutations)
            
            if np.random.rand() < args.mutation_rate:
                mutation = num_mutations + 1
                num_mutations += 1

                daughter_cell2 = Cell(cell.anatomical_site, cell.mutations | frozenset([mutation]))
                if np.random.rand() < args.driver_prob:
                    mutation_type_map[mutation] = MutationType.DRIVER
                else:
                    mutation_type_map[mutation] = MutationType.PASSENGER
            else:
                daughter_cell2 = Cell(cell.anatomical_site, cell.mutations)

            next_generation.append(daughter_cell1)
            next_generation.append(daughter_cell2)

            T.add_node(daughter_cell1)
            T.add_node(daughter_cell2)

            T.add_edge(cell, daughter_cell1)
            T.add_edge(cell, daughter_cell2)

        next_generation = simulate_migration(T, next_generation, mutation_type_map, args.migration_rate, args.structure)
        generations[i+1] = next_generation

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate metastatic cancer evolution along a phylogenetic tree.")
    parser.add_argument("-r", "--random-seed", help="Random seed", type=int, default=0)
    parser.add_argument("-o", "--output", help="Output prefix", default="result")
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
    simulate_evolution(args)
