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

def compute_birth_probability(cell, base_growth_prob, driver_fitness, passenger_fitness, carrying_capacity, mutation_type_map, generation):
    N_c = 0
    for other_cell in generation:
        if phenotype(cell, mutation_type_map) != phenotype(other_cell, mutation_type_map):
            continue

        if other_cell.anatomical_site != cell.anatomical_site:
            continue

        N_c += 1

    K_c = carrying_capacity * len(phenotype(cell, mutation_type_map))
    log_p = np.log(base_growth_prob)
    for mutation in cell.mutations:
        mut_fitness = driver_fitness if mutation_type_map[mutation] == MutationType.DRIVER else passenger_fitness
        log_p += np.log((1+mut_fitness)*(1 - (N_c/K_c)))

    return np.exp(log_p)

def simulate_polyclonal_tree_migration(generation, migration_rate):
    pass

def simulate_polyclonal_dag_migration(generation, migration_rate):
    pass

def simulate_monoclonal_tree_migration(generation, migration_rate):
    pass

def simulate_monoclonal_dag_migration(generation, migration_rate):
    pass

def simulate_migration(generation, migration_rate, structure):
    anatomical_sites = set(cell.anatomical_site for cell in generation)
    for site in anatomical_sites:
        phenotypes = frozenset(phenotype(cell, mutation_type_map) for cell in generation if cell.anatomical_site == site)
        

    if structure == "polyclonal_tree":
        return simulate_polyclonal_tree_migration(generation, migration_rate)
    elif structure == "polyclonal_dag":
        return simulate_polyclonal_dag_migration(generation, migration_rate)
    elif structure == "monoclonal_tree":
        return simulate_monoclonal_tree_migration(generation, migration_rate)
    elif structure == "monoclonal_dag":
        return simulate_monoclonal_dag_migration(generation, migration_rate)
    else:
        raise ValueError(f"Unknown structure {structure}")

def simulate_evolution(args):
    T = nx.DiGraph()
    founder = Cell(0, frozenset())
    T.add_node(founder)

    generations = {}
    generations[0] = [founder]
    mutation_type_map = {}
    num_mutations = 0
    for i in range(args.generations):
        logger.info(f"Creating generation {i}")

        current_generation = generations[i]
        next_generation = []
        for cell in current_generation:
            base_growth_prob = 1/2 + (1/(2 + i))
            p = compute_birth_probability(
                    cell, base_growth_prob, args.driver_fitness, args.passenger_fitness, 
                    args.carrying_capacity, mutation_type_map, current_generation
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

        next_generation = simulate_migration(next_generation, args.migration_rate, args.structure)
        generations[i+1] = next_generation

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate metastatic cancer evolution along a phylogenetic tree.")
    parser.add_argument("-m", help="Number of labels", type=int, default=6)
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
    character_set = range(args.m)
    simulate_evolution(args)

