import argparse
import pandas as pd
import gurobipy as gp
import networkx as nx

from collections import defaultdict
from gurobipy import GRB
from tqdm import tqdm

def build_model(character_matrix, missing_character, cell_to_idx):
    model = gp.Model("Large star homoplasy parsimony problem")

    idx_to_cell = {idx: cell for cell, idx in cell_to_idx.items()}

    characters = character_matrix.columns
    character_to_mutations = defaultdict(list)
    character_to_states = {}
    for character in characters:
        character_states = set(character_matrix[character].unique()) - {missing_character, "0"}
        for state in character_states:
            character_to_mutations[character].append((state, state))
            character_to_mutations[character].append(("0", state))
            character_to_mutations[character].append((state, missing_character))
        character_to_mutations[character].append(("0", "0"))
        character_to_mutations[character].append(("0", missing_character))
        character_to_states[character] = set(character_matrix[character].unique()) | {"0", missing_character}

    N = len(character_matrix)
    character_mutations = [(c, s1, s2) for c in characters for s1, s2 in character_to_mutations[c]]
    x = model.addVars(2 * N, 2 * N, character_mutations, vtype=gp.GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(2 * N, 2 * N, vtype=gp.GRB.BINARY, name="y", lb=0, ub=1)
    print(character_mutations)
    print(len(character_mutations))

    internal_nodes = list(range(N, 2 * N))
    leaf_nodes = list(range(N))

    # LEAF CONSTRAINTS
    for character in characters:
        for v in leaf_nodes:
            states = character_to_states[character]
            for state in states:
                parent_states = [s1 for (s1, s2) in character_to_mutations[character] if s2 == state]
                if character_matrix.at[idx_to_cell[v], character] == state:
                    model.addConstr(
                        gp.quicksum(
                            x[w, v, character, parent_state, state] for w in internal_nodes for parent_state in parent_states 
                        ) == 1
                    )
                else:
                    model.addConstr(
                        gp.quicksum(
                            x[w, v, character, parent_state, state] for w in internal_nodes for parent_state in parent_states
                        ) == 0
                    )

    for character in characters:
        for u in internal_nodes:
            for v in range(2*N):
                if u == v: continue
                model.addConstr(gp.quicksum(x[u, v, character, s1, s2] for s1, s2 in character_to_mutations[character]) == y[u, v])

    # INTERNAL NODE CONSTRAINTS
    if False:
        for character in tqdm(characters):
            for u in internal_nodes:
                for v in internal_nodes:
                    for w in range(2*N):
                        if u == v or u == w or v == w:
                            continue
                        for s1 in character_to_states[character]:
                            parent_states = [s2 for (s2, s3) in character_to_mutations[character] if s3 == s1]
                            child_states = [s3 for (s2, s3) in character_to_mutations[character] if s2 == s1]
                            lhs_expr = gp.quicksum(x[u, v, character, s2, s1] for s2 in parent_states)
                            rhs_expr = gp.quicksum(x[v, w, character, s1, s3] for s3 in child_states)
                            c1 = model.addConstr(lhs_expr <= rhs_expr + (1 - y[u, v]) + (1 - y[v, w]))
                            c2 = model.addConstr(lhs_expr >= rhs_expr - (1 - y[u, v]) - (1 - y[v, w]))

    # CONNECTING CONSTRAINTS
    for u in range(2 * N):
        for v in range(2 * N):
            for (c, s1, s2) in character_mutations:
                model.addConstr(x[u,v,c,s1,s2] <= y[u,v])

    # ROOT CONSTRAINTS
    DUMMY_ROOT = 2 * N - 1
    ROOT = 2 * N - 2
    for u in range(2 * N - 2):
        model.addConstr(y[u, ROOT] == 0)
        model.addConstr(y[DUMMY_ROOT, u] == 0)
        model.addConstr(y[u, DUMMY_ROOT] == 0)

    model.addConstr(y[DUMMY_ROOT, ROOT] == 1)

    # LEAF CONSTRAINTS
    for v in range(N):
        for u in range(2 * N):
            model.addConstr(y[v, u] == 0)

    # TREE CONSTRAINTS
    for v in range(2 * N - 1):
        model.addConstr(gp.quicksum(y[u, v] for u in range(2 * N) if u != v) == 1)

    for v in internal_nodes:
        if v == DUMMY_ROOT: continue
        model.addConstr(gp.quicksum(y[v, u] for u in range(2 * N) if u != v) == 2) # make it 2 for binary trees

    # SYMMETRY BREAKING CONSTRAINTS
    for u in range(2 * N):
        for v in range(2 * N):
            if u > v: continue
            model.addConstr(y[u, v] == 0)

    # TREE CONSTRAINTS ARE ADDED AS LAZY CONSTRAINTS
    model.addConstr(gp.quicksum(y[u, v] for u in range(2 * N) for v in range(2 * N) if u != v) == 2 * N - 1)

    # OBJECTIVE FUNCTION
    model.setObjective(
        gp.quicksum(x[u, v, c, s1, s2] for u in range(2 * N - 1) for v in range(2 * N - 1) for c, s1, s2 in character_mutations 
                    if u != v and s1 != s2 and s1 != missing_character and s2 != missing_character), 
        gp.GRB.MINIMIZE
    )

    def lazy_constraint_callback(model, where):
        if where != GRB.Callback.MIPSOL:
            return

        x_sol = model.cbGetSolution(x)
        y_sol = model.cbGetSolution(y)

        labeling = defaultdict(dict)
        violated_constraint = None
        for u in range(2 * N):
            for v in range(2 * N):
                if u == v: continue
                if y_sol[u, v] < 0.5: continue
                if violated_constraint is not None:
                    break

                for c, s1, s2 in character_mutations:
                    if x_sol[u, v, c, s1, s2] < 0.5:
                        continue

                    if labeling[u].get(c, s1) != s1:
                        violated_constraint = (u, c, s1, labeling[u][c])
                        break

                    if labeling[v].get(c, s2) != s2:
                        violated_constraint = (v, c, s2, labeling[v][c])
                        break

                    labeling[u][c] = s1
                    labeling[v][c] = s2

        if not violated_constraint:
            print("No violated constraints found.")
            return

        u, c, s1, s2 = violated_constraint
        print(f"Adding constraint for vertex {u} double labeled with '{s1}' and '{s2}' for character {c}.")
        parents = [v for v in range(2 * N) if y_sol[v, u] > 0.5 and v != u]
        assert len(parents) == 1
        p = parents[0]
        children = [v for v in range(2 * N) if y_sol[u, v] > 0.5 and v != u]

        for state in [s1, s2]:
            for child in children:
                parent_states = [parent_state for (parent_state, child_state) in character_to_mutations[c] if child_state == state]
                child_states = [child_state for (parent_state, child_state) in character_to_mutations[c] if parent_state == state]
                lhs_expr = gp.quicksum(x[p, u, c, parent_state, state] for parent_state in parent_states)
                rhs_expr = gp.quicksum(x[u, child, c, state, child_state] for child_state in child_states)
                model.cbLazy(lhs_expr <= rhs_expr + (1 - y[p, u]) + (1 - y[u, child]))
                model.cbLazy(lhs_expr >= rhs_expr - (1 - y[p, u]) - (1 - y[u, child]))

        return

        assert False

    model.Params.LazyConstraints = 1
    model.optimize(lazy_constraint_callback)

    # print y
    for u in range(2 * N):
        for v in range(2 * N):
            if u != v and y[u, v].x > 0.5:
                print(f"y[{u},{v}] = {y[u, v].x}")

    # print x
    for u in range(2 * N):
        for v in range(2 * N):
            for c, s1, s2 in character_mutations:
                if x[u, v, c, s1, s2].x > 0.5:
                    print(f"x[{u},{v},{c},{s1},{s2}] = {x[u, v, c, s1, s2].x}")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Solves the star homoplasy large parsimony problem using the tree labeling polytope."
    )

    parser.add_argument(
        "character_matrix",
        type=str,
        help="Path to the character matrix file."
    )

    parser.add_argument(
        "-m", "--missing",
        type=str,
        help="Missing data character.",
        default="-"
    )

    parser.add_argument(
        "-s", "--separator",
        type=str,
        help="Character matrix separator.",
        default="\t"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    character_matrix = pd.read_csv(args.character_matrix, sep=args.separator, index_col=0).astype(str)
    missing_character = args.missing
    cell_to_idx = {cell: idx for idx, cell in enumerate(character_matrix.index)}

    model = build_model(character_matrix, missing_character, cell_to_idx)


