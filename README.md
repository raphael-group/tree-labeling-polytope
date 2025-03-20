## The Tree Labeling Polytope

The tree labeling polytope is a technical tool for fixed topology 
ancestral reconstruction problems. Here, we provide an implementation
of the tree labeling polytope for the following three problems:
- the parsimonious migration history problem
- the softwired small parsimony problem 
- the convex recoloring problem

Each problem takes as input either a tree topology or a phylogenetic
network and a leaf labeling, and outputs a labeling of the internal
nodes of the tree or network that is optimal with respect to the
problem-specific objective function.

### Usage

The tool offers three modes, each corresponding to one of the three
problems. 

Each mode offers a different set of options,
which are described in the help message for each mode. The top-level
help message is shown below for the parsimonious migration
history and convex recoloring problem solvers:
```bash
$ python scripts/tlp.py --help
usage: tlp.py [-h] {fast_machina,convex_recoloring} ...

Constrained tree labeling using the tree labeling polytope.

positional arguments:
  {fast_machina,convex_recoloring}
                        Methods
    fast_machina        fastMACHINA
    convex_recoloring   Convex Recoloring

options:
  -h, --help            show this help message and exit
```

The top-level help message is shown below for the softwired
small parsimony problem solvers:
```bash
python scripts/softwired_tlp.py --help
usage: softwired_tlp.py [-h] [--mode {scornavacca,tlp}] network sequences output

Softwired parsimony problem solver using TLP

positional arguments:
  network               Phylogenetic network in edgelist format
  sequences             Leaf sequences in CSV format
  output                Output prefix

options:
  -h, --help            show this help message and exit
  --mode {scornavacca,tlp}
                        Mode to use
```

The tool requires the following dependencies:
- Python 3.11 or higher
- NetworkX
- NumPy
- loguru
- Pyomo
- Gurobi

We note that a Gurobi license is currently required to use the tool,
but it can be replaced with another solver supported by Pyomo
(e.g., CPLEX, GLPK) upon request. A Gurobi license is free for academic
use and can be obtained from the Gurobi website.

#### Example 1: The Parsimonious Migration History Problem

We demonstrate the use of the tree labeling polytope on tumor clone
`CP28` from a single-cell lineage tracing dataset modeling metastatic 
dissemination in mice from Yang et al. (2021). The tree topology
and branch lengths were inferred by applying Startle followed
by LAML, which are specialized tools for inferring phylogenies
from lineage tracing data. The tree topology and branch lengths
are stored in the file `examples/CP28_tree_edgelist.tsv`. The
leaf labeling is stored in the file `examples/CP28_leaf_labeling.csv`
and describes the anatomical location of each cell in the tree.

We apply the tree labeling polytope to infer the most parsimonious 
ancestral labeling such that the induced migration graph is
*i)* unconstrained, *ii)* acyclic, and *iii)* a tree. We do
this by executing the following commands:

```bash
python scripts/tlp.py fast_machina examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -l LL -c none -o examples/CP28_unconstrained
python scripts/tlp.py fast_machina examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -l LL -c dag -o examples/CP28_dag
python scripts/tlp.py fast_machina examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -l LL -c tree -o examples/CP28_tree
```

Here, the `-l` flag specifies the anatomical location of the root of the tree,
which in this case is the left lung (LL). If this flag is not provided, the
root is automatically inferred. The `-c` flag specifies the constraint on the
induced migration graph, which can be `none`, `dag`, or `tree`. This procedure 
results in the following output files:
- `examples/CP28_unconstrained_migration_graph.csv`
- `examples/CP28_dag_migration_graph.csv`
- `examples/CP28_tree_migration_graph.csv`
which describe the induced migration multigraphs for the unconstrained,
acyclic, and tree constraints, respectively.

#### Example 2: The Convex Recoloring Problem

To solve the convex recoloring problem, we apply the tree labeling polytope
to the tumor clone `CP28` from the same single-cell lineage tracing dataset.
We use the same tree topology and leaf labeling as in the previous example.

```bash
python scripts/tlp.py convex_recoloring examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -o examples/CP28_convex_recoloring
```

This command results in the output file `examples/CP28_convex_recoloring_vertex_labeling.csv`
which describes the convex recoloring of the tree with the minimum number of leaf recolorings.
