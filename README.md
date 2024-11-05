# The Tree Labeling Polytope

The tree labeling polytope is a technical tool for fixed topology 
ancestral reconstruction problems. Here, we provide the implementation
of the application of the tree labeling polytope 
to the following three problems:
- the parsimonious migration history problem
- the softwired small parsimony problem 
- the convex recoloring problem

## Usage

### Example 1: The Parsimonious Migration History Problem

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
python src/tlp.py fast_machina examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -l LL -c none -o examples/CP28_unconstrained
python src/tlp.py fast_machina examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -l LL -c dag -o examples/CP28_dag
python src/tlp.py fast_machina examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -l LL -c tree -o examples/CP28_tree
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

### Example 2: The Convex Recoloring Problem

```bash
python src/tlp.py convex_recoloring examples/CP28_tree_edgelist.tsv examples/CP28_leaf_labeling.csv -o examples/CP28_convex_recoloring
```
