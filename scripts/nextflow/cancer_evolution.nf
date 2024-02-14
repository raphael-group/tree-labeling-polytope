params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.ngesh      = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/ngesh"

params.ncells   = [10, 15, 20, 25, 50]
params.mlabels  = [3, 5, 10]
params.settings = ['monoclonal_tree', 'monoclonal_dag', 'polyclonal_tree', 'polyclonal_dag']
params.seeds    = [1, 2, 3]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'

    input:
        tuple val(cells), val(labels), val(setting), val(seed)

    output:
        tuple path("sim_labeling.csv"), path("sim_leaf_labeling.csv"), path("sim_migration_graph.csv"), 
        path("sim_tree_edgelist.tsv"), path("sim_tree.newick"), 
        val(cells), val(labels), val(setting), val(seed), val("n${cells}_m${labels}_s${seed}_${setting}")

    """
    ${params.ngesh} -l ${cells} -L ${cells} -x enum -r ${seed} > sim_tree.newick
    ${params.python} ${params.scripts_dir}/newick_to_edgelist.py sim_tree.newick > sim_tree_edgelist.tsv
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -m ${labels} -r ${seed}\
                     sim_tree_edgelist.tsv root -s ${setting}
    """
}

process fast_machina {
    cpus 16
    memory '8 GB'
    time '48h'

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), val(id)

    """
    ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -c ${setting} -o inferred
    """
}

workflow {
    create_sim([10, 5, 'polyclonal_tree', 1]) | map {[it[1], it[3], it[7], it[9]]} | fast_machina
}
