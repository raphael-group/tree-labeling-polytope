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
        tuple file("sim_labeling.csv"), file("sim_leaf_labeling.csv"), file("sim_migration_graph.csv"), 
        file("sim_tree_edgelist.tsv"), file("sim_tree.newick"), 
        val(cells), val(labels), val(setting), val(seed), val("n${cells}_m${labels}_s${seed}_${setting}")

    """
    ${params.ngesh} -l ${cells} -L ${cells} -x enum -r ${seed} > sim_tree.newick
    ${params.python} ${params.scripts_dir}/newick_to_edgelist.py sim_tree.newick > sim_tree_edgelist.tsv
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -m ${labels} -r ${seed}\
                     sim_tree_edgelist.tsv root -s ${setting}
    """
}

workflow {
    create_sim([10, 5, 'polyclonal_tree', 1]) 
}
