params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"

params.ncells   = [100, 200, 500, 1000]                     // number of sampled cells
params.mrate    = [5e-4]                                    // migration rate
params.settings = ['polyclonal_dag']                        // structure
params.seeds    = 1..10                                     // random seed
params.error    = [5]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'
    stageInMode 'copy'

    publishDir "${params.output_dir}/ground_truth/n${cells}_m${mrate}_e${error}_s${seed}_${setting}", mode: 'copy', overwrite: true

    input:
        tuple val(cells), val(mrate), val(setting), val(error), val(seed)

    output:
        tuple path("sim_labeling.csv"), path("sim_leaf_labeling.csv"), path("sim_leaf_labeling.tsv"), 
        path("sim_migration_graph.csv"), path("sim_tree_edgelist.tsv"), path("sim_perturbed_tree_edgelist.tsv"),
        val(cells), val(mrate), val(setting), val(seed), val("n${cells}_m${mrate}_e${error}_s${seed}_${setting}")

    """
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -n ${cells} --migration-rate ${mrate} -r ${seed} -s ${setting} --generations 44 -e ${error}
    tail -n +2 sim_leaf_labeling.csv | sed 's/,/\t/' > sim_leaf_labeling.tsv
    """
}

process convex_recoloring {
    cpus 16
    memory '16 GB'
    time '4h'
    stageInMode 'copy'

    publishDir "${params.output_dir}/convex_recoloring/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("inferred_results.json"), path("timing.txt"), val(id)

    """
    module load gurobi
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tlp.py convex_recoloring ${edgelist} ${leaf_labeling} -o inferred &> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.ncells)
                               .combine(channel.fromList(params.mrate))
                               .combine(channel.fromList(params.settings))
                               .combine(channel.fromList(params.error))
                               .combine(channel.fromList(params.seeds))

    simulation = parameter_channel | create_sim 
    convex_recoloring_res = simulation | map {[it[1], it[5], it[8], it[10]]} | convex_recoloring 
}
