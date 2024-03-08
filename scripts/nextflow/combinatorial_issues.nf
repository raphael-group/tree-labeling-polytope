params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.machina    = "/n/fs/ragr-data/bin/pmh"

params.ncells   = [500]                                   // number of sampled cells
params.mrate    = [5e-3]                                          // migration rate
params.settings = ['polyclonal_dag']
params.seeds    = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]                                
params.error    = [0,1]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'
    stageInMode 'copy'

    publishDir "${params.output_dir}/ground_truth/n${cells}_m${mrate}_s${seed}_e${error}_${setting}", mode: 'copy', overwrite: true

    input:
        tuple val(cells), val(mrate), val(error), val(setting), val(seed)

    output:
        tuple path("sim_labeling.csv"), path("sim_leaf_labeling.csv"), path("sim_leaf_labeling.tsv"), 
        path("sim_migration_graph.csv"), path("sim_tree_edgelist.tsv"), path("sim_perturbed_tree_edgelist.tsv"),
        val(cells), val(mrate), val(setting), val(seed), val("n${cells}_m${mrate}_s${seed}_e${error}_${setting}")

    """
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -n ${cells} --migration-rate ${mrate} -r ${seed} -s ${setting} --generations 44 -e ${error} 
    tail -n +2 sim_leaf_labeling.csv | sed 's/,/\t/' > sim_leaf_labeling.tsv
    """
    //${params.python} ${params.scripts_dir}/plots/draw_colored_tree.py sim_tree_edgelist.tsv sim_labeling.csv -o sim --svg
}

process fitch {
    cpus 1
    memory '4 GB'
    time '59m'
    stageInMode 'copy'
    //errorStrategy 'ignore'

    publishDir "${params.output_dir}/combinatorics/fitch/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        path("result.json")

    """
    ${params.python} ${params.scripts_dir}/fitch.py ${edgelist} ${leaf_labeling} -o result.json --samples 10000
    """
}

workflow {
    parameter_channel = channel.fromList(params.ncells)
                               .combine(channel.fromList(params.mrate))
                               .combine(channel.fromList(params.error))
                               .combine(channel.fromList(params.settings))
                               .combine(channel.fromList(params.seeds))

    simulation = parameter_channel | create_sim 
    simulation | map {[it[1], it[5], it[8], it[10]]} | fitch 
}
