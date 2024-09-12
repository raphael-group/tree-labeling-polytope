params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/hypothesis_testing/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.machina    = "/n/fs/ragr-data/bin/pmh"

params.ncells   = [1000]                                                                                    // number of sampled cells
params.mrate    = [5e-4]                                                                                   // migration rate
params.settings = ['monoclonal_tree', 'polyclonal_tree']                                                   // structure
params.seeds    = 1..200                                                                                   // number of random seeds
params.error    = [0, 5]                                                                                   // number of errors

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'
    stageInMode 'copy'

    input:
        tuple val(cells), val(mrate), val(setting), val(seed), val(error)

    output:
        tuple path("sim_labeling.csv"), path("sim_leaf_labeling.csv"), path("sim_leaf_labeling.tsv"), 
        path("sim_migration_graph.csv"), path("sim_tree_edgelist.tsv"), path("sim_perturbed_tree_edgelist.tsv"),
        val(cells), val(mrate), val(setting), val(seed), 
        val("n${cells}_m${mrate}_s${seed}_e${error}_${setting}")


    """
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -n ${cells} --migration-rate ${mrate} -r ${seed} -s ${setting} --generations 44 -e ${error}
    tail -n +2 sim_leaf_labeling.csv | sed 's/,/\t/' > sim_leaf_labeling.tsv
    """
}

process parsimonious_relabeling {
    cpus 16
    memory '16 GB'
    time '59m'
    stageInMode 'copy'
    errorStrategy 'ignore'

    publishDir "${params.output_dir}/parsimonious_relabeling/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("inferred_results.json"), 
              path("timing.txt"), val(id)

    """
    module load gurobi
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tlp.py parsimonious_relabeling ${edgelist} ${leaf_labeling} -o inferred &> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.ncells)
                               .combine(channel.fromList(params.mrate))
                               .combine(channel.fromList(params.settings))
                               .combine(channel.fromList(params.seeds))
                               .combine(channel.fromList(params.error))

    simulation = parameter_channel | create_sim 
    parsimonious_relabeling_results = simulation | map {[it[1], it[5], it[8], it[10]]} | parsimonious_relabeling 
  }
