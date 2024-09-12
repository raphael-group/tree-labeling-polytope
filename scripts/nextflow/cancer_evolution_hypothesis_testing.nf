params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/hypothesis_testing/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.machina    = "/n/fs/ragr-data/bin/pmh"

params.ncells   = [500]                                                                                    // number of sampled cells
params.mrate    = [5e-4]                                                                                   // migration rate
params.settings = ['polyclonal_tree', 'polyclonal_dag', 'none']                                            // structure
params.seeds    = 1..2//00                                                                                   // number of random seeds
params.error    = [0, 5]                                                                                   // number of errors

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'
    stageInMode 'copy'

    publishDir "${params.output_dir}/ground_truth/n${cells}_m${mrate}_s${seed}_e${error}_${input_setting}_vs_${output_setting}", mode: 'copy', overwrite: true

    input:
        tuple val(cells), val(mrate), val(input_setting), val(output_setting), val(seed), val(error)

    output:
        tuple path("sim_labeling.csv"), path("sim_leaf_labeling.csv"), path("sim_leaf_labeling.tsv"), 
        path("sim_migration_graph.csv"), path("sim_tree_edgelist.tsv"), path("sim_perturbed_tree_edgelist.tsv"),
        val(cells), val(mrate), val(input_setting), val(output_setting), val(seed), 
        val("n${cells}_m${mrate}_s${seed}_e${error}_${input_setting}_vs_${output_setting}")

    //path("sim_colored_tree.svg"), path("sim_color_graph.svg")

    """
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -n ${cells} --migration-rate ${mrate} -r ${seed} -s ${input_setting} --generations 44 -e ${error}
    tail -n +2 sim_leaf_labeling.csv | sed 's/,/\t/' > sim_leaf_labeling.tsv
    """
    // ${params.python} ${params.scripts_dir}/plots/draw_colored_tree.py sim_tree_edgelist.tsv sim_labeling.csv -o sim --svg
}

process fast_machina {
    cpus 16
    memory '16 GB'
    time '59m'
    stageInMode 'copy'
    errorStrategy 'ignore'

    publishDir "${params.output_dir}/fast_machina/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(input_setting), val(output_setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("inferred_results.json"), 
              path("transition_probs.csv"), path("timing.txt"), val(id)

    """
    module load gurobi
    ${params.python} ${params.scripts_dir}/processing/branch_lengths_to_probs.py ${edgelist} ${leaf_labeling} --format edgelist > transition_probs.csv
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -c ${output_setting} -o inferred 2>> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.ncells)
                               .combine(channel.fromList(params.mrate))
                               .combine(channel.fromList(params.settings))
                               .combine(channel.fromList(params.settings))
                               .combine(channel.fromList(params.seeds))
                               .combine(channel.fromList(params.error))

    simulation = parameter_channel | create_sim 
    fast_machina_results = simulation | map {[it[1], it[5], it[8], it[9], it[11]]} | fast_machina 
  }
