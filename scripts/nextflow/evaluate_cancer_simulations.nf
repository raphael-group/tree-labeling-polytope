nextflow.enable.dsl=2

params.python           = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.root_dir         = "/n/fs/ragr-research/projects/tree-labeling-polytope/"
params.output_dir       = "${params.root_dir}/nextflow_results/"
params.ground_truth_dir = "${params.root_dir}/nextflow_results/ground_truth/"

params.algorithms = [
    ['fast_machina', 'nextflow_results/fast_machina/', '_labeling.csv'],
    //['machina', 'nextflow_results/machina/', '_labeling.csv'],
    //['exact_tnet', 'nextflow_results/exact_tnet/', '_labeling.csv'],
    //['tnet', 'nextflow_results/tnet/', '_labeling.csv'],
    ['parsimony', 'nextflow_results/parsimony/', '_labeling.csv'],
]

process EvaluateLabelings {
    stageInMode 'copy'
    errorStrategy 'ignore'

    publishDir "${params.output_dir}/evaluation/${algo}_${name}/", mode: 'copy', overwrite: true

    input:
        tuple val(algo), val(name), path(perturbed_tree), path(tree), 
        path(true_labeling, stageAs: 'true_labeling.csv'), 
        path(inferred_labeling, stageAs: 'inferred_labeling.csv'), path(timing)

    output:
        path("results.json")
        //, path("inferred_colored_tree.svg"), path("inferred_color_graph.svg")

    """
    ${params.python} ${params.root_dir}/scripts/processing/score_result.py ${perturbed_tree} ${tree} ${true_labeling} ${inferred_labeling} ${timing} -o results.json
    """
    //${params.python} ${params.root_dir}/scripts/plots/draw_colored_tree.py ${perturbed_tree} ${inferred_labeling} -o inferred --svg
}

workflow {
    instances = Channel
        .fromPath(params.ground_truth_dir + '/*/sim_tree_edgelist.tsv')
        .map { file -> 
            words = file.toString().split('/')
            words[words.size() - 2]
        }
        .unique()

    algorithms_ch = Channel.fromList(params.algorithms)
    eval_channel = algorithms_ch.combine(instances).map { 
        algo = it[0]
        tree              = "${params.ground_truth_dir}/${it[3]}/sim_tree_edgelist.tsv"
        perturbed_tree    = "${params.ground_truth_dir}/${it[3]}/sim_perturbed_tree_edgelist.tsv"
        true_labeling     = "${params.ground_truth_dir}/${it[3]}/sim_labeling.csv"
        inferred_labeling = "${params.root_dir}/${it[1]}/${it[3]}/inferred_vertex_labeling.csv"
        timing            = "${params.root_dir}/${it[1]}/${it[3]}/timing.txt"
        [algo, it[3], perturbed_tree, tree, true_labeling, inferred_labeling, timing]
    }

    eval_channel | filter { file(it[4]).exists() } | EvaluateLabelings
}
