nextflow.enable.dsl=2

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.root_dir = '/n/fs/ragr-research/projects/tree-labeling-polytope/'
params.ground_truth_dir = 'nextflow_results/ground_truth/'
params.algorithms = [
    ['fast_machina', 'nextflow_results/fast_machina/', '_labeling.csv'],
    ['machina', 'nextflow_results/machina/', '_labeling.csv'],
    ['exact_tnet', 'nextflow_results/exact_tnet/', '_labeling.csv'],
]

def trimSuffix(original, suffix) {
	if(original.endsWith(suffix)) {
	    return original.substring(0, original.length() - suffix.length())
	}
	return original
}

process EvaluateLabelings {
    errorStrategy 'ignore'

    input:
        tuple val(algo), val(name), path(tree), path(true_labeling, stageAs: 'true_labeling.csv'), path(inferred_labeling, stageAs: 'inferred_labeling.csv'), path(timing)

    output:
        tuple path("results.json"), val(algo), val(name)

    """
    ${params.python} ${params.root_dir}/scripts/processing/score_result.py ${tree} ${true_labeling} ${inferred_labeling} ${timing} -o results.json
    """
}

workflow {
    ground_truth_trees_ch = Channel
        .fromPath(params.ground_truth_dir + '*_tree.newick')
        .map { file -> trimSuffix(file.baseName, '_tree') }
        .unique()

    algorithms_ch = Channel.fromList(params.algorithms)
    eval_channel = algorithms_ch.combine(ground_truth_trees_ch).map { 
        algo = it[0]
        tree              = "${params.root_dir}/${params.ground_truth_dir}/${it[3]}_tree_edgelist.csv"
        true_labeling     = "${params.root_dir}/${params.ground_truth_dir}/${it[3]}_labeling.csv"
        inferred_labeling = "${params.root_dir}/${it[1]}/${it[3]}${it[2]}"
        timing            = "${params.root_dir}/${it[1]}/${it[3]}_timing.txt"
        [algo, it[3], tree, true_labeling, inferred_labeling, timing]
    }

    eval_channel | filter { file(it[4]).exists() } | EvaluateLabelings | map { result, algo, name ->
        result.moveTo("nextflow_results/evaluation/${algo}_${name}.json")
    }
}
