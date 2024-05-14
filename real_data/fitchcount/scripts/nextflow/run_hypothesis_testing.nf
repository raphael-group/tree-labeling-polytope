params.proj_dir           = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.scripts_dir        = "${params.proj_dir}/scripts/"
params.real_data_dir      = "${params.proj_dir}/real_data/fitchcount/"
params.real_data_scripts  = "${params.proj_dir}/real_data/fitchcount/scripts/"

params.python  = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"

params.tumor_coloring = "/n/fs/ragr-research/projects/tree-labeling-polytope/real_data/fitchcount/input/tumor_coloring.csv"

process fitch {
    stageInMode 'copy'
    publishDir "results/fitch/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(id)

    output:
        path("results.json")

    """
    ${params.python} ${params.scripts_dir}/fitch.py ${edgelist} ${leaf_labeling} -l LL -s 10000 -o results.json
    """
}

process parsimonious_relabeling {
    stageInMode 'copy'
    publishDir "results/parsimonious_relabeling/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("inferred_results.json"), 
        path("log.txt"), val(id)

    """
    ${params.python} ${params.scripts_dir}/tlp.py parsimonious_relabeling ${edgelist} ${leaf_labeling} -o inferred &> log.txt
    """
}

process fast_machina {
    stageInMode 'copy'

    publishDir "results/fast_machina/${id}/${weighted}/${setting}/", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), path(newick), val(setting), val(weighted), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("inferred_results.json"), 
        path("log.txt"), path("*.svg"), path("*.dot"), val(id)

    script:
    if (weighted == 'weighted') {
        """
        ${params.python} ${params.scripts_dir}/processing/branch_lengths_to_probs.py -f edgelist ${edgelist} ${leaf_labeling} > probabilities.csv
        ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -l LL -c ${setting} -o inferred -w probabilities.csv > log.txt
        ${params.python} ${params.scripts_dir}/plots/draw_colored_tree.py ${edgelist} inferred_vertex_labeling.csv -p ${params.tumor_coloring} -f edgelist --svg -m
        """
    } else if (weighted == 'unweighted') {
        """
        ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -l LL -c ${setting} -o inferred > log.txt
        ${params.python} ${params.scripts_dir}/plots/draw_colored_tree.py ${edgelist} inferred_vertex_labeling.csv -p ${params.tumor_coloring} -f edgelist --svg -m
        """
    } else {
        error "Invalid weighted parameter."
    }
}

workflow {
    inferred_trees = Channel.fromPath("results/laml/*", type: "dir") | map {
        ["${it}/laml_leaf_labeling.csv", "${it}/laml_tree_edgelist.tsv", "${it}/laml_tree.newick", it.baseName]
    } 

    settings = ['none', 'polyclonal_tree', 'polyclonal_dag']

    inferred_trees | combine(settings) | combine(['weighted', 'unweighted']) | map {
        leaf_labeling, laml_edgelist, laml_newick, id, setting, weighted ->
        [leaf_labeling, laml_edgelist, laml_newick, setting, weighted, id]
    } | fast_machina

    inferred_trees | map {
        leaf_labeling, laml_edgelist, laml_newick, id ->
        [leaf_labeling, laml_edgelist, id]
    } | fitch

    inferred_trees | map {
        leaf_labeling, laml_edgelist, laml_newick, id ->
        [leaf_labeling, laml_edgelist, id]
    } | parsimonious_relabeling
}
