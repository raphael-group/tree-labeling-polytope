params.proj_dir           = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir         = "/n/fs/ragr-research/projects/tree-labeling-polytope/real_data/fitchcount/character_matrices/"
params.scripts_dir        = "${params.proj_dir}/scripts/"
params.real_data_dir      = "${params.proj_dir}/real_data/fitchcount/"
params.real_data_scripts  = "${params.proj_dir}/real_data/fitchcount/scripts/"

params.python  = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.scmail  = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/run_scmail"
params.laml    = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/run_laml"

params.tumor_coloring = "/n/fs/ragr-research/projects/tree-labeling-polytope/real_data/fitchcount/input/tumor_coloring.csv"
params.instances      = "/n/fs/ragr-research/projects/tree-labeling-polytope/real_data/fitchcount/input/instances.txt"

process preprocess_fitchcount {
    cpus 1

    stageInMode 'copy'
    publishDir "results/preprocessing/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(raw_character_matrix), path(seed_tree), val(id)
    output:
        tuple path("laml_character_matrix.csv"), path("laml_mutation_priors.csv"), path("laml_leaf_labeling.csv"), 
        path("seed_tree.newick"), val(id)

    """
    ${params.python} ${params.real_data_scripts}/preprocess_character_matrix.py ${raw_character_matrix} laml
    ${params.python} ${params.real_data_scripts}/preprocess_newick.py ${seed_tree} > tree_edgelist.csv
    ${params.python} ${params.scripts_dir}/processing/edgelist_to_newick.py tree_edgelist.csv > seed_tree.newick
    """
}

process laml {
    cpus 32
    memory '8 GB'
    time '7d'
    stageInMode 'copy'

    publishDir "results/laml/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(character_matrix), path(mutation_priors), path(leaf_labeling), path(seed_tree), val(id)

    output:
        tuple path(leaf_labeling), path("laml_tree_edgelist.tsv"), path("laml_tree.newick"), path("laml_params.txt"), path("laml_annotations.txt"),
              path("laml.log"), path("seed_params.txt"), val(id)

    """
    export MOSEKLM_LICENSE_FILE=/n/fs/grad/hs2435
    ${params.laml} -t ${seed_tree} -c ${character_matrix} -p ${mutation_priors} -o seed --delimiter comma --nInitials 1 -m -1 
    ${params.laml} -t ${seed_tree} -c ${character_matrix} -p ${mutation_priors} -o laml --delimiter comma --nInitials 1 -m -1 --topology_search --parallel --maxIters 50 --randomreps 1
    mv laml_trees.nwk laml_tree.newick
    ${params.python} ${params.scripts_dir}/processing/newick_to_edgelist.py -b laml_tree.newick > laml_tree_edgelist.tsv
    ${params.python} ${params.scripts_dir}/processing/edgelist_to_newick.py -b laml_tree_edgelist.tsv > laml_tree.newick
    """
}

workflow {
    instances = Channel.from(file(params.instances).readLines()) | map {"m5k_lg${it}"} 
    inferred_trees = instances | filter {!file("results/laml/${it}/laml_tree.newick").exists()} | map {id -> 
        ["${params.real_data_dir}/input/character_matrices/${id}_character_matrix.alleleThresh.txt", 
         "${params.real_data_dir}/input/trees/${id}_tree_hybrid_priors.alleleThresh.processed.txt", 
         id]
    } | preprocess_fitchcount | laml
}
