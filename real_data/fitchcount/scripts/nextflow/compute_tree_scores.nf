params.proj_dir           = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir         = "/n/fs/ragr-research/projects/tree-labeling-polytope/real_data/fitchcount/character_matrices/"
params.scripts_dir        = "${params.proj_dir}/scripts/"
params.real_data_dir      = "${params.proj_dir}/real_data/fitchcount/"
params.real_data_scripts  = "${params.proj_dir}/real_data/fitchcount/scripts/"

params.python  = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.startle = "/n/fs/ragr-research/projects/startle/build/src/startle"
params.laml    = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/run_laml"

params.tumor_coloring = "/n/fs/ragr-research/projects/tree-labeling-polytope/real_data/fitchcount/tumor_coloring.csv"

process startle {
    cpus 1
    memory '8 GB'
    time '59m'
    stageInMode 'copy'
    errorStrategy 'ignore'

    publishDir "${method}/${id}/", mode: 'copy', overwrite: true

    input:
        tuple val(method), val(id), path(character_matrix), path(mutation_priors), path(tree)

    output:
        path("parsimony_score.txt")

    """
    startle small ${character_matrix} ${mutation_priors} ${tree} -o test > parsimony_score.txt
    cat parsimony_score.txt | tail -n 1 | rev | cut -d' ' -f1 | rev > tmp_parsimony_score.txt
    mv tmp_parsimony_score.txt parsimony_score.txt
    """
}

process laml {
    cpus 1
    memory '8 GB'
    time '59m'
    stageInMode 'copy'
    errorStrategy 'ignore'

    publishDir "${method}/${id}/", mode: 'copy', overwrite: true

    input:
        tuple val(method), val(id), path(character_matrix), path(mutation_priors), path(tree)

    output:
        path("laml_params.txt")

    """
    export MOSEKLM_LICENSE_FILE=/n/fs/grad/hs2435
    ${params.laml} -c ${character_matrix} -p ${mutation_priors} -t ${tree} -m -1 -o laml --nInitials 1 --delimiter comma
    """
}

workflow {
    instances = [
        'm5k_lg100', 'm5k_lg10', 'm5k_lg11', 'm5k_lg12', 'm5k_lg13', 'm5k_lg14', 'm5k_lg15', 'm5k_lg17',
        'm5k_lg19', 'm5k_lg1', 'm5k_lg20', 'm5k_lg21', 'm5k_lg22', 'm5k_lg23', 'm5k_lg24',
        'm5k_lg26', 'm5k_lg27', 'm5k_lg28', 'm5k_lg29', 'm5k_lg2', 'm5k_lg30', 'm5k_lg31', 'm5k_lg32', 
        'm5k_lg34', 'm5k_lg35', 'm5k_lg36', 'm5k_lg37', 'm5k_lg3', 'm5k_lg40', 
        'm5k_lg42', 'm5k_lg43', 'm5k_lg44', 'm5k_lg45', 'm5k_lg46', 'm5k_lg47', 'm5k_lg48', 'm5k_lg49', 'm5k_lg4',
        'm5k_lg50', 'm5k_lg51', 'm5k_lg52','m5k_lg54', 'm5k_lg55', 'm5k_lg56', 'm5k_lg57', 'm5k_lg58',
        'm5k_lg59', 'm5k_lg5', 'm5k_lg61', 'm5k_lg62', 'm5k_lg63', 'm5k_lg64', 'm5k_lg66',
        'm5k_lg67', 'm5k_lg68', 'm5k_lg6', 'm5k_lg70', 'm5k_lg71', 'm5k_lg72', 'm5k_lg73', 'm5k_lg74',
        'm5k_lg76', 'm5k_lg77', 'm5k_lg78', 'm5k_lg79', 'm5k_lg7', 'm5k_lg80', 'm5k_lg82',
        'm5k_lg83', 'm5k_lg84', 'm5k_lg85', 'm5k_lg86', 'm5k_lg89', 'm5k_lg8', 'm5k_lg90',
        'm5k_lg91', 'm5k_lg92', 'm5k_lg94', 'm5k_lg95', 'm5k_lg96', 'm5k_lg97', 'm5k_lg98', 'm5k_lg99',
        'm5k_lg9'
    ] 

    raw_data = Channel.from(instances) | map {id -> 
        [id,
         "${params.real_data_dir}/preprocessing/${id}/startle_character_matrix.csv",
         "${params.real_data_dir}/preprocessing/${id}/startle_mutation_priors.csv"
        ]
    } 

    laml_trees = raw_data | map {
        id, character_matrix, mutation_priors -> 
        ["startle", id, character_matrix, mutation_priors, "${params.real_data_dir}/startle/${id}/laml_tree.newick"]
    }

    cassiopeia_trees = raw_data | map {
        id, character_matrix, mutation_priors -> 
        ["cassiopeia", id, character_matrix, mutation_priors, 
        "${params.real_data_dir}/preprocessing/${id}/seed_tree.newick"]
    }

    laml_trees.concat(cassiopeia_trees) | laml
}
