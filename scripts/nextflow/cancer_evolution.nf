params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.machina    = "/n/fs/ragr-data/bin/pmh"

params.ncells   = [250, 500, 750, 1000]                                     // number of sampled cells
params.mrate    = [1e-3]                                                    // migration rate
params.settings = ['polyclonal_tree', 'polyclonal_dag']                     // structure
params.seeds    = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]       // random parameter

// params.methods  = ['tnet', 'fast_machina', 'parsimony']
params.methods = ['fast_machina']

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'
    stageInMode 'copy'

    publishDir "${params.output_dir}/ground_truth/n${cells}_m${mrate}_s${seed}_${setting}", mode: 'copy', overwrite: true

    input:
        tuple val(cells), val(mrate), val(setting), val(seed)

    output:
        tuple path("sim_labeling.csv"), path("sim_leaf_labeling.csv"), path("sim_leaf_labeling.tsv"), 
        path("sim_migration_graph.csv"), path("sim_tree_edgelist.tsv"), path("sim_perturbed_tree_edgelist.tsv"),
        val(cells), val(mrate), val(setting), val(seed), val("n${cells}_m${mrate}_s${seed}_${setting}"), 
        path("sim_colored_tree.svg"), path("sim_color_graph.svg")

    """
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -n ${cells} --migration-rate ${mrate} -r ${seed} -s ${setting} --generations 44 -e 10
    tail -n +2 sim_leaf_labeling.csv | sed 's/,/\t/' > sim_leaf_labeling.tsv
    ${params.python} ${params.scripts_dir}/plots/draw_colored_tree.py sim_tree_edgelist.tsv sim_labeling.csv -o sim --svg
    """
}

process fast_machina {
    cpus 16
    memory '16 GB'
    time '24h'
    stageInMode 'copy'
    errorStrategy 'ignore'

    publishDir "${params.output_dir}/fast_machina/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("timing.txt"), val(id)

    """
    module load gurobi
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -c ${setting} -o inferred 2>> timing.txt
    """
}

process parsimony {
    cpus 8
    memory '16 GB'
    time '4h'
    stageInMode 'copy'

    publishDir "${params.output_dir}/parsimony/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("timing.txt"), val(id)

    """
    module load gurobi
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -c none -o inferred 2>> timing.txt
    """
}

process exact_tnet {
    cpus 8
    memory '16 GB'
    time '4h'
    stageInMode 'copy'

    publishDir "${params.output_dir}/exact_tnet/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("timing.txt"), val(id)

    """
    module load gurobi
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tlp.py exact_tnet ${edgelist} ${leaf_labeling} -o inferred 2>> timing.txt
    """
}

process create_machina_input {
    cpus 1
    memory '4 GB'
    time '59m'
    stageInMode 'copy'

    input: 
        tuple path(labeling), path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("machina_coloring.txt"), path("root_label.txt"), path(leaf_labeling), path(edgelist), val(setting), val(id)

    """
    printf "dummy_root\ts0\ndummy_root\tdummy_leaf\n" >> ${edgelist} 
    printf "dummy_leaf\tdummy_label\n" >> ${leaf_labeling}
    echo "dummy_label" >> root_label.txt
    ${params.python} ${params.scripts_dir}/processing/create_machina_input.py ${labeling}
    """
}

process machina {
    cpus 8
    memory '32 GB'
    time '24h'
    errorStrategy 'ignore'
    stageInMode 'copy'

    publishDir "${params.output_dir}/machina/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(leaf_labeling), path(edgelist), path(coloring), val(root_label), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("timing.txt"), val(id)
    
    """
    module load gurobi
    mkdir machina/
    /usr/bin/time -v ${params.machina} ${edgelist} ${leaf_labeling} -p ${root_label} -c ${coloring} -m ${setting} -o machina 2>> timing.txt
    mv machina/*.labeling inferred_vertex_labeling.txt
    echo 'vertex,label' | cat - inferred_vertex_labeling.txt | sed 's/ /,/g' > inferred_vertex_labeling.csv
    """
}

process create_tnet_input {
    cpus 1
    memory '4 GB'
    time '59m'
    stageInMode 'copy'

    input: 
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("tnet_input.newick"), path(leaf_labeling), path(edgelist), val(setting), val(id)

    """
    ${params.python} ${params.scripts_dir}/processing/create_tnet_input.py ${edgelist} ${leaf_labeling} > tnet_input.newick
    """
}

process tnet {
    cpus 8
    memory '16 GB'
    time '59m'
    stageInMode 'copy'
    
    publishDir "${params.output_dir}/tnet/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(tnet_input), path(leaf_labeling), path(edgelist), val(setting), val(id)

    output: 
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("timing.txt"), val(id)

    """
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tnet.py ${tnet_input} tnet_output_network.tsv 2>> timing.txt
    ${params.python} ${params.scripts_dir}/processing/create_weights_from_network.py tnet_output_network.tsv > network_weights.csv
    module load gurobi
    ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -w network_weights.csv -c none -o inferred
    """
}

workflow {
    parameter_channel = channel.fromList(params.ncells)
                               .combine(channel.fromList(params.mrate))
                               .combine(channel.fromList(params.settings))
                               .combine(channel.fromList(params.seeds))

    simulation = parameter_channel | create_sim 

    if (params.methods.contains('fast_machina')) {
        fast_machina_results = simulation | map {[it[1], it[5], it[8], it[10]]} | fast_machina 
    }
   
    if (params.methods.contains('parsimony')) {
        parsimony_results = simulation | map {[it[1], it[5], it[8], it[10]]} | parsimony 
    }

    if (params.methods.contains('exact_tnet')) { 
        exact_tnet_results = simulation | map {[it[1], it[5], it[8], it[10]]} | exact_tnet 
    }

    if (params.methods.contains('tnet')) { 
        tnet_results = simulation | map {[it[1], it[5], it[8], it[10]]} | create_tnet_input | tnet 
    }

    if (params.methods.contains('machina')) {
        machina_input = simulation | map {[it[0], it[2], it[5], it[8], it[10]]} | create_machina_input | map {
            root_label = it[1].text.trim()
            setting_map = ["polyclonal_tree": 1, "polyclonal_dag": 2]
            setting = setting_map[it[4]]
            [it[2], it[3], it[0], root_label, setting, it[5]]
        }

        machina_results = machina_input | machina
    }
}
