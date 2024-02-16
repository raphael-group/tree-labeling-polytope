params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"
params.ngesh      = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/ngesh"
params.machina    = "/n/fs/ragr-data/bin/pmh"

params.ncells   = [10, 15, 20, 25, 50, 100]
params.mlabels  = [3, 5, 10, 20]
params.settings = ['polyclonal_tree', 'polyclonal_dag']
params.seeds    = [1, 2, 3]

params.methods  = ['exact_tnet']

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    stageInMode 'copy'

    input:
        tuple val(cells), val(labels), val(setting), val(seed)

    output:
        tuple path("sim_labeling.csv"), path("sim_leaf_labeling.csv"), path("sim_leaf_labeling.tsv"), 
        path("sim_migration_graph.csv"), path("sim_tree_edgelist.tsv"), path("sim_tree.newick"), 
        val(cells), val(labels), val(setting), val(seed), val("n${cells}_m${labels}_s${seed}_${setting}")

    """
    ${params.ngesh} -l ${cells} -L ${cells} -x enum -r ${seed} > sim_tree.newick
    ${params.python} ${params.scripts_dir}/processing/newick_to_edgelist.py sim_tree.newick > sim_tree_edgelist.tsv
    ${params.python} ${params.scripts_dir}/simulations/cancer_evolution.py -o sim -m ${labels} -r ${seed}\
                     sim_tree_edgelist.tsv root -s ${setting}
    tail -n +2 sim_leaf_labeling.csv | sed 's/,/\t/' > sim_leaf_labeling.tsv
    """
}

process fast_machina {
    cpus 8
    memory '8 GB'
    time '4h'
    stageInMode 'copy'

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("timing.txt"), val(id)

    """
    /usr/bin/time -v ${params.python} ${params.scripts_dir}/tlp.py fast_machina ${edgelist} ${leaf_labeling} -c ${setting} -o inferred 2>> timing.txt
    """
}

process exact_tnet {
    cpus 8
    memory '8 GB'
    time '4h'
    stageInMode 'copy'

    input:
        tuple path(leaf_labeling), path(edgelist), val(setting), val(id)

    output:
        tuple path("inferred_vertex_labeling.csv"), path("inferred_migration_graph.csv"), path("timing.txt"), val(id)

    """
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
    printf "dummy_root\troot\ndummy_root\tdummy_leaf\n" >> ${edgelist} 
    printf "dummy_leaf\tdummy_label\n" >> ${leaf_labeling}
    echo "dummy_label" >> root_label.txt
    ${params.python} ${params.scripts_dir}/processing/create_machina_input.py ${labeling}
    """
}

process machina {
    cpus 8
    memory '32 GB'
    time '24h'
    stageInMode 'copy'

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

workflow {
    parameter_channel = channel.fromList(params.ncells)
                               .combine(channel.fromList(params.mlabels))
                               .combine(channel.fromList(params.settings))
                               .combine(channel.fromList(params.seeds))

    simulation = parameter_channel | create_sim 

    // create directories
    file("${params.output_dir}/exact_tnet/").mkdirs()
    file("${params.output_dir}/fast_machina/").mkdirs()
    file("${params.output_dir}/machina/").mkdirs()
    file("${params.output_dir}/ground_truth/").mkdirs()

    // move input files to nextflow_results
    simulation | map {
        labeling, leaf_labeling_csv, leaf_labeling_tsv, migration_graph, edgelist, newick,
        cells, labels, setting, seed, id ->

        output_prefix = "${params.output_dir}/ground_truth/${id}"
        labeling.copyTo("${output_prefix}_labeling.csv")
        leaf_labeling_csv.copyTo("${output_prefix}_leaf_labeling.csv")
        leaf_labeling_tsv.copyTo("${output_prefix}_leaf_labeling.tsv")
        migration_graph.copyTo("${output_prefix}_migration_graph.csv")
        edgelist.copyTo("${output_prefix}_tree_edgelist.csv")
        newick.copyTo("${output_prefix}_tree.newick")
    }

    // run all methods
   
    if (params.methods.contains('fast_machina')) {
        fast_machina_results = simulation | map {[it[1], it[4], it[8], it[10]]} | fast_machina 
        fast_machina_results | map {
            inferred_labeling, inferred_migration_graph, timing, id ->

            output_prefix = "${params.output_dir}/fast_machina/${id}"
            inferred_labeling.toFile().copyTo("${output_prefix}_labeling.csv")
            inferred_migration_graph.toFile().copyTo("${output_prefix}_migration_graph.csv")
            timing.toFile().copyTo("${output_prefix}_timing.txt")
        }
    }

     if (params.methods.contains('exact_tnet')) { 
        exact_tnet_results = simulation | map {[it[1], it[4], it[8], it[10]]} | exact_tnet 
        exact_tnet_results | map {
            inferred_labeling, inferred_migration_graph, timing, id ->

            output_prefix = "${params.output_dir}/exact_tnet/${id}"
            inferred_labeling.toFile().copyTo("${output_prefix}_labeling.csv")
            inferred_migration_graph.toFile().copyTo("${output_prefix}_migration_graph.csv")
            timing.toFile().copyTo("${output_prefix}_timing.txt")
        }
    }

    if (params.methods.contains('machina')) {
        machina_input = simulation | map {[it[0], it[2], it[4], it[8], it[10]]} | create_machina_input | map {
            root_label = it[1].text.trim()
            setting_map = ["polyclonal_tree": 1, "polyclonal_dag": 2]
            setting = setting_map[it[4]]
            [it[2], it[3], it[0], root_label, setting, it[5]]
        }

        machina_results      = machina_input | machina
        machina_results | map {
            inferred_labeling, timing, id ->
            output_prefix = "${params.output_dir}/machina/${id}"
            inferred_labeling.toFile().copyTo("${output_prefix}_labeling.csv")
            timing.toFile().copyTo("${output_prefix}_timing.txt")
        }
    }
}
