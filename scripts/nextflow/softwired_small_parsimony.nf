params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"

params.nreticulations = [5, 10, 15, 20, 25, 30, 35, 40, 80, 160, 320]
params.seeds          = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

process simulate {
    cpus 1
    memory '4 GB'
    time '59m'

    input:
        tuple val(ntaxa), val(nloci), val(nreticulations), val(seed)
    output:
        tuple path("sim_network.edgelist"), path("sim_labeling.csv"), val("n${ntaxa}_m${nloci}_r${nreticulations}_s${seed}")

    publishDir "${params.output_dir}/softwired/ground_truth/n${ntaxa}_m${nloci}_r${nreticulations}_s${seed}", mode: 'copy', overwrite: true

    """
    ${params.python} ${params.scripts_dir}/simulations/phylogenetic_networks.py -m ${nloci} -n ${ntaxa} -r ${nreticulations} -s ${seed} -o sim --uniform
    """
}

process softwired_tlp {
    cpus 1
    memory '16 GB'
    time '59m'

    input:
        tuple path(edgelist), path(labeling), val(id)

    output:
        tuple path("inferred_results.json"), path("log.txt")

    publishDir "${params.output_dir}/softwired/inferred/${id}", mode: 'copy', overwrite: true
    
    """
    ${params.python} ${params.scripts_dir}/softwired_tlp.py ${edgelist} ${labeling} inferred &> log.txt
    """
}

workflow {
    simulations = Channel.fromList([2000])
                         .combine(Channel.fromList([1]))
                         .combine(Channel.fromList(params.nreticulations))
                         .combine(Channel.fromList(params.seeds))
    simulations | simulate | softwired_tlp
}
