params.proj_dir    = "/n/fs/ragr-research/projects/tree-labeling-polytope"
params.output_dir  = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/"
params.scripts_dir = "${params.proj_dir}/scripts/"

params.python     = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"

params.nreticulations = [0, 3, 5, 10, 20, 40]
params.nleaves        = [100, 500, 1000, 2000]
params.alphabet       = [20]
params.seeds          = 0..10

process simulate {
    cpus 1
    memory '4 GB'
    time '59m'

    input:
        tuple val(ntaxa), val(nloci), val(nreticulations), val(alphabet), val(seed)
    output:
        tuple path("sim_network.edgelist"), path("sim_labeling.csv"), val("n${ntaxa}_m${nloci}_r${nreticulations}_s${seed}_a${alphabet}")

    publishDir "${params.output_dir}/softwired/ground_truth/n${ntaxa}_m${nloci}_r${nreticulations}_s${seed}_a${alphabet}", mode: 'copy', overwrite: true

    """
    ${params.python} ${params.scripts_dir}/simulations/phylogenetic_networks.py -m ${nloci} -n ${ntaxa} -r ${nreticulations} -s ${seed} -a ${alphabet} -o sim --uniform
    """
}

process softwired_tlp {
    cpus 1
    memory '8 GB'
    time '59m'

    input:
        tuple path(edgelist), path(labeling), val(formulation), val(id)

    output:
        tuple path("inferred_results.json"), path("log.txt")

    publishDir "${params.output_dir}/softwired/${formulation}/${id}", mode: 'copy', overwrite: true
    
    """
    ${params.python} ${params.scripts_dir}/softwired_tlp.py ${edgelist} ${labeling} inferred --mode ${formulation} &> log.txt
    """
}

workflow {
    simulations = Channel.fromList(params.nleaves)
                         .combine(Channel.fromList([1]))
                         .combine(Channel.fromList(params.nreticulations))
                         .combine(Channel.fromList(params.alphabet))
                         .combine(Channel.fromList(params.seeds))
    simulations | simulate | flatMap {
        [[it[0], it[1], "scornavacca", it[2]], 
         [it[0], it[1], "tlp", it[2]]]
    } | softwired_tlp
}
