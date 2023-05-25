params.output = "output.txt"
params.ref = "$baseDir/reference.fasta"
params.sequences = "$baseDir/sequences.fasta"

process split_sequences {
    input:
    path x

    output:
    path 'seq_*'
    
    script:
    """
    split -l 2 $x seq_
    """

}

process clean_nucleotides {
    input:
    path x

    output:
    stdout

    script:
    """
    python $baseDir/clean_nucleotides.py $x
    """
}

process identify_region {
    input:
    stdin

    output:
    stdout

    script:
    """
    python $baseDir/identify_hiv_region_and_direction.py 
    """
}

process check_prrt_reversal {
    input:
    stdin

    output:
    stdout

    script:
    """
    python $baseDir/check_prrt_reversal.py
    """
}

process write_output {
    publishDir "$baseDir"

    input:
    stdin

    output:
    path params.output

    script:
    """
    cat >> $params.output
    """
}

workflow {
    def seqs = Channel.fromPath(params.sequences)
    split_sequences(seqs) | flatten | clean_nucleotides | identify_region | check_prrt_reversal | collect | write_output
}
