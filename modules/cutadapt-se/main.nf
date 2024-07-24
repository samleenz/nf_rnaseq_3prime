process RUN_TRIM_READS_SE {
    cpus 4
    memory "8.G"
    time "15.m"
    container "quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1"
    tag "$sample"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("tr_*")

    script:
    out1 = "tr_${reads[0]}"
    """
    cutadapt \
        -m 20 \
        -O 20 \
        -a "polyA=A{20}" \
        -a "QUALITY=G{20}" \
        -n 2 \
        ${reads[0]} | \
    cutadapt \
        -m 20 \
        -O 3 \
        --nextseq-trim=10 \
        -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" – | \
    cutadapt \
        -m 20 \
        -O 20 \
        -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" \
        --discard-trimmed \
        -o $out1 –


    cutadapt \
        -a AGATCGGAAGAG \
        -A AGATCGGAAGAG \
        -j $task.cpus \
        -q 20 \
        -o $out1 \
        -p $out2 \
        -m 15 \
        ${reads[0]} \
        ${reads[1]} 


    """

}
