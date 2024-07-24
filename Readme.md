# Quantseq 3' pipeline

Sam Lee
email: samuel.lee@onjcri.org.au


## To run

```bash
# start a tmux / screen session

module load nextflow

nextflow run main.nf --params-file params.yml -resume -with-dag flowchart.mmd
```

## Analysis steps

Using (https://github.com/samleenz/nf_modules) for analysis code as well as local versions (for single end reads) under `modules/`. 

- Quality control
  - fastqc
  - multiQC
- Read Processing
  - Trimming with cutadapt
- Alignment
  - STAR two-pass alignment
- Quantification  
  - featureCounts 

## Versions


