#!/bin/bash -ue
mkdir fastqc_liver_logs
fastqc -o fastqc_liver_logs -f fastq -q liver_1.fq liver_2.fq
