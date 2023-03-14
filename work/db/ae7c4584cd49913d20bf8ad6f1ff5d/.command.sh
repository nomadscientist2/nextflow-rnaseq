#!/bin/bash -ue
mkdir fastqc_lung_logs
fastqc -o fastqc_lung_logs -f fastq -q lung_1.fq lung_2.fq
