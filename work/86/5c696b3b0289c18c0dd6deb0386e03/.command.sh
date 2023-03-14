#!/bin/bash -ue
salmon quant --threads 1 --libType=U -i salmon_index -1 liver_1.fq -2 liver_2.fq -o liver
