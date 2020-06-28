#!/bin/bash

rsem-generate-data-matrix *[0-9][HF][1-9]_RNA.genes.results \
*i[HF][1-9]_RNA.genes.results > RNA.counts.matrix

rsem-generate-data-matrix *w[HF][1-9]_RNA.genes.results > wRNA.counts.matrix
