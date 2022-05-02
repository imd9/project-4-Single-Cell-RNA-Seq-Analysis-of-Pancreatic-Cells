# Project 4: Single Cell RNA-Seq Analysis of Pancreatic Cells

The pancreas is a complex organ comprised of a diverse set of cell types. Proper function of the pancreas is required to maintain healthy metabolism, and pancreatic dysfunction leads to serious illnesses, including type 1 diabetes. Baron et al performed single cell RNA sequencing in a set of post-mortem human donor pancreatic cells from four subjects and two mouse models to better understand the cellular diversity in the pancreas. Analysis of the data identified previously known cell types as well as rare and novel cell type subpopulations, and created a more detailed characterization of the diversity of those cell types. In this project, we will attempt to replicate their primary findings using current analytical methodology and software packages. 

Reference:
Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. 2016. “A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure.” Cell Systems 3 (4): 346–60.e4.

# Contributors

Monica Roberts: Data Curator (monicapr@bu.edu)

Preshita Dave: Analyst (preshita@bu.edu)

Italo Duran: Programmer (duran01@bu.edu)

# Repository Contents

1. Data Curator:

    barcode.py: Takes zipped fasta file as input and outputs txt file that contains unique barcodes with number of reads for each.
    
    combine.py: Takes text files of filtered barcodes from each sample and combines them into single list in text file for salmon alevin whitelist input.
    
    dist_plot.py: Plots cumulative distribution of reads per distinct barcode.
    
    filter_barcode.py: Takes barcode list text file, filters for barcodes that have > 1000 reads, and outputs list of filtered barcodes in txt file.
    
    index.qsub: Creates index of human transcriptome.
    
    salmon.qsub: Runs salmon alevin to map and count single-cell RNA-seq reads.
    
2. Programmer:

3. Analyst:

    analyst.R: This script analyses and identifies marker genes for clusters that were created. 

