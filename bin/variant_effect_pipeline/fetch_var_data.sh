#!/usr/bin/env bash
# Script to extract specific accesions from SIFT and FoldX bulk data

# Define variables
species=$1 # one of yeast/human/ecoli
name=$2 # gene name in output files
acc=$3 # uniprot accession to search for
proj_dir=$HOME/Projects/mutations/data

case $species in
    yeast)
        sift_file=$proj_dir/mutfunc/yeast/conservation/sift_filtered.tab
        foldx_file=$proj_dir/mutfunc/yeast/structure/mod_ddg1.tab
        ;;
    human)
        sift_file=$proj_dir/mutfunc/human/conservation/sift_parsed_5.tab
        foldx_file=$proj_dir/mutfunc/human/structure/exp_ddg1.tab
        ;;
    ecoli)
        sift_file=$proj_dir/mutfunc/ecoli/conservation/sift_filtered.tab
        foldx_file=$proj_dir/mutfunc/ecoli/structure/mod_ddg1.tab
        ;;
    *)
        echo "Incorrect species ($1) given, must be one of yeasy/human/ecoli"
        exit -1
        ;;
esac

# Extract sift
cat <(head -n 1 $sift_file) <(grep $acc $sift_file) > $proj_dir/$species'_'$name'_sift.tsv'

# Extract FoldX
cat <(head -n 1 $foldx_file) <(grep $acc $foldx_file) > $proj_dir/$species'_'$name'_foldx.tsv'
