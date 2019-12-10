#!/usr/bin/env bash

source ~/.bash_profile

dir_repo="/Users/kef74yk/Dropbox_w/repos"

cd ${dir_repo}/RADTE/test

# R -q -e 'devtools::install_local(path="/Users/kef74yk/Dropbox_w/repos/rkftools", reload=TRUE, quick=FALSE, local=TRUE, dep=FALSE)'

out_dir="../data/sample00/output"
#out_dir="../data/OG0000567/output"
#out_dir="../data/OG0001076/output"
#out_dir="../data/OG0002332/output"
#out_dir="../data/PLANT_OG0000099/output"
#out_dir="../data/PLANT_IDA/output"

if [ -e ${out_dir} ]; then
    rm -r ${out_dir}
fi
mkdir ${out_dir}
cd ${out_dir}

java -jar -Xmx2g /Applications/Notung-2.9.1.3.jar \
-s ../input/species_tree.nwk \
-g ../input/gene_tree.nwk \
--reconcile \
--infertransfers "false" \
--treeoutput newick \
--parsable \
--speciestag prefix \
--maxtrees 1 \
--nolosses \
--outputdir .

if [ -e "gene_tree.nwk.reconciled.0" ]; then
    notung_reconciled="gene_tree.nwk.reconciled.0"
    notung_parsable="gene_tree.nwk.reconciled.0.parsable.txt"
else
    notung_reconciled="gene_tree.nwk.reconciled"
    notung_parsable="gene_tree.nwk.reconciled.parsable.txt"
fi

../../../radte \
--species_tree=../input/species_tree.nwk \
--gene_tree=${notung_reconciled} \
--notung_parsable=${notung_parsable} \
--max_age=1000 \
--chronos_lambda=1 \
--chronos_model=discrete \
--pad_short_edge=0.001
