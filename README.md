## Overview
**R**econciliation-**A**ssisted **D**ivergence **T**ime **E**stimation (**RADTE**) is a method to date gene trees with the aid of dated species trees.
This program can handle a rooted gene tree containing duplication/loss events.
The divergence time of duplication nodes are estimated while constraining speciation nodes.
Currently, only point estimates are supported to constrain speciation events.
Confidence intervals of species divergence time will be incorporated in future. 
**RADTE** first attempts to constrain all available calibration points transferred from the species tree (**R**, root node; **D**, duplication node; **S**, speciation node) for the divergence time estimation by `chronos` from the **ape** package.
If the estimation fails, the constraints are gradually relaxed until successful estimation is obtained.
**RADTE** is under development. Any feedback is welcomed.

![](img/radte_method.svg)

## Dependency
* [R 3.x](https://www.r-project.org/)
* [ape](http://ape-package.ird.fr/)
* [treeio](https://github.com/YuLab-SMU/treeio)
* [rkftools](https://github.com/kfuku52/rkftools)

In addition to the dependencies above, RADTE needs an output from a phylogeny reconciliation program. Currently, **NOTUNG** and **GeneRax** are supported.
* [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/)
* [GeneRax](https://github.com/BenoitMorel/GeneRax)

## Installation
After installing the above dependencies, please download the `radte.r` script by, for example, `git` or `svn`, and change the file permission. You can also download a zipped repository by `Clone or dowonload` above.
```
# With git
git clone https://github.com/kfuku52/RADTE
cd RADTE

# With svn
svn export https://github.com/kfuku52/RADTE/trunk/radte.r
chmod +x ./radte.r
```

## Options
#### `--species_tree`
Species tree with estimated divergence time.
Leaves (species) should be labeled as `GENUS_SPECIES` (e.g., Homo_sapiens).
The tree is expected to be ultrametric and branch lengths should represent evolutionary time (e.g., million years).
Internal nodes must be uniquely labeled and the same file should be consistently used for **NOTUNG/GeneRax** and **RADTE**.
Don't know how to label internal nodes? Try this R one-liner.
```
R -q -e "library(ape); t=read.tree('species_tree_noLabel.nwk'); \
t[['node.label']]=paste0('s',1:Nnode(t)); \
write.tree(t, 'species_tree.nwk')"
```
#### `--gene_tree`
Rooted newick tree. Leaves (genes) should be labeled as `GENUS_SPECIES_GENEID` (e.g., Homo_sapiens_ENSG00000102144). The tree is expected to be non-ultrametric and branch lengths should represent substitutions per site. 
Use the tree that **NOTUNG** produces because its internal nodes are correctly labeled in accordance with the **NOTUNG parsable file**, another input for this program.
#### `--notung_parsable`
An output file from **NOTUNG** (tested with version 2.9) can be used to acquire the species–gene relationships in phylogeny reconciliation. See **Examples** for details.
#### `--generax_nhx`
Instead of the **NOTUNG** output, the NHX tree from **GeneRax** can also be used as an input. If specified, `--gene_tree` and `--notung_parsable` will be ignored. See **Examples** for details.
#### `--max_age`
If duplication nodes are deeper than the root node of the species tree, this value will be used as an upper limit of the root node.
#### `--chronos_lambda`
Passed to `chronos` for divergence time estimation. See `chronos` in the [**ape** documentation](https://www.rdocumentation.org/packages/ape/versions/5.2/topics/chronos).
#### `--chronos_model`
Passed to `chronos` for divergence time estimation. See `chronos` in the [**ape** documentation](https://www.rdocumentation.org/packages/ape/versions/5.2/topics/chronos).
#### `--pad_short_edge`
Prohibit dated branches shorter than this value. If detected, the branch length is readjusted by transferring a small portion of branch length from the parent branch.

## Example 1: RADTE after NOTUNG
```
# Run NOTUNG in the reconciliation mode
# Don't forget to specify --parsable
java -jar -Xmx2g Notung-2.9.jar \
-s species_tree.nwk \
-g gene_tree_input.nwk \
--reconcile \
--infertransfers "false" \
--treeoutput newick \
--parsable \
--speciestag prefix \
--maxtrees 1 \
--nolosses \
--outputdir .

# Run RADTE
./radte.r \
--species_tree=species_tree.nwk \
--gene_tree=gene_tree_input.nwk.reconciled \
--notung_parsable=gene_tree_input.nwk.reconciled.parsable.txt \
--max_age=1000 \
--chronos_lambda=1 \
--chronos_model=difscrete \
--pad_short_edge=0.001
```
#### species_tree.nwk
![](img/notung_radte_species_tree.svg)

#### gene_tree_input.nwk.reconciled
![](img/notung_radte_gene_tree_input.svg)

#### gene_tree_output.nwk
![](img/notung_radte_gene_tree_output.svg)

## Example 2: RADTE after GeneRax
```
./radte.r \
--species_tree=species_tree.nwk \
--generax_nhx=gene_tree_input.nhx \
--max_age=1000 \
--chronos_lambda=1 \
--chronos_model=discrete \
--pad_short_edge=0.001
```

#### species_tree.nwk
![](img/generax_radte_species_tree.svg)

#### gene_tree_input.nhx
![](img/generax_radte_gene_tree_input.svg)

#### gene_tree_output.nwk
![](img/generax_radte_gene_tree_output.svg)

## Update
Although the above tree images are not painted, the latest version of RADTE generates pdf files in which branches are colored. Red means the branches connected to the unconstrained nodes, whereas blue indicates the branches connected to the constrained nodes.

## Citation
Fukushima K, Pollock DD. 2020. Amalgamated cross-species transcriptomes reveal organ-specific propensity in gene expression evolution. **Nature Communications 11**: 4459 ([DOI: 10.1038/s41467-020-18090-8]( https://doi.org/10.1038/s41467-020-18090-8 ))

## Licensing
This program is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.
