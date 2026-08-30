# R CMD check runs this against the installed package (not sourced modules).
library(RADTE)

local({
    dir = tempfile('radte-package-')
    dir.create(dir)
    on.exit(unlink(dir, recursive=TRUE))
    species = file.path(dir, 'species.nwk')
    gene = file.path(dir, 'gene.nwk')
    events = file.path(dir, 'events.txt')
    writeLines('((C_sp:7,D_sp:7)cd:3,(A_sp:5,B_sp:5)ab:5)root;', species)
    writeLines('((A_sp_g1:1,B_sp_g1:1)gn_ab:1,(C_sp_g1:1,D_sp_g1:1)gn_cd:1)gn_root;', gene)
    writeLines(character(), events)
    result = RADTE::radte_main(c(
        paste0('--species_tree=', species), paste0('--gene_tree=', gene),
        paste0('--notung_parsable=', events), paste0('--outdir=', dir),
        '--max_age=20', '--chronos_model=discrete', '--chronos_lambda=1'
    ))
    output = ape::read.tree(file.path(dir, 'radte_gene_tree_output.nwk'))
    input = ape::read.tree(gene)
    stopifnot(identical(output$edge, input$edge), identical(output$node.label, input$node.label))
    ages = max(ape::node.depth.edgelength(output)) - ape::node.depth.edgelength(output)
    used = read.delim(file.path(dir, 'radte_calibration_used.tsv'))
    stopifnot(all(abs(ages[used$node] - used$age.min) < 1e-7))
    manifest = read.delim(result$manifest, colClasses='character')
    stopifnot(manifest$value[manifest$key == 'run_status'] == 'complete')
})
