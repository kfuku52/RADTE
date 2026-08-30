write_radte_outputs = function(args, gene_context, species_context, calibration_tables,
    calibration_table, result, transaction, run_started_at, command_args, input_hashes) {
    gn_tree = gene_context$tree
    gn_node_table = gene_context$nodes
    sp_tree = species_context$tree
    sp_node_table = species_context$nodes
    root_num = get_root_num(gn_tree)
    calibration_table_R = calibration_tables$R
    chronos_out = result$tree
    current_calibration_table = result$calibration
    calibrated_node = result$calibrated_nodes
    mcmctree_mirror_table = result$mirror_table
    shared_speciation_age_summary = result$shared_ages
    dating_backend = args$dating_backend
    outdir = args$outdir
    prefix = args$prefix
    output_path = function(suffix) file.path(transaction$stage, paste0(prefix, suffix))
    cat('Writing output files.\n')
    chronos_out2 = chronos_out
    padding_applied = FALSE
    if (dating_backend == 'chronos' && !is.null(args$pad_short_edge) &&
        any(chronos_out2$edge.length < args$pad_short_edge)) {
        chronos_out2 = build_dated_tree_without_chronos(
            chronos_out2, current_calibration_table, root_num, min_edge=args$pad_short_edge, require_root=FALSE)
        if (inherits(chronos_out2, 'try-error')) stop(as.character(chronos_out2))
        padding_applied = TRUE
    }
    validate_dated_tree(chronos_out2, gn_tree,
        calibration=if (dating_backend == 'chronos') current_calibration_table else NULL,
        tolerance=if (dating_backend == 'mcmctree') 1e-5 else 1e-7)

    atomic_write_lines(calibrated_node, file=output_path('_calibrated_nodes.txt'))

    atomic_write_tree(chronos_out2, file=output_path('_gene_tree_output.nwk'))
    emitted_calibration_nodes = as.integer(current_calibration_table[['node']])
    current_calibration_table = annotate_calibration_output(
        calibration_table=current_calibration_table,
        gn_node_table=gn_node_table,
        mirror_table=mcmctree_mirror_table,
        emitted_nodes=emitted_calibration_nodes
    )
    current_calibration_table[current_calibration_table$node==calibration_table_R$node,'event'] = 'R'
    current_calibration_table = annotate_calibration_fit(current_calibration_table, chronos_out2, dating_backend, emitted_calibration_nodes)
    atomic_write_table(
        current_calibration_table,
        file=output_path('_calibration_used.tsv'),
        sep='\t', quote=FALSE, row.names=FALSE
    )

    calibration_table = annotate_calibration_output(
        calibration_table=calibration_table,
        gn_node_table=gn_node_table,
        mirror_table=mcmctree_mirror_table,
        emitted_nodes=emitted_calibration_nodes
    )
    calibration_table = annotate_calibration_fit(calibration_table, chronos_out2, dating_backend, emitted_calibration_nodes)
    atomic_write_table(
        calibration_table,
        file=output_path('_calibration_all.tsv'),
        sep='\t', quote=FALSE, row.names=FALSE
    )

    if (dating_backend == 'mcmctree') {
        atomic_write_table(
            shared_speciation_age_summary,
            file=output_path('_shared_speciation_ages.tsv'),
            sep='\t',
            quote=FALSE,
            row.names=FALSE
        )
    }

    gn_node_table$spp = NULL
    atomic_write_table(
        gn_node_table,
        file=output_path('_gene_tree.tsv'),
        sep='\t', quote=FALSE, row.names=FALSE
    )

    sp_node_table$spp = NULL
    atomic_write_table(
        sp_node_table,
        file=output_path('_species_tree.tsv'),
        sep='\t', quote=FALSE, row.names=FALSE
    )

    if ((dating_backend == 'mcmctree') && (!is.na(result$metadata$mcmctree_workdir))) {
        mcmctree_artifacts = c(
            'mcmctree.ctl',
            'input.trees',
            'out.txt',
            'mcmc.txt',
            'FigTree.tre',
            'mcmctree.stdout.log',
            'mcmctree.stderr.log'
        )
        for (artifact in mcmctree_artifacts) {
            source_artifact = file.path(result$metadata$mcmctree_workdir, artifact)
            if (file.exists(source_artifact)) {
                target_artifact = output_path(paste0('_mcmctree_', artifact))
                atomic_copy_file(source_artifact, target_artifact)
            }
        }
    }

    node_nums = (length(chronos_out2[['tip.label']])+1):max(chronos_out2[['edge']])
    noncalibrated_nodes = node_nums[!node_nums %in% current_calibration_table[['node']]]
    ec = list('red'=noncalibrated_nodes, 'blue'=current_calibration_table[['node']])
    atomic_save_tree_pdf(phy=gn_tree, file=output_path('_gene_tree_input.pdf'), show.age=FALSE, edge_colors=ec)
    atomic_save_tree_pdf(phy=chronos_out2, file=output_path('_gene_tree_output.pdf'), show.age=TRUE, edge_colors=ec)
    atomic_save_tree_pdf(phy=sp_tree, file=output_path('_species_tree.pdf'), show.age=TRUE)

    option_entries = lapply(args, function(value) {
        if (is.infinite(value)) 'inf' else value
    })
    names(option_entries) = paste0('option.', names(option_entries))
    manifest_entries = c(
        list(
            radte_version=radte_version,
            run_id=transaction$run_id,
            run_status="complete",
            rng_kind=paste(RNGkind(), collapse=", "),
            branch_padding_applied=padding_applied,
            calibration_policy=if (dating_backend == "chronos") "hard bounds retained; drops explicit" else "PAML soft priors; structural validation",
            run_started_at=run_started_at,
            run_completed_at=format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'),
            command=paste(c('radte', command_args), collapse=' '),
            r_version=R.version.string,
            ape_version=as.character(utils::packageVersion('ape')),
            treeio_version=if (requireNamespace('treeio', quietly=TRUE)) as.character(utils::packageVersion('treeio')) else NA_character_,
            output_directory=outdir,
            output_prefix=prefix
        ),
        result$metadata,
        option_entries,
        input_hashes
    )
    manifest_file = output_path('_run_manifest.tsv')
    artifacts = list.files(transaction$stage, full.names=TRUE)
    output_hashes = lapply(artifacts, compute_file_sha256)
    names(output_hashes) = paste0('output_sha256.', basename(artifacts))
    manifest_entries = c(manifest_entries, output_hashes)
    write_run_manifest(manifest_entries, manifest_file)
    commit_output_transaction(transaction, basename(manifest_file))
    manifest_file = file.path(outdir, basename(manifest_file))

    cat('dating backend used:', dating_backend, '\n')
    if (!is.na(result$metadata$chronos_model_used)) {
        cat('chronos model used:', result$metadata$chronos_model_used, '\n')
    }
    if (!is.na(result$metadata$chronos_lambda_used)) {
        cat('chronos lambda used:', result$metadata$chronos_lambda_used, '\n')
    }
    if (!is.na(result$metadata$chronos_seed_used)) {
        cat('chronos seed used:', result$metadata$chronos_seed_used, '\n')
    }
    if (!is.na(result$metadata$mcmctree_workdir)) {
        cat('MCMCTree workdir:', result$metadata$mcmctree_workdir, '\n')
    }
    cat('Calibrated nodes:', calibrated_node, '\n')
    cat('Tree height:', max(ape::node.depth.edgelength(sp_tree)), '\n')
    constrained_sp_nodes = unique(gn_node_table[
        grepl('^S', gn_node_table[['event']]) &
        !is.na(gn_node_table[['constraint_sp_node']]) &
        (gn_node_table[['constraint_sp_node']] != ''),
        'constraint_sp_node'
    ])
    num_spnode_used_for_constraint = length(constrained_sp_nodes)
    cat('Number of species tree node used for the gene tree constraint:', num_spnode_used_for_constraint, '\n')
    cat('Completed: RADTE divergence time estimation', '\n')
    return(invisible(list(
        output_directory=outdir,
        output_prefix=prefix,
        manifest=manifest_file,
        dating_backend=dating_backend,
        calibrated_nodes=calibrated_node
    )))
}
