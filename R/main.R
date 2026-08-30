radte_main <- function(command_args = commandArgs(trailingOnly = TRUE)) {
    if (length(command_args) == 1L && command_args[[1]] %in% c('--help', '-h')) {
        cat(render_radte_help(), '\n')
        return(invisible(0L))
    }
    if (length(command_args) == 1L && command_args[[1]] %in% c('--version', '-V')) {
        cat(radte_version, '\n', sep='')
        return(invisible(0L))
    }
    if (!requireNamespace('ape', quietly=TRUE)) {
        stop('The R package "ape" is required. Install it before running RADTE.')
    }
    run_started_at = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z')

    cat(paste(version[['version.string']], '\n'))
    cat(paste('ape version:', utils::packageVersion('ape'), '\n'))
    cat(paste('RADTE version:', radte_version, '\n'))
    cat(paste('RADTE command:', paste(c('radte', command_args), collapse=' '), '\n'))
    cat(paste('RADTE bug report:', 'https://github.com/kfuku52/RADTE/issues', '\n'))

    cat('arguments:\n')
    args = resolve_radte_config(parse_radte_cli_args(command_args, print=TRUE))
    mode = attr(args, 'mode')
    dating_backend = args$dating_backend
    outdir = args$outdir
    dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
    outdir = normalizePath(outdir, mustWork=TRUE)
    args$outdir = outdir
    prefix = args$prefix
    if (dating_backend == 'mcmctree' && is.null(args$mcmctree_workdir)) {
        args$mcmctree_workdir = tempfile(pattern=paste0(prefix, '_mcmctree_run_'), tmpdir=outdir)
    }
    input_names = c('species_tree', 'gene_tree', 'notung_parsable', 'generax_nhx',
                    'species_map_tsv', 'species_node_bounds_tsv', 'mcmctree_seqfile')
    input_files = args[intersect(input_names, names(args))]
    attr(args, "input_files") = input_files
    input_hashes = lapply(input_files, compute_file_sha256)
    names(input_hashes) = paste0('input_sha256.', names(input_hashes))
    transaction = begin_output_transaction(outdir, prefix, unlist(input_files, use.names=FALSE))
    on.exit(cleanup_output_transaction(transaction), add=TRUE)
    species_parser = build_species_parser(args$species_parser, args$species_regex, args$species_map_tsv)
    cat('dating backend:', dating_backend, '\nspecies parser:', species_parser$name, '\n')
    species_context = read_radte_species_context(args, species_parser)
    sp_tree = species_context$tree
    sp_node_table = species_context$nodes
    gene_context = read_radte_gene_context(args, mode, sp_tree, sp_node_table, species_parser)
    calibrations = prepare_radte_calibrations(gene_context)
    result = run_radte_dating(args, gene_context, species_context, species_parser, calibrations$tables)
    write_radte_outputs(args, gene_context, species_context, calibrations$tables,
                        calibrations$all, result, transaction, run_started_at, command_args, input_hashes)
}
