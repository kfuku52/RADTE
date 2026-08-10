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
    args = parse_radte_cli_args(command_args, print=TRUE)

    outdir = if ('outdir' %in% names(args)) args[['outdir']] else '.'
    if (!is_nonempty_scalar_string(outdir)) {
        stop('--outdir should be a non-empty path.')
    }
    dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
    outdir = normalizePath(outdir, mustWork=TRUE)
    prefix = if ('prefix' %in% names(args)) args[['prefix']] else 'radte'
    if (!is_nonempty_scalar_string(prefix) || grepl('[/\\\\]', prefix)) {
        stop('--prefix should be a non-empty filename prefix without path separators.')
    }
    random_seed = if ('seed' %in% names(args)) args[['seed']] else 1L
    output_path = function(suffix) file.path(outdir, paste0(prefix, suffix))

    if (('generax_nhx' %in% names(args))&('notung_parsable' %in% names(args))) {
        stop('Only one of --notung_parsable and --generax_nhx should be specified. Exiting.\n')
    } else if ('generax_nhx' %in% names(args)) {
        cat('--generax_nhx was detected. GeneRax mode.', '\n')
        mode = 'generax'
        generax_file = args[['generax_nhx']]
    } else if ('notung_parsable' %in% names(args)) {
        cat('--notung_parsable was detected. Notung mode.', '\n')
        mode = 'notung'
        gn_file = args[['gene_tree']]
        parsable_file = args[['notung_parsable']]
    } else {
        stop('--notung_parsable or --generax_nhx should be specified. Exiting.\n')
    }

    dating_backend = 'chronos'
    if ('dating_backend' %in% names(args)) {
        dating_backend = parse_choice_arg(args[['dating_backend']], '--dating_backend', c('chronos', 'mcmctree'))
    }
    cat('dating backend:', dating_backend, '\n')

    required_args = c('species_tree', 'max_age')
    if (dating_backend == 'chronos') {
        required_args = c(required_args, 'chronos_lambda', 'chronos_model')
    } else if (dating_backend == 'mcmctree') {
        required_args = c(required_args, 'mcmctree_seqfile')
    }
    missing_required_args = required_args[!required_args %in% names(args)]
    if (length(missing_required_args) > 0) {
        stop(
            'Missing required argument(s): ',
            paste(paste0('--', missing_required_args), collapse=', ')
        )
    }
    if ((mode=='notung') && !('gene_tree' %in% names(args))) {
        stop('--gene_tree should be specified in Notung mode.')
    }

    sp_file = args[['species_tree']]
    max_age = suppressWarnings(as.numeric(args[['max_age']]))
    if (is.na(max_age) || (!is.finite(max_age)) || (max_age <= 0)) {
        stop('--max_age should be a positive finite number.')
    }
    chronos_lambda = NA_real_
    chronos_model = NA_character_
    chronos_iter_max = 250L
    chronos_eval_max = 250L
    chronos_dual_iter_max = 20L
    chronos_high_control_fallback = TRUE
    if (dating_backend == 'chronos') {
        chronos_lambda = suppressWarnings(as.numeric(args[['chronos_lambda']]))
        if (is.na(chronos_lambda) || (!is.finite(chronos_lambda)) || (chronos_lambda < 0)) {
            stop('--chronos_lambda should be a non-negative finite number.')
        }
        chronos_model = args[['chronos_model']]
        if ((!is.character(chronos_model)) || (length(chronos_model) != 1) || is.na(chronos_model) || (nchar(chronos_model)==0)) {
            stop('--chronos_model should be a non-empty string.')
        }
        supported_chronos_models = c('discrete', 'relaxed', 'correlated')
        if (!chronos_model %in% supported_chronos_models) {
            suggestion = ''
            if (tolower(chronos_model)=='difscrete') {
                suggestion = ' Did you mean "discrete"?'
            }
            stop(
                '--chronos_model should be one of: ',
                paste(supported_chronos_models, collapse=', '),
                '. Received: ',
                chronos_model,
                '.',
                suggestion
            )
        }
        if ('chronos_iter_max' %in% names(args)) {
            chronos_iter_max = parse_integer_arg(args[['chronos_iter_max']], '--chronos_iter_max', min_value=1L)
        }
        if ('chronos_eval_max' %in% names(args)) {
            chronos_eval_max = parse_integer_arg(args[['chronos_eval_max']], '--chronos_eval_max', min_value=1L)
        }
        if ('chronos_dual_iter_max' %in% names(args)) {
            chronos_dual_iter_max = parse_integer_arg(
                args[['chronos_dual_iter_max']],
                '--chronos_dual_iter_max',
                min_value=0L
            )
        }
        if ('chronos_high_control_fallback' %in% names(args)) {
            chronos_high_control_fallback = parse_bool_arg(
                args[['chronos_high_control_fallback']],
                '--chronos_high_control_fallback'
            )
        }
    }
    if ('pad_short_edge' %in% names(args)) {
        pad_short_edge = suppressWarnings(as.numeric(args[['pad_short_edge']]))
        if (is.na(pad_short_edge) || (!is.finite(pad_short_edge)) || (pad_short_edge <= 0)) {
            stop('--pad_short_edge should be a positive finite number.')
        }
        args[['pad_short_edge']] = pad_short_edge
    }
    allow_constraint_drop = TRUE
    if ('allow_constraint_drop' %in% names(args)) {
        allow_constraint_drop = parse_bool_arg(args[['allow_constraint_drop']], '--allow_constraint_drop')
    }
    chronos_attempt_timeout_sec = 60
    chronos_total_timeout_sec = 300
    if (dating_backend == 'chronos') {
        if ('chronos_attempt_timeout_sec' %in% names(args)) {
            chronos_attempt_timeout_sec = parse_timeout_arg(args[['chronos_attempt_timeout_sec']], '--chronos_attempt_timeout_sec')
        }
        if ('chronos_total_timeout_sec' %in% names(args)) {
            chronos_total_timeout_sec = parse_timeout_arg(args[['chronos_total_timeout_sec']], '--chronos_total_timeout_sec')
        }
        if (is.finite(chronos_total_timeout_sec) && is.finite(chronos_attempt_timeout_sec) && (chronos_total_timeout_sec < chronos_attempt_timeout_sec)) {
            cat(
                'Adjusting chronos attempt timeout from',
                chronos_attempt_timeout_sec,
                'to',
                chronos_total_timeout_sec,
                'sec to respect total time budget.\n'
            )
            chronos_attempt_timeout_sec = chronos_total_timeout_sec
        }
        chronos_timeout_label = if (is.finite(chronos_attempt_timeout_sec)) chronos_attempt_timeout_sec else 'inf'
        chronos_budget_label = if (is.finite(chronos_total_timeout_sec)) chronos_total_timeout_sec else 'inf'
        cat('chronos timeout settings: attempt_sec=', chronos_timeout_label, ', total_sec=', chronos_budget_label, '\n', sep='')
    }

    mcmctree_seqfile = NULL
    mcmctree_bin = 'mcmctree'
    mcmctree_workdir = NULL
    mcmctree_timeout_sec = Inf
    mcmctree_usedata = 1L
    mcmctree_seqtype = 0L
    mcmctree_clock = 2L
    mcmctree_model = 0L
    mcmctree_burnin = 2000L
    mcmctree_sampfreq = 10L
    mcmctree_nsample = 20000L
    mcmctree_ncatG = 5L
    if (dating_backend == 'mcmctree') {
        mcmctree_seqfile = as.character(args[['mcmctree_seqfile']])
        if (!is_nonempty_scalar_string(mcmctree_seqfile)) {
            stop('--mcmctree_seqfile should be a non-empty path.')
        }
        if ('mcmctree_bin' %in% names(args)) {
            mcmctree_bin = as.character(args[['mcmctree_bin']])
        }
        if ('mcmctree_workdir' %in% names(args)) {
            mcmctree_workdir = as.character(args[['mcmctree_workdir']])
            if (!is_nonempty_scalar_string(mcmctree_workdir)) {
                stop('--mcmctree_workdir should be a non-empty path.')
            }
        }
        if (is.null(mcmctree_workdir)) {
            run_stamp = paste0(format(Sys.time(), '%Y%m%dT%H%M%S'), '_', Sys.getpid())
            mcmctree_workdir = file.path(outdir, paste0(prefix, '_mcmctree_run_', run_stamp))
        }
        if ('mcmctree_timeout_sec' %in% names(args)) {
            mcmctree_timeout_sec = args[['mcmctree_timeout_sec']]
        }
        if ('mcmctree_usedata' %in% names(args)) {
            mcmctree_usedata = parse_integer_arg(args[['mcmctree_usedata']], '--mcmctree_usedata', min_value=1)
        }
        if (!mcmctree_usedata %in% c(1L)) {
            stop('--mcmctree_usedata currently supports only 1 in RADTE.')
        }
        if ('mcmctree_seqtype' %in% names(args)) {
            mcmctree_seqtype = parse_integer_arg(args[['mcmctree_seqtype']], '--mcmctree_seqtype', min_value=0)
        }
        if ('mcmctree_clock' %in% names(args)) {
            mcmctree_clock = parse_integer_arg(args[['mcmctree_clock']], '--mcmctree_clock', min_value=1)
        }
        if ('mcmctree_model' %in% names(args)) {
            mcmctree_model = parse_integer_arg(args[['mcmctree_model']], '--mcmctree_model', min_value=0)
        }
        if ('mcmctree_burnin' %in% names(args)) {
            mcmctree_burnin = parse_integer_arg(args[['mcmctree_burnin']], '--mcmctree_burnin', min_value=1)
        }
        if ('mcmctree_sampfreq' %in% names(args)) {
            mcmctree_sampfreq = parse_integer_arg(args[['mcmctree_sampfreq']], '--mcmctree_sampfreq', min_value=1)
        }
        if ('mcmctree_nsample' %in% names(args)) {
            mcmctree_nsample = parse_integer_arg(args[['mcmctree_nsample']], '--mcmctree_nsample', min_value=1)
        }
        if ('mcmctree_ncatG' %in% names(args)) {
            mcmctree_ncatG = parse_integer_arg(args[['mcmctree_ncatG']], '--mcmctree_ncatG', min_value=1)
        }
        cat(
            'MCMCTree settings: usedata=', mcmctree_usedata,
            ', seqtype=', mcmctree_seqtype,
            ', clock=', mcmctree_clock,
            ', model=', mcmctree_model,
            ', burnin=', mcmctree_burnin,
            ', sampfreq=', mcmctree_sampfreq,
            ', nsample=', mcmctree_nsample,
            ', ncatG=', mcmctree_ncatG,
            ', seed=', random_seed,
            ', timeout_sec=', if (is.finite(mcmctree_timeout_sec)) mcmctree_timeout_sec else 'inf',
            '\n',
            sep=''
        )
    }

    species_parser_name = 'legacy'
    if ('species_parser' %in% names(args)) {
        species_parser_name = args[['species_parser']]
    }
    species_regex = NULL
    if ('species_regex' %in% names(args)) {
        species_regex = args[['species_regex']]
    }
    species_map_tsv = NULL
    if ('species_map_tsv' %in% names(args)) {
        species_map_tsv = args[['species_map_tsv']]
    }
    species_node_bounds_tsv = NULL
    if ('species_node_bounds_tsv' %in% names(args)) {
        species_node_bounds_tsv = as.character(args[['species_node_bounds_tsv']])
        if (!is_nonempty_scalar_string(species_node_bounds_tsv)) {
            stop('--species_node_bounds_tsv should be a non-empty path.')
        }
    }
    species_parser = build_species_parser(
        parser_name=species_parser_name,
        species_regex=species_regex,
        species_map_tsv=species_map_tsv
    )
    cat('species parser:', species_parser[['name']], '\n')

    cat('\nStart: species tree processing', '\n')
    tree_text0 = scan(sp_file, what=character(), sep="\n", blank.lines.skip=FALSE)
    tree_text1 = gsub("'([0-9]+)'", "PLACEHOLDER\\1", tree_text0)
    sp_tree = ape::read.tree(text=tree_text1)
    if (!is.null(sp_tree[['node.label']])) {
        sp_tree[['node.label']] = sub('^PLACEHOLDER', '', sp_tree[['node.label']])
    }
    sp_tree = validate_species_tree_labels(sp_tree)
    sp_tip_species = validate_species_tip_parser_labels(sp_tree, species_parser)
    if (contains_polytomy(sp_tree)) {
        stop('Input species tree contains polytomy. A completely bifurcated tree is expected as input.')
    }
    validate_tree_edge_lengths(sp_tree, 'species tree')
    if (length(args[['pad_short_edge']])) {
        sp_tree = pad_short_edges(sp_tree, threshold=args[['pad_short_edge']], external_only=FALSE)
    }
    sp_tree = force_ultrametric(sp_tree, stop_if_larger_change=0.01)
    sp_node_depths = ape::node.depth.edgelength(sp_tree)
    root_depth = max(sp_node_depths)
    sp_node_ages = abs(sp_node_depths - root_depth)
    sp_node_names = c(sp_tree[['tip.label']], sp_tree[['node.label']])
    sp_node_table = data.frame(
        node=sp_node_names,
        age=sp_node_ages,
        age_min=sp_node_ages,
        age_max=sp_node_ages,
        spp=NA,
        stringsAsFactors=FALSE
    )
    sp_tip_rows = seq_len(ape::Ntip(sp_tree))
    max_tip_age = max(sp_node_table[sp_tip_rows,'age'])
    if (max_tip_age!=0) {
        cat(paste0('Nonzero tip age(s) were detected (max=', max_tip_age, '). Coercing to 0.\n'))
        sp_node_table[sp_tip_rows,'age'] = 0
        sp_node_table[sp_tip_rows,'age_min'] = 0
        sp_node_table[sp_tip_rows,'age_max'] = 0
    }
    if (!is.null(species_node_bounds_tsv)) {
        cat('Reading species node bounds TSV.\n')
        species_node_bounds = read_species_node_bounds_tsv(species_node_bounds_tsv)
        sp_node_table = merge_species_node_bounds(sp_node_table, species_node_bounds)
    }
    has_species_node_interval_input = any(abs(sp_node_table[['age_max']] - sp_node_table[['age_min']]) > 1e-12)
    cat('End: species tree processing', '\n\n')

    cat('Start: gene tree processing', '\n')
    if (mode=='generax') {
        cat('Reading GeneRax tree.\n')
        nhxtree = read_generax_nhx(generax_file)

        gn_tree = nhxtree@phylo
        if (contains_polytomy(gn_tree)) {
            stop('Input tree contains polytomy. A completely bifurcated tree is expected as input.')
        }
        gn_tree = validate_gene_tree_labels(gn_tree)
        validate_gene_tip_labels(gn_tree[['tip.label']], species_parser=species_parser, species_tree_labels=sp_tree[['tip.label']])
        validate_gene_species_coverage(gn_tree, sp_tree, species_parser=species_parser, species_tree_labels=sp_tree[['tip.label']])
        gn_tree = pad_branch_length(gn_tree, pad_size=args[['pad_short_edge']])
        #gn_tree = adjust_branch_length_order(gn_tree, min_bl=args[['pad_short_edge']])
        cat('Minimum branch length in gene tree:', min(gn_tree[['edge.length']]), '\n')
        cols = c('event', 'gn_node', 'gn_node_num', 'lower_sp_node', 'upper_sp_node', 'lower_age', 'upper_age')
        gn_node_table = data.frame(nhxtree@data, stringsAsFactors=FALSE)
        expected_node_rows = length(gn_tree[['tip.label']]) + gn_tree[['Nnode']]
        if (nrow(gn_node_table) != expected_node_rows) {
            stop(
                'GeneRax NHX metadata row count does not match the number of nodes in the gene tree. ',
                'Expected ', expected_node_rows, ' row(s) but found ', nrow(gn_node_table), '.'
            )
        }
        if (!('S' %in% colnames(gn_node_table))) {
            stop('GeneRax NHX metadata is missing required species annotation tag: S')
        }
        if (any(is.na(gn_node_table[['S']]) | (gn_node_table[['S']] == ''))) {
            stop('GeneRax NHX metadata contains missing value(s) in required species annotation tag: S')
        }
        if (!('D' %in% colnames(gn_node_table))) {
            gn_node_table[['D']] = 'N'
        }
        gn_node_table[,'event'] = 'S'
        d_raw = gn_node_table[['D']]
        d_chr = toupper(trimws(as.character(d_raw)))
        d_chr[is.na(d_raw)] = 'N'
        is_missing_d = d_chr == ''
        d_chr[is_missing_d] = 'N'
        dup_markers = c('Y', 'YES', 'TRUE', 'T', '1')
        nondup_markers = c('N', 'NO', 'FALSE', 'F', '0')
        valid_markers = c(dup_markers, nondup_markers)
        unknown_markers = unique(d_chr[!d_chr %in% valid_markers])
        if (length(unknown_markers) > 0) {
            stop(
                'GeneRax NHX metadata contains unsupported duplication tag value(s): ',
                paste(unknown_markers, collapse=', '),
                '. Supported values are: ',
                paste(valid_markers, collapse=', ')
            )
        }
        is_dup = d_chr %in% dup_markers
        gn_node_table[is_dup,'event'] = 'D'
        colnames(gn_node_table) = sub('^S$', 'lower_sp_node', colnames(gn_node_table))
        gn_node_table[,'upper_sp_node'] = gn_node_table[['lower_sp_node']]
        gn_node_table = gn_node_table[order(gn_node_table[['node']]),]
        gn_node_table[,'gn_node'] = c(gn_tree[['tip.label']], gn_tree[['node.label']])
        duplication_rows = which(gn_node_table[['event']]=='D')
        gn_node_table[duplication_rows,'upper_sp_node'] = NA_character_
        if (length(duplication_rows) > 0) {
            lower_node_nums = match(
                gn_node_table[duplication_rows,'lower_sp_node'],
                sp_node_table[['node']]
            )
            sp_parent_map = get_parent_map(sp_tree)
            parent_nums = sp_parent_map[lower_node_nums]
            parent_names = rep(NA_character_, length(parent_nums))
            has_parent = !is.na(parent_nums)
            parent_names[has_parent] = sp_node_table[['node']][parent_nums[has_parent]]
            gn_node_table[duplication_rows,'upper_sp_node'] = parent_names
        }
        lower_sp_rows = match(gn_node_table[['lower_sp_node']], sp_node_table[['node']])
        is_exact_species_node = (!is.na(lower_sp_rows)) &
            (!is.na(gn_node_table[['upper_sp_node']])) &
            (gn_node_table[['lower_sp_node']] == gn_node_table[['upper_sp_node']])
        gn_node_table[,'lower_age'] = NA_real_
        gn_node_table[,'upper_age'] = NA_real_
        gn_node_table[is_exact_species_node,'lower_age'] = sp_node_table[['age_min']][lower_sp_rows[is_exact_species_node]]
        gn_node_table[is_exact_species_node,'upper_age'] = sp_node_table[['age_max']][lower_sp_rows[is_exact_species_node]]
        gn_node_table[,'gn_node_num'] = get_node_num_by_name(gn_tree, gn_node_table[['gn_node']])
        gn_node_table = data.frame(gn_node_table[,cols], stringsAsFactors=FALSE)
        validate_gn_species_nodes(gn_node_table, sp_tree)
    } else if (mode=='notung') {
        cat('Reading NOTUNG tree.\n')
        gn_tree = ape::read.tree(gn_file)
        gn_tree[['node.label']] = gsub("\\'", "",gn_tree[['node.label']])
        if (contains_polytomy(gn_tree)) {
            stop('Input tree contains polytomy. A completely bifurcated tree is expected as input.')
        }
        gn_tree = validate_gene_tree_labels(gn_tree)
        validate_gene_tip_labels(gn_tree[['tip.label']], species_parser=species_parser, species_tree_labels=sp_tree[['tip.label']])
        validate_gene_species_coverage(gn_tree, sp_tree, species_parser=species_parser, species_tree_labels=sp_tree[['tip.label']])
        gn_tree = pad_branch_length(gn_tree, pad_size=args[['pad_short_edge']])

        gn_node_table = read_notung_parsable(file=parsable_file, mode='D')
        gn_node_table = merge(gn_node_table, data.frame(lower_age=NA, upper_age=NA, spp=NA), all=TRUE)
        check_gn_node_name_uniqueness(gn_node_table, gn_tree)
        if (nrow(gn_node_table) > 0) {
            gn_node_table$gn_node_num = get_node_num_by_name(
                gn_tree,
                gn_node_table[['gn_node']]
            )
            lower_sp_rows = match(gn_node_table[['lower_sp_node']], sp_node_table[['node']])
            upper_sp_rows = match(gn_node_table[['upper_sp_node']], sp_node_table[['node']])
            gn_node_table[['lower_age']] = sp_node_table[['age_min']][lower_sp_rows]
            gn_node_table[['upper_age']] = sp_node_table[['age_max']][upper_sp_rows]
        } else {
            gn_node_table = gn_node_table[0,]
        }

        gene_species_nodes = map_gene_nodes_to_species_nodes(
            gn_tree=gn_tree,
            sp_tree=sp_tree,
            species_parser=species_parser,
            species_tree_labels=sp_tree[['tip.label']]
        )
        internal_gene_nodes = seq.int(ape::Ntip(gn_tree) + 1L, ape::Ntip(gn_tree) + gn_tree[['Nnode']])
        speciation_gene_nodes = setdiff(internal_gene_nodes, gn_node_table[['gn_node_num']])
        sp_parent_map = get_parent_map(sp_tree)
        for (root_num in speciation_gene_nodes) {
            root_node = get_node_name_by_num(gn_tree, root_num)
            sp_node_num = gene_species_nodes[[root_num]]
            # Historical RADTE behavior maps a single-species internal gene
            # node to the smallest internal species-tree clade, not to a tip.
            if (sp_node_num <= ape::Ntip(sp_tree)) {
                sp_node_num = sp_parent_map[[sp_node_num]]
            }
            sp_node = get_node_name_by_num(sp_tree, sp_node_num)
            if (length(sp_node) != 1 || is.na(sp_node)) {
                stop('Could not map gene node to a species-tree clade: ', root_node)
            }
            ind = nrow(gn_node_table)+1
            gn_node_table[ind,'event'] = "S"
            gn_node_table[ind,'gn_node'] = root_node
            gn_node_table[ind,'gn_node_num'] = root_num
            gn_node_table[ind,'lower_sp_node'] = sp_node
            gn_node_table[ind,'upper_sp_node'] = sp_node
            gn_node_table[ind,'lower_age'] = sp_node_table[['age_min']][[sp_node_num]]
            gn_node_table[ind,'upper_age'] = sp_node_table[['age_max']][[sp_node_num]]
            gn_node_table[ind,'spp'] = NA_character_
        }
    }

    cat('End: gene tree processing', '\n\n')
    gn_node_table = annotate_species_constraint_groups(gn_node_table, gn_tree)
    has_shared_interval_speciation_input = gn_node_table_has_shared_interval_speciation(gn_node_table)
    validate_gn_node_table(gn_node_table)
    validate_duplication_nodes_internal(gn_node_table, gn_tree)

    # Calibration node check
    if ((sum(gn_node_table[['event']]=="D") > 0)&(any(is.na(gn_node_table[['upper_age']])))) {
        gn_spp = unique(
            get_species_names(
                phy=gn_tree,
                species_parser=species_parser,
                species_tree_labels=sp_tree[['tip.label']]
            )
        )
        gn_sp_tips = resolve_species_tree_tips(species_parser, gn_spp, sp_tree[['tip.label']])
        if (any(is.na(gn_sp_tips))) {
            stop(
                'Species tree tip mapping failed for gene-tree species: ',
                paste(gn_spp[is.na(gn_sp_tips)], collapse=', ')
            )
        }
        num_sp = length(gn_spp)
        cat('# species in the gene tree:', num_sp, '\n')
        cat('Species in the gene tree:', paste(gn_spp, collapse=', '), '\n')
        if (length(gn_sp_tips) == 1) {
            num_sp_gntree = get_node_num_by_name(sp_tree, gn_sp_tips)
        } else {
            num_sp_gntree = ape::getMRCA(sp_tree, gn_sp_tips)
        }
        num_sp_gntree = max(1, num_sp_gntree)
        if (num_sp_gntree==get_root_num(sp_tree)) {
            divtime_max = max_age
            divtime_min = max(ape::node.depth.edgelength(sp_tree))
        } else {
            if (length(gn_spp)==1) {
                num_mrca = get_node_num_by_name(sp_tree, gn_sp_tips)
            } else {
                num_mrca = ape::getMRCA(sp_tree, gn_sp_tips)
            }
            num_parent = sp_tree$edge[,1][sp_tree$edge[,2]==num_mrca]
            label_mrca = get_node_name_by_num(phy=sp_tree, node_num=num_mrca)
            label_parent = get_node_name_by_num(phy=sp_tree, node_num=num_parent)
            divtime_max = sp_node_table[sp_node_table$node==label_parent,'age']
            divtime_min = sp_node_table[sp_node_table$node==label_mrca,'age']
            mrca_species = species_parser_get_species_tip_labels(
                species_parser,
                get_descendant_tip_labels(sp_tree, num_mrca),
                strict=TRUE
            )
            parent_species = species_parser_get_species_tip_labels(
                species_parser,
                get_descendant_tip_labels(sp_tree, num_parent),
                strict=TRUE
            )
            cat('Species in the MRCA species tree clade:', paste(mrca_species, collapse=', '), '\n')
            cat('Species in the parent species tree clade:', paste(parent_species, collapse=', '), '\n')
        }
        cat('Divergence time of the parent species tree clade:', divtime_max, '\n')
        cat('Divergence time of the MRCA species tree clade:', divtime_min, '\n')
        is_upper_na = is.na(gn_node_table$upper_age)
        gn_node_table$lower_age[is_upper_na] = divtime_min
        gn_node_table$upper_age[is_upper_na] = divtime_max
    }
    is_missing_age = is.na(gn_node_table[['lower_age']]) | is.na(gn_node_table[['upper_age']])
    if (any(is_missing_age)) {
        unresolved_nodes = unique(gn_node_table[is_missing_age, 'gn_node'])
        stop(
            'Gene node table contains unresolved age bound(s) for node(s): ',
            paste(unresolved_nodes, collapse=', '),
            '. Check lower_sp_node/upper_sp_node annotations in the reconciled metadata.'
        )
    }
    root_num = get_root_num(gn_tree)
    is_root_row = (gn_node_table$gn_node_num==root_num)
    if (sum(is_root_row) != 1) {
        stop('Gene node table root event mapping failed: expected exactly one row for root node.')
    }
    gn_node_table[is_root_row,'event'] = ensure_root_event_tag(gn_node_table[is_root_row,'event'])

    if (dating_backend == 'chronos') {
    constraint_conflicts_before = find_descendant_constraint_conflicts(gn_node_table, gn_tree, root_num)
    if (nrow(constraint_conflicts_before) > 0) {
        cat('Potential chronos failure risk was detected: descendant constraint is identical to or older than an ancestor constraint.\n')
        max_conflict_report = min(60, nrow(constraint_conflicts_before))
        for (i in seq_len(max_conflict_report)) {
            gn_node_num = constraint_conflicts_before$node[i]
            gn_node_name = get_node_name_by_num(gn_tree, gn_node_num)
            cat(
                paste(
                    c(
                        gn_node_name,
                        gn_node_num,
                        constraint_conflicts_before$child_upper[i],
                        constraint_conflicts_before$child_lower[i],
                        constraint_conflicts_before$ancestor_upper[i]
                    ),
                    collapse='/'
                ),
                '\n'
            )
        }
        if (nrow(constraint_conflicts_before) > max_conflict_report) {
            cat('... and ', nrow(constraint_conflicts_before) - max_conflict_report, ' more conflicting nodes.\n', sep='')
        }
    }
    stabilized_constraints = stabilize_descendant_constraints(gn_node_table, gn_tree, root_num)
    gn_node_table = stabilized_constraints$gn_node_table
    if (nrow(stabilized_constraints$adjusted_nodes) > 0) {
        cat('Calibration constraints were stabilized to avoid chronos failure (name/id/upper_before/lower_before/upper_after/lower_after):\n')
        max_stabilized_report = min(80, nrow(stabilized_constraints$adjusted_nodes))
        for (i in seq_len(max_stabilized_report)) {
            node_i = stabilized_constraints$adjusted_nodes$node[i]
            node_name = get_node_name_by_num(gn_tree, node_i)
            cat(
                paste(
                    c(
                        node_name,
                        node_i,
                        stabilized_constraints$adjusted_nodes$upper_before[i],
                        stabilized_constraints$adjusted_nodes$lower_before[i],
                        stabilized_constraints$adjusted_nodes$upper_after[i],
                        stabilized_constraints$adjusted_nodes$lower_after[i]
                    ),
                    collapse='/'
                ),
                '\n'
            )
        }
        if (nrow(stabilized_constraints$adjusted_nodes) > max_stabilized_report) {
            cat('... and ', nrow(stabilized_constraints$adjusted_nodes) - max_stabilized_report, ' more stabilized nodes.\n', sep='')
        }
    }
    constraint_conflicts_after = find_descendant_constraint_conflicts(gn_node_table, gn_tree, root_num)
    if (nrow(constraint_conflicts_after) > 0) {
        max_constraint_age = suppressWarnings(max(gn_node_table[['upper_age']], na.rm=TRUE))
        if ((!is.finite(max_constraint_age)) || is.na(max_constraint_age)) {
            max_constraint_age = 1
        }
        aggressive_min_deltas = unique(sort(c(
            max(1e-6, max_constraint_age * 1e-6),
            max(1e-4, max_constraint_age * 1e-4),
            max(1e-3, max_constraint_age * 1e-3)
        )))
        for (min_delta_i in aggressive_min_deltas) {
            if (nrow(constraint_conflicts_after) == 0) {
                break
            }
            cat(
                'Retrying calibration stabilization with stronger minimum margin: ',
                signif(min_delta_i, digits=4),
                '\n',
                sep=''
            )
            stabilized_retry = stabilize_descendant_constraints(
                gn_node_table=gn_node_table,
                gn_tree=gn_tree,
                root_num=root_num,
                min_delta=min_delta_i
            )
            gn_node_table = stabilized_retry$gn_node_table
            constraint_conflicts_after = find_descendant_constraint_conflicts(gn_node_table, gn_tree, root_num)
        }
    }
    if (nrow(constraint_conflicts_after) > 0) {
        unresolved_nodes = unique(get_node_name_by_num(gn_tree, constraint_conflicts_after$node))
        stop(
            'Calibration constraint stabilization failed for node(s): ',
            format_limited_values(unresolved_nodes, max_items=50)
        )
    }
    } else {
        cat('MCMCTree backend: preserving input calibration bounds without chronos stabilization.\n')
    }
    cat('\n')
    gn_node_table_calibration = gn_node_table[(gn_node_table[,'gn_node_num']>ape::Ntip(gn_tree)),]
    num_constrained_speciation = sum(grepl('^S', gn_node_table_calibration[,'event']))
    num_constrained_duplication = sum(grepl('^D', gn_node_table_calibration[,'event']))
    cat('Number of constrained speciation nodes:', num_constrained_speciation, '\n')
    cat('Number of constrained duplication nodes:', num_constrained_duplication, '\n')
    if ('shared_speciation_group' %in% colnames(gn_node_table)) {
        shared_groups = unique(gn_node_table[['shared_speciation_group']])
        shared_groups = shared_groups[!is.na(shared_groups) & (shared_groups != '')]
        if (length(shared_groups) > 0) {
            cat('Number of shared speciation age groups:', length(shared_groups), '\n')
        }
    }

    # Calibration table
    calibration_table = data.frame(
        node=as.integer(gn_node_table_calibration$gn_node_num),
        age.min=as.numeric(gn_node_table_calibration$lower_age),
        age.max=as.numeric(gn_node_table_calibration$upper_age),
        soft.bounds=NA,
        stringsAsFactors=FALSE
    )

    calibration_table_R = calibration_table[(calibration_table$node==root_num),]
    if (any(grepl('^S', gn_node_table_calibration$event))) {
        S_nodes = gn_node_table_calibration[
            grepl('^S', gn_node_table_calibration$event) &
            (gn_node_table_calibration$gn_node_num != root_num),
            'gn_node_num'
        ]
        calibration_table_S = calibration_table[calibration_table$node %in% S_nodes,]
    } else {
        calibration_table_S = calibration_table[0,]
    }

    calibration_tables = list(
        'RS' = rbind(calibration_table_R, calibration_table_S),
        'S' = calibration_table_S,
        'R' = calibration_table_R
    )

    # dating
    chronos_out = NULL
    current_calibration_table = calibration_tables[['RS']]
    dating_backend_used = dating_backend
    mcmctree_workdir_used = NA_character_
    mcmctree_executable_used = NA_character_
    mcmctree_banner = NA_character_
    mcmctree_mirror_table = NULL
    shared_speciation_age_summary = empty_shared_speciation_age_summary()
    chronos_model_used = NA_character_
    chronos_lambda_used = NA_real_
    chronos_seed_used = NA_integer_
    has_duplication_event = any(grepl('^D', gn_node_table$event))
    has_interval_calibration = has_species_node_interval_input
    has_shared_interval_speciation = has_shared_interval_speciation_input

    if (dating_backend == 'chronos') {
        if (has_shared_interval_speciation) {
            cat(
                'chronos backend note: mirrored speciation nodes share the same interval bounds, ',
                'but chronos cannot force separate internal nodes to have identical estimated ages unless the bounds are exact.\n',
                sep=''
            )
        }
        chronos_control_profiles = make_chronos_control_profiles(
            iter_max=chronos_iter_max,
            eval_max=chronos_eval_max,
            dual_iter_max=chronos_dual_iter_max,
            enable_high_fallback=chronos_high_control_fallback
        )

        if ((!has_duplication_event) && (!has_interval_calibration)) {
            # Gene tree without duplication nodes
            calibrated_node = "allS"
            cat("Constrained nodes:", calibrated_node, '\n')
            cat("All nodes are speciation nodes. Transferring node ages from species tree without age inference by chronos.", '\n')
            gn_spp = unique(
                get_species_names(
                    phy=gn_tree,
                    species_parser=species_parser,
                    species_tree_labels=sp_tree[['tip.label']]
                )
            )
            gn_sp_tips = resolve_species_tree_tips(species_parser, gn_spp, sp_tree[['tip.label']])
            if (any(is.na(gn_sp_tips))) {
                stop(
                    'Species tree tip mapping failed for gene-tree species: ',
                    paste(gn_spp[is.na(gn_sp_tips)], collapse=', ')
                )
            }
            drop_spp = sp_tree$tip.label[! sp_tree$tip.label %in% gn_sp_tips]
            if (length(drop_spp) > 0) {
                chronos_out = ape::drop.tip(phy=sp_tree, tip=drop_spp, trim.internal = TRUE)
            } else {
                chronos_out = sp_tree
            }
            gn_tip_species = species_parser_get_gene_species(
                species_parser=species_parser,
                tip_labels=gn_tree$tip.label,
                species_tree_labels=sp_tree[['tip.label']],
                strict=TRUE
            )
            chronos_tip_species = species_parser_get_species_tip_labels(
                species_parser=species_parser,
                tip_labels=chronos_out$tip.label,
                strict=TRUE
            )
            gn_tip_index = c()
            for (sp in chronos_tip_species) {
                tip_matches = which(gn_tip_species == sp)
                if (length(tip_matches) != 1) {
                    stop(
                        'Gene tree tip mapping failed for species ',
                        sp,
                        '. Expected exactly one matching tip but found ',
                        length(tip_matches),
                        '.'
                    )
                }
                gn_tip_index = c(gn_tip_index, tip_matches)
            }
            chronos_out$tip.label = gn_tree$tip.label[gn_tip_index]
            chronos_out = transfer_node_labels(phy_from=gn_tree, phy_to=chronos_out)
            current_calibration_table = rbind(calibration_table_R, calibration_table_S)
        } else {
            if (!has_duplication_event) {
                cat(
                    'All nodes are speciation nodes, but interval age bounds are present. ',
                    'Running chronos to honor age.min/age.max constraints.\n',
                    sep=''
                )
            }
            validate_tree_edge_lengths(gn_tree, 'gene tree for chronos')
            # Normalize edge lengths to prevent numerical overflow in chronos
            gn_tree = normalize_edge_length_range(gn_tree, max_edge = 1000, min_edge = 1e-8)
            chronos_out = structure('PLACEHOLDER', class='try-error')
            calibrated_node = 'RS'
            current_calibration_table = calibration_tables[['RS']]

            chronos_model_used = chronos_model
            chronos_lambda_used = chronos_lambda
            chronos_seed_used = NA_integer_
            max_restarts_main = 3
            max_restarts_fallback = 2
            soft_attempts = list(list(model=chronos_model, lambda=chronos_lambda, label='requested'))
            if (chronos_model != 'discrete') {
                soft_attempts[[length(soft_attempts)+1]] = list(model='discrete', lambda=chronos_lambda, label='model-discrete')
            }
            if (chronos_model == 'discrete') {
                if (!isTRUE(all.equal(chronos_lambda, 0.1))) {
                    soft_attempts[[length(soft_attempts)+1]] = list(model='discrete', lambda=0.1, label='lambda0.1')
                }
                if (!isTRUE(all.equal(chronos_lambda, 0))) {
                    soft_attempts[[length(soft_attempts)+1]] = list(model='discrete', lambda=0, label='lambda0')
                }
            }
            if (chronos_model != 'relaxed') {
                soft_attempts[[length(soft_attempts)+1]] = list(model='relaxed', lambda=chronos_lambda, label='model-relaxed')
            }
            if (chronos_model != 'correlated') {
                soft_attempts[[length(soft_attempts)+1]] = list(model='correlated', lambda=chronos_lambda, label='model-correlated')
            }
            chronos_time_budget = create_chronos_time_budget(chronos_total_timeout_sec)

            calibration_sequence = 'RS'
            if (allow_constraint_drop) {
                calibration_sequence = c('RS', 'S', 'R')
            }
            for (control_profile_name in names(chronos_control_profiles)) {
                if (!("try-error" %in% class(chronos_out))) {
                    break
                }
                chronos_control = chronos_control_profiles[[control_profile_name]]
                cat(
                    '\nchronos control profile: ',
                    control_profile_name,
                    ' (iter.max=',
                    chronos_control$iter.max,
                    ', eval.max=',
                    chronos_control$eval.max,
                    ', dual.iter.max=',
                    chronos_control$dual.iter.max,
                    ')\n',
                    sep=''
                )
                if (control_profile_name == 'high-fallback') {
                    cat('Fast chronos control profile was exhausted; enabling high-cost fallback.\n')
                }
                seed_cursor = random_seed
                for (cn in calibration_sequence) {
                    if (!("try-error" %in% class(chronos_out))) {
                        break
                    }
                    stage_calibration = calibration_tables[[cn]]
                    if (nrow(stage_calibration) == 0) {
                        cat("\nchronos, calibrated nodes: ", cn, " (skipped; no calibration nodes)\n", sep='')
                        next
                    }
                    if (cn != 'RS') {
                        stage_index = match(cn, calibration_sequence)
                        prev_cn = calibration_sequence[stage_index - 1]
                        cat(
                            "\nchronos constraint-drop stage: ",
                            prev_cn,
                            " retries were exhausted; retrying with ",
                            cn,
                            " constraints.\n",
                            sep=''
                        )
                    }
                    calibrated_node = cn
                    stage_out = run_chronos_retry_pipeline(
                        phy=gn_tree,
                        calibration_table=stage_calibration,
                        root_num=root_num,
                        chronos_control=chronos_control,
                        chronos_lambda=chronos_lambda,
                        chronos_model=chronos_model,
                        soft_attempts=soft_attempts,
                        calibration_label=cn,
                        max_restarts_main=max_restarts_main,
                        max_restarts_fallback=max_restarts_fallback,
                        seed_cursor=seed_cursor,
                        attempt_timeout_sec=chronos_attempt_timeout_sec,
                        time_budget=chronos_time_budget
                    )
                    seed_cursor = stage_out$seed_cursor
                    chronos_out = stage_out$chronos_out
                    current_calibration_table = stage_out$calibration_table
                    if (stage_out$success) {
                        gn_tree = stage_out$phy
                        chronos_model_used = stage_out$used_model
                        chronos_lambda_used = stage_out$used_lambda
                        chronos_seed_used = stage_out$used_seed
                    }
                }
            }

            if (("try-error" %in% class(chronos_out)) && (!allow_constraint_drop)) {
                cat("\nchronos, calibrated nodes: RS (deterministic no-drop fallback)\n")
                deterministic_out = build_dated_tree_without_chronos(
                    phy=gn_tree,
                    calibration_table=current_calibration_table,
                    root_num=root_num
                )
                if (!("try-error" %in% class(deterministic_out))) {
                    chronos_out = deterministic_out
                    chronos_model_used = 'deterministic-fallback'
                    chronos_lambda_used = NA_real_
                    chronos_seed_used = NA_integer_
                } else {
                    cat('deterministic no-drop fallback: failed -> ', as.character(deterministic_out), '\n', sep='')
                }
            }
        }
    } else if (dating_backend == 'mcmctree') {
        calibrated_node = 'RS'
        current_calibration_table = calibration_tables[['RS']]
        if (nrow(current_calibration_table) == 0) {
            stop('No root/speciation calibration nodes are available for --dating_backend=mcmctree.')
        }
        if ('allow_constraint_drop' %in% names(args)) {
            cat('--allow_constraint_drop is ignored when --dating_backend=mcmctree.\n')
        }
        chronos_only_options = c(
            'chronos_attempt_timeout_sec',
            'chronos_total_timeout_sec',
            'chronos_iter_max',
            'chronos_eval_max',
            'chronos_dual_iter_max',
            'chronos_high_control_fallback'
        )
        if (any(chronos_only_options %in% names(args))) {
            cat('chronos options are ignored when --dating_backend=mcmctree.\n')
        }
        cat('Running MCMCTree backend.\n')
        mcmctree_out = run_mcmctree_backend(
            phy=gn_tree,
            gn_node_table=gn_node_table,
            calibration_table=current_calibration_table,
            root_num=root_num,
            seqfile=mcmctree_seqfile,
            bin=mcmctree_bin,
            workdir=mcmctree_workdir,
            usedata=mcmctree_usedata,
            seqtype=mcmctree_seqtype,
            clock=mcmctree_clock,
            model=mcmctree_model,
            burnin=mcmctree_burnin,
            sampfreq=mcmctree_sampfreq,
            nsample=mcmctree_nsample,
            ncatG=mcmctree_ncatG,
            seed=random_seed,
            timeout_sec=mcmctree_timeout_sec
        )
        chronos_out = mcmctree_out$tree
        mcmctree_workdir_used = mcmctree_out$workdir
        mcmctree_mirror_table = mcmctree_out$mirror_table
        shared_speciation_age_summary = mcmctree_out$shared_speciation_ages
        mcmctree_executable_used = mcmctree_out$executable
        mcmctree_banner = mcmctree_out$banner
    }

    if ("try-error" %in% class(chronos_out)) {
        stop('All attempts for divergence time estimation were failed. Exiting.')
    } else {
        cat('Writing output files.\n')
        chronos_out2 = chronos_out
        num_neg = 1
        counter = 1
        if ((dating_backend == 'chronos') && length(args[['pad_short_edge']])) {
            while ((num_neg>0)&(counter<100)) {
                cat(paste0('Branch length padding round ', counter, ' started.\n'))
                chronos_out2 = pad_short_edges(chronos_out2, threshold=args[['pad_short_edge']], external_only=FALSE)
                chronos_out2 = force_ultrametric(chronos_out2, stop_if_larger_change=0.01)
                num_neg = sum(chronos_out2[['edge.length']]<0)
                cat(num_neg, 'negative value(s) were detected in estimated branch length.\n\n')
                counter = counter + 1
            }
        }

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

        if ((dating_backend == 'mcmctree') && (!is.na(mcmctree_workdir_used))) {
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
                source_artifact = file.path(mcmctree_workdir_used, artifact)
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

            input_files = list(
                species_tree=sp_file,
                gene_tree=if (exists('gn_file', inherits=FALSE)) gn_file else NULL,
                notung_parsable=if (exists('parsable_file', inherits=FALSE)) parsable_file else NULL,
                generax_nhx=if (exists('generax_file', inherits=FALSE)) generax_file else NULL,
                species_map_tsv=species_map_tsv,
                species_node_bounds_tsv=species_node_bounds_tsv,
                mcmctree_seqfile=mcmctree_seqfile
            )
            input_files = input_files[!vapply(input_files, is.null, logical(1))]
            input_hashes = lapply(input_files, compute_file_sha256)
            names(input_hashes) = paste0('input_sha256.', names(input_hashes))
            option_entries = lapply(args, function(value) {
                if (is.infinite(value)) 'inf' else value
            })
            names(option_entries) = paste0('option.', names(option_entries))
            manifest_entries = c(
                list(
                    radte_version=radte_version,
                    run_started_at=run_started_at,
                    run_completed_at=format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'),
                    command=paste(c('radte', command_args), collapse=' '),
                    r_version=R.version.string,
                    ape_version=as.character(utils::packageVersion('ape')),
                    treeio_version=if (requireNamespace('treeio', quietly=TRUE)) as.character(utils::packageVersion('treeio')) else NA_character_,
                    dating_backend=dating_backend_used,
                    calibrated_nodes=calibrated_node,
                    seed_requested=random_seed,
                    effective_seed=if (!is.na(chronos_seed_used)) chronos_seed_used else random_seed,
                    chronos_model_used=chronos_model_used,
                    chronos_lambda_used=chronos_lambda_used,
                    chronos_seed_used=chronos_seed_used,
                    mcmctree_executable=mcmctree_executable_used,
                    mcmctree_banner=mcmctree_banner,
                    mcmctree_workdir=mcmctree_workdir_used,
                    output_directory=outdir,
                    output_prefix=prefix
                ),
                option_entries,
                input_hashes
            )
            manifest_file = output_path('_run_manifest.tsv')
            write_run_manifest(manifest_entries, manifest_file)

            cat('dating backend used:', dating_backend_used, '\n')
            if (!is.na(chronos_model_used)) {
                cat('chronos model used:', chronos_model_used, '\n')
            }
            if (!is.na(chronos_lambda_used)) {
                cat('chronos lambda used:', chronos_lambda_used, '\n')
            }
            if (!is.na(chronos_seed_used)) {
                cat('chronos seed used:', chronos_seed_used, '\n')
            }
            if (!is.na(mcmctree_workdir_used)) {
                cat('MCMCTree workdir:', mcmctree_workdir_used, '\n')
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
                dating_backend=dating_backend_used,
                calibrated_nodes=calibrated_node
            )))
    }
}
