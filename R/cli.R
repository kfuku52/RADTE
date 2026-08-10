get_parsed_args = function(args, print=TRUE) {
    parsed = list()
    if (length(args)==0) {
        return(parsed)
    }
    split = strsplit(sub("^--", "", args), "=", fixed=TRUE)
    for (i in seq_along(split)) {
        if (!startsWith(args[[i]], "--")) {
            stop('Argument is not in --key=value format: ', args[[i]])
        }
        if (length(split[[i]]) < 2) {
            stop('Argument is not in --key=value format: ', args[[i]])
        }
        param = gsub("-", "_", split[[i]][1], fixed=TRUE)
        if (nchar(param)==0) {
            stop('Argument key is empty in --key=value format: ', args[[i]])
        }
        if (param %in% names(parsed)) {
            stop('Argument key is duplicated: --', param)
        }
        parsed[[param]] = paste(split[[i]][2:length(split[[i]])], collapse='=')
    }
    if (print) {
        for (name in names(parsed)) {
            cat(name, '=', parsed[[name]], '\n')
        }
        cat('\n')
    }
    return(parsed)
}

radte_option = function(type, default=NULL, help, choices=NULL, min_value=NULL, backend='all') {
    list(
        type=type,
        default=default,
        help=help,
        choices=choices,
        min_value=min_value,
        backend=backend
    )
}

get_radte_option_schema = function() {
    list(
        allow_constraint_drop=radte_option('boolean', TRUE, 'Allow S/R constraint-drop fallback stages.', backend='chronos'),
        chronos_attempt_timeout_sec=radte_option('timeout', 60, 'Wall-time limit per chronos attempt; 0/inf/none/off disables it.', backend='chronos'),
        chronos_dual_iter_max=radte_option('integer', 20L, 'Initial chronos dual.iter.max.', min_value=0L, backend='chronos'),
        chronos_eval_max=radte_option('integer', 250L, 'Initial chronos eval.max.', min_value=1L, backend='chronos'),
        chronos_high_control_fallback=radte_option('boolean', TRUE, 'Enable the high-cost chronos control fallback.', backend='chronos'),
        chronos_iter_max=radte_option('integer', 250L, 'Initial chronos iter.max.', min_value=1L, backend='chronos'),
        chronos_lambda=radte_option('number', NULL, 'Chronos smoothing parameter (required for chronos).', min_value=0, backend='chronos'),
        chronos_model=radte_option('choice', NULL, 'Chronos model (required for chronos).', choices=c('discrete', 'relaxed', 'correlated'), backend='chronos'),
        chronos_total_timeout_sec=radte_option('timeout', 300, 'Total wall-time budget for chronos retries; 0/inf/none/off disables it.', backend='chronos'),
        dating_backend=radte_option('choice', 'chronos', 'Dating engine.', choices=c('chronos', 'mcmctree')),
        gene_tree=radte_option('path', NULL, 'Rooted Newick gene tree (required in NOTUNG mode).'),
        generax_nhx=radte_option('path', NULL, 'GeneRax NHX reconciliation input.'),
        max_age=radte_option('number', NULL, 'Upper age for duplication nodes above the species-tree root.', min_value=0),
        mcmctree_bin=radte_option('path', 'mcmctree', 'MCMCTree executable name or path.', backend='mcmctree'),
        mcmctree_burnin=radte_option('integer', 2000L, 'MCMCTree burn-in samples.', min_value=1L, backend='mcmctree'),
        mcmctree_clock=radte_option('integer', 2L, 'MCMCTree clock model.', min_value=1L, backend='mcmctree'),
        mcmctree_model=radte_option('integer', 0L, 'MCMCTree substitution model.', min_value=0L, backend='mcmctree'),
        mcmctree_ncatG=radte_option('integer', 5L, 'MCMCTree number of gamma categories.', min_value=1L, backend='mcmctree'),
        mcmctree_nsample=radte_option('integer', 20000L, 'MCMCTree number of posterior samples.', min_value=1L, backend='mcmctree'),
        mcmctree_sampfreq=radte_option('integer', 10L, 'MCMCTree sampling frequency.', min_value=1L, backend='mcmctree'),
        mcmctree_seqfile=radte_option('path', NULL, 'MCMCTree alignment file (required for mcmctree).', backend='mcmctree'),
        mcmctree_seqtype=radte_option('integer', 0L, 'MCMCTree sequence type.', min_value=0L, backend='mcmctree'),
        mcmctree_timeout_sec=radte_option('timeout', Inf, 'Wall-time limit for the external MCMCTree process.', backend='mcmctree'),
        mcmctree_usedata=radte_option('integer', 1L, 'MCMCTree usedata value; RADTE currently supports only 1.', min_value=1L, backend='mcmctree'),
        mcmctree_workdir=radte_option('path', NULL, 'Explicit MCMCTree staging directory; default is a unique directory under outdir.', backend='mcmctree'),
        notung_parsable=radte_option('path', NULL, 'NOTUNG parsable reconciliation input.'),
        outdir=radte_option('path', '.', 'Directory for RADTE output artifacts.'),
        pad_short_edge=radte_option('number', NULL, 'Minimum dated branch length.', min_value=0),
        prefix=radte_option('string', 'radte', 'Output filename prefix.'),
        seed=radte_option('integer', 1L, 'Base random seed for chronos retries and MCMCTree.', min_value=1L),
        species_map_tsv=radte_option('path', NULL, 'Species mapping TSV for --species_parser=map.'),
        species_node_bounds_tsv=radte_option('path', NULL, 'Species-tree node age bounds TSV.'),
        species_parser=radte_option('choice', 'legacy', 'Species-label parser.', choices=c('legacy', 'taxonomic', 'regex', 'map')),
        species_regex=radte_option('string', NULL, 'Species extraction regex for --species_parser=regex.'),
        species_tree=radte_option('path', NULL, 'Rooted, dated, fully bifurcating species tree.')
    )
}

get_radte_allowed_args = function() {
    names(get_radte_option_schema())
}

validate_parsed_args = function(args, allowed_args=get_radte_allowed_args()) {
    unknown_args = setdiff(names(args), allowed_args)
    if (length(unknown_args) > 0) {
        stop(
            'Unknown argument(s): ',
            paste(paste0('--', gsub('_', '-', unknown_args, fixed=TRUE)), collapse=', '),
            '. Allowed arguments are: ',
            paste(paste0('--', gsub('_', '-', sort(allowed_args), fixed=TRUE)), collapse=', '),
            '.'
        )
    }
    return(invisible(args))
}

coerce_radte_args = function(args, schema=get_radte_option_schema()) {
    for (name in names(args)) {
        spec = schema[[name]]
        arg_name = paste0('--', name)
        value = args[[name]]
        if (spec$type == 'boolean') {
            value = parse_bool_arg(value, arg_name)
        } else if (spec$type == 'integer') {
            value = parse_integer_arg(value, arg_name, min_value=spec$min_value)
        } else if (spec$type == 'timeout') {
            value = parse_timeout_arg(value, arg_name)
        } else if (spec$type == 'choice') {
            # Preserve the dedicated chronos typo suggestion in radte_main.
            if (name != 'chronos_model') {
                value = parse_choice_arg(value, arg_name, spec$choices)
            }
        } else if (spec$type == 'number') {
            value_num = suppressWarnings(as.numeric(value))
            if (length(value_num) != 1 || is.na(value_num) || !is.finite(value_num)) {
                value = value
            } else {
                value = value_num
            }
        } else if (spec$type %in% c('path', 'string')) {
            value = as.character(value)
        }
        args[[name]] = value
    }
    args
}

parse_radte_cli_args = function(command_args, print=TRUE) {
    parsed = get_parsed_args(command_args, print=FALSE)
    validate_parsed_args(parsed)
    parsed = coerce_radte_args(parsed)
    if (print) {
        for (name in names(parsed)) {
            value = if (is.infinite(parsed[[name]])) 'inf' else parsed[[name]]
            cat(name, '=', value, '\n')
        }
        cat('\n')
    }
    parsed
}

format_radte_default = function(value) {
    if (is.null(value)) {
        return('required/conditional')
    }
    if (is.logical(value)) {
        return(tolower(as.character(value)))
    }
    if (is.infinite(value)) {
        return('inf')
    }
    as.character(value)
}

render_radte_help = function() {
    schema = get_radte_option_schema()
    lines = c(
        paste0('RADTE ', radte_version),
        '',
        'Usage:',
        '  radte --species_tree=FILE (--notung_parsable=FILE --gene_tree=FILE | --generax_nhx=FILE) [options]',
        '',
        'Meta options:',
        '  --help, -h       Show this help and exit.',
        '  --version, -V    Show the RADTE version and exit.',
        '',
        'Options use --key=value syntax:'
    )
    for (name in names(schema)) {
        spec = schema[[name]]
        choice_text = if (is.null(spec$choices)) '' else paste0(' Choices: ', paste(spec$choices, collapse=', '), '.')
        lines = c(
            lines,
            paste0('  --', name, '=VALUE'),
            paste0('      ', spec$help, ' [default: ', format_radte_default(spec$default), '; backend: ', spec$backend, '].', choice_text)
        )
    }
    paste(lines, collapse='\n')
}

parse_bool_arg = function(value, arg_name) {
    if ((is.logical(value)) && (length(value)==1) && (!is.na(value))) {
        return(value)
    }
    if ((is.numeric(value)) && (length(value)==1) && (!is.na(value)) && (value %in% c(0, 1))) {
        return(as.logical(value))
    }
    if ((is.character(value)) && (length(value)==1) && (!is.na(value))) {
        value_lower = tolower(trimws(value))
        if (value_lower %in% c('true', 't', 'yes', 'y', '1')) {
            return(TRUE)
        }
        if (value_lower %in% c('false', 'f', 'no', 'n', '0')) {
            return(FALSE)
        }
    }
    stop(arg_name, ' should be boolean (true/false, yes/no, or 1/0).')
}

parse_timeout_arg = function(value, arg_name) {
    if ((is.character(value)) && (length(value)==1) && (!is.na(value))) {
        value_lower = tolower(trimws(value))
        if (value_lower %in% c('inf', 'infinity', 'none', 'off', 'disable', 'disabled')) {
            return(Inf)
        }
    }
    value_num = suppressWarnings(as.numeric(value))
    if ((!is.na(value_num)) && is.infinite(value_num)) {
        if (value_num > 0) {
            return(Inf)
        }
        stop(arg_name, ' should be a non-negative number or one of: inf, none, off.')
    }
    if (is.na(value_num) || (value_num < 0)) {
        stop(arg_name, ' should be a non-negative number or one of: inf, none, off.')
    }
    if (value_num == 0) {
        return(Inf)
    }
    return(value_num)
}

parse_choice_arg = function(value, arg_name, supported_values) {
    if ((!is.character(value)) || (length(value) != 1) || is.na(value) || (nchar(trimws(value)) == 0)) {
        stop(arg_name, ' should be one of: ', paste(supported_values, collapse=', '))
    }
    value_norm = tolower(trimws(as.character(value)))
    if (!value_norm %in% supported_values) {
        stop(arg_name, ' should be one of: ', paste(supported_values, collapse=', '), '. Received: ', value)
    }
    return(value_norm)
}

parse_integer_arg = function(value, arg_name, min_value=NULL) {
    value_num = suppressWarnings(as.numeric(value))
    if (is.na(value_num) || (!is.finite(value_num)) || (value_num != round(value_num))) {
        stop(arg_name, ' should be an integer.')
    }
    value_int = as.integer(round(value_num))
    if (!is.null(min_value) && (value_int < min_value)) {
        stop(arg_name, ' should be >= ', min_value, '.')
    }
    return(value_int)
}

format_limited_values = function(values, max_items=50) {
    values = as.character(values)
    if (length(values)==0) {
        return('')
    }
    if (length(values) <= max_items) {
        return(paste(values, collapse=', '))
    }
    shown = paste(values[seq_len(max_items)], collapse=', ')
    return(paste0(shown, ', ... (', length(values)-max_items, ' more)'))
}
