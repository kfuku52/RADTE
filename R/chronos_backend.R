validate_tree_edge_lengths = function(tree, tree_name) {
    edge_lengths = tree[['edge.length']]
    if (is.null(edge_lengths) || (length(edge_lengths)==0)) {
        stop('Input ', tree_name, ' does not contain branch length(s).')
    }
    if (length(edge_lengths) != nrow(tree[['edge']])) {
        stop(
            'Input ', tree_name,
            ' has inconsistent branch length and edge table lengths (',
            length(edge_lengths), ' vs ', nrow(tree[['edge']]), ').'
        )
    }
    if (any(is.na(edge_lengths) | (!is.finite(edge_lengths)))) {
        stop('Input ', tree_name, ' contains missing or non-finite branch length(s).')
    }
    if (any(edge_lengths < 0)) {
        stop('Input ', tree_name, ' contains negative branch length(s).')
    }
}

# Small indirection keeps external numerical work mockable in unit tests while
# avoiding search-path attachment in the package and generated CLI.
radte_chronos = function(...) {
    ape::chronos(...)
}

# copied from rkftools https://github.com/kfuku52/rkftools
transfer_node_labels = function(phy_from, phy_to) {
    if (is.null(phy_to$node.label) || is.null(phy_from$node.label)) {
        return(phy_to)
    }
    from_index = match_internal_nodes_by_clade(phy_from, phy_to)
    matched = !is.na(from_index)
    if (any(matched)) {
        phy_to$node.label[matched] = phy_from$node.label[from_index[matched]]
    }
    return(phy_to)
}

check_gn_node_name_uniqueness = function(gn_node_table, gn_tree)
for (gn_node_name in gn_node_table[,'gn_node']) {
    n = get_node_num_by_name(gn_tree, gn_node_name)
    if (length(n) == 0) {
        stop(paste('Input gene tree does not contain node name:', gn_node_name))
    }
    if (length(n) > 1) {
        stop(paste('Input gene tree contains multiple nodes with the identical name:', gn_node_name))
    }
}

pad_branch_length = function(tree, pad_size=1e-6) {
  is_bl_zero = (tree[['edge.length']]<pad_size)
    if (any(is_bl_zero)) {
        txt = paste0(sum(is_bl_zero), ' out of ', length(is_bl_zero))
        txt = paste0(txt, ' branches have small length. Padding with ', pad_size, '.\n')
        cat(txt)
        tree[['edge.length']][is_bl_zero] = pad_size
    }
    return(tree)
}

adjust_branch_length_order = function(tree, min_bl=1e-6) {
    is_bl_zero = (tree[['edge.length']]<=0)
    if (any(is_bl_zero)) {
        stop('The tree contains branch(es) with zero or negative length.')
    }
    min_tree_bl = min(tree[['edge.length']])
    for (i in 1:20) {
        if (min_tree_bl<min_bl) {
            tree[['edge.length']] = tree[['edge.length']] * 10
            min_tree_bl =  min(tree[['edge.length']])
        } else {
            break
        }
    }
    return(tree)
}

normalize_edge_length_range = function(tree, max_edge = 1000, min_edge = 1e-8) {
    # Scale down edge lengths if max is too large to prevent numerical overflow in chronos.
    # Chronos fails with NaN errors when edge lengths exceed ~6000-7000.
    # Using max_edge=1000 provides a safe margin.
    edges <- tree$edge.length
    if (max(edges) > max_edge) {
        scale_factor <- max_edge / max(edges)
        tree$edge.length <- edges * scale_factor
        cat('Edge lengths scaled by factor:', scale_factor, 'to prevent numerical overflow.\n')
    }
    # Ensure minimum edge length is not too small (avoids underflow after scaling)
    tree$edge.length[tree$edge.length < min_edge] <- min_edge
    return(tree)
}

find_descendant_constraint_conflicts = function(gn_node_table, gn_tree, root_num) {
    conflicts = data.frame(
        node=integer(),
        child_lower=numeric(),
        child_upper=numeric(),
        ancestor_lower=numeric(),
        ancestor_upper=numeric(),
        stringsAsFactors=FALSE
    )
    if (nrow(gn_node_table)==0) {
        return(conflicts)
    }
    node_nums = as.integer(gn_node_table[['gn_node_num']])
    lower_context = get_ancestor_constraint_minima(
        gn_tree, node_nums, gn_node_table[['lower_age']]
    )[['minimum']]
    upper_context = get_ancestor_constraint_minima(
        gn_tree, node_nums, gn_node_table[['upper_age']]
    )[['minimum']]
    ancestor_lower = lower_context[node_nums]
    ancestor_upper = upper_context[node_nums]
    is_conflict = (node_nums != root_num) & is.finite(ancestor_lower) &
        is.finite(ancestor_upper) &
        (gn_node_table[['lower_age']] >= ancestor_upper)
    if (!any(is_conflict)) {
        return(conflicts)
    }
    data.frame(
        node=node_nums[is_conflict],
        child_lower=gn_node_table[['lower_age']][is_conflict],
        child_upper=gn_node_table[['upper_age']][is_conflict],
        ancestor_lower=ancestor_lower[is_conflict],
        ancestor_upper=ancestor_upper[is_conflict],
        stringsAsFactors=FALSE
    )
}

stabilize_descendant_constraints = function(gn_node_table, gn_tree, root_num, min_delta=NULL) {
    if (nrow(gn_node_table)) {
        calibration_feasible_ranges(gn_tree, data.frame(
            node=gn_node_table$gn_node_num,
            age.min=gn_node_table$lower_age,
            age.max=gn_node_table$upper_age
        ), min_edge=min_delta)
    }
    list(gn_node_table=gn_node_table, adjusted_nodes=data.frame(), min_delta=0)
}

detect_chronos_failure_risks = function(tree, calibration_table, root_num, edge_ratio_threshold=1e8, max_edge_threshold=5000, min_span_threshold=NULL) {
    edge_lengths = tree[['edge.length']]
    max_edge = max(edge_lengths)
    positive_edges = edge_lengths[edge_lengths>0]
    if (length(positive_edges)==0) {
        min_positive_edge = NA
        edge_ratio = Inf
    } else {
        min_positive_edge = min(positive_edges)
        edge_ratio = max_edge / min_positive_edge
    }

    if (is.null(min_span_threshold)) {
        if (nrow(calibration_table)==0) {
            max_cal_age = 1
        } else {
            max_cal_age = suppressWarnings(max(calibration_table[['age.max']], na.rm=TRUE))
            if ((!is.finite(max_cal_age)) || is.na(max_cal_age)) {
                max_cal_age = 1
            }
        }
        min_span_threshold = max(1e-8, max_cal_age * 1e-8)
    }

    spans = calibration_table$age.max - calibration_table$age.min
    non_root = calibration_table$node != root_num
    tight_nodes = calibration_table$node[non_root & (spans < min_span_threshold)]
    nonpositive_nodes = calibration_table$node[non_root & (spans <= 0)]
    risk_flags = c(
        extreme_edge_ratio = is.finite(edge_ratio) && (edge_ratio > edge_ratio_threshold),
        large_edge_values = (max_edge > max_edge_threshold),
        tight_nonroot_constraints = (length(tight_nodes) > 0),
        nonpositive_nonroot_constraints = (length(nonpositive_nodes) > 0)
    )

    return(
        list(
            max_edge=max_edge,
            min_positive_edge=min_positive_edge,
            edge_ratio=edge_ratio,
            tight_nodes=tight_nodes,
            nonpositive_nodes=nonpositive_nodes,
            min_span_threshold=min_span_threshold,
            risk_flags=risk_flags
        )
    )
}

validate_chronos_output = function(chronos_out) {
    if ("try-error" %in% class(chronos_out)) {
        return(list(valid=FALSE, reason=as.character(chronos_out)))
    }
    if (!inherits(chronos_out, 'chronos')) {
        return(list(valid=FALSE, reason='Output does not inherit class "chronos".'))
    }
    edge_lengths = chronos_out[['edge.length']]
    if (is.null(edge_lengths) || (length(edge_lengths)==0)) {
        return(list(valid=FALSE, reason='No edge lengths were returned.'))
    }
    if (any(is.na(edge_lengths) | (!is.finite(edge_lengths)))) {
        return(list(valid=FALSE, reason='Edge lengths contain non-finite values.'))
    }
    if (any(edge_lengths<=0)) {
        return(list(valid=FALSE, reason='Edge lengths contain zero or negative values.'))
    }
    ploglik = attr(chronos_out, 'ploglik')
    if (!is.null(ploglik)) {
        if ((!is.finite(ploglik)) || (ploglik <= -1e50)) {
            return(list(valid=FALSE, reason=paste('Penalized log-likelihood is invalid:', ploglik)))
        }
    }
    phiic = attr(chronos_out, 'PHIIC')
    if (!is.null(phiic)) {
        if (is.list(phiic)) {
            if ((!is.null(phiic[['logLik']])) && (!is.finite(phiic[['logLik']]))) {
                return(list(valid=FALSE, reason='PHIIC logLik is non-finite.'))
            }
            if ((!is.null(phiic[['PHIIC']])) && (!is.finite(phiic[['PHIIC']]))) {
                return(list(valid=FALSE, reason='PHIIC score is non-finite.'))
            }
        } else if (!is.finite(phiic)) {
            return(list(valid=FALSE, reason='PHIIC is non-finite.'))
        }
    }
    convergence = attr(chronos_out, 'convergence')
    if (!is.null(convergence) && identical(convergence, FALSE)) {
        return(list(valid=FALSE, reason='chronos reported non-convergence.'))
    }
    rates = attr(chronos_out, 'rates')
    if (!is.null(rates) && any(!is.finite(rates))) {
        return(list(valid=FALSE, reason='Estimated rates contain non-finite values.'))
    }
    frequencies = attr(chronos_out, 'frequencies')
    if (!is.null(frequencies) && any(!is.finite(frequencies))) {
        return(list(valid=FALSE, reason='Estimated frequencies contain non-finite values.'))
    }
    return(list(valid=TRUE, reason='OK'))
}

make_chronos_control_profiles = function(
    iter_max=250L,
    eval_max=250L,
    dual_iter_max=20L,
    enable_high_fallback=TRUE,
    fallback_iter_max=100000L,
    fallback_eval_max=100000L,
    fallback_dual_iter_max=200L
) {
    fast_control = ape::chronos.control()
    fast_control$iter.max = as.integer(iter_max)
    fast_control$eval.max = as.integer(eval_max)
    fast_control$dual.iter.max = as.integer(dual_iter_max)
    profiles = list(fast=fast_control)
    if (isTRUE(enable_high_fallback)) {
        fallback_control = ape::chronos.control()
        fallback_control$iter.max = as.integer(max(iter_max, fallback_iter_max))
        fallback_control$eval.max = as.integer(max(eval_max, fallback_eval_max))
        fallback_control$dual.iter.max = as.integer(max(dual_iter_max, fallback_dual_iter_max))
        is_distinct = !identical(fast_control$iter.max, fallback_control$iter.max) ||
            !identical(fast_control$eval.max, fallback_control$eval.max) ||
            !identical(fast_control$dual.iter.max, fallback_control$dual.iter.max)
        if (is_distinct) {
            profiles[['high-fallback']] = fallback_control
        }
    }
    return(profiles)
}

create_chronos_time_budget = function(total_timeout_sec=Inf) {
    budget = new.env(parent=emptyenv())
    budget$enabled = is.finite(total_timeout_sec) && (total_timeout_sec > 0)
    budget$total_timeout_sec = total_timeout_sec
    budget$start_time = Sys.time()
    return(budget)
}

get_chronos_budget_remaining = function(time_budget=NULL) {
    if (is.null(time_budget) || (!isTRUE(time_budget$enabled))) {
        return(Inf)
    }
    elapsed_sec = as.numeric(difftime(Sys.time(), time_budget$start_time, units='secs'))
    return(max(0, time_budget$total_timeout_sec - elapsed_sec))
}

run_chronos_once = function(phy, calibration, chronos_control, chronos_lambda, chronos_model, timeout_sec=Inf) {
    if ((!is.finite(timeout_sec)) || (timeout_sec <= 0)) {
        return(
            try(
                radte_chronos(
                    phy=phy,
                    lambda=chronos_lambda,
                    model=chronos_model,
                    calibration=calibration,
                    control=chronos_control,
                    quiet=TRUE
                ),
                silent=TRUE
            )
        )
    }

    # setTimeLimit does not always interrupt long C-level loops inside chronos.
    # Run in a child process on Unix to enforce hard wall-time termination.
    if (.Platform$OS.type == 'unix') {
        child = parallel::mcparallel(
            expr=try(
                radte_chronos(
                    phy=phy,
                    lambda=chronos_lambda,
                    model=chronos_model,
                    calibration=calibration,
                    control=chronos_control,
                    quiet=TRUE
                ),
                silent=TRUE
            ),
            silent=TRUE,
            mc.set.seed=FALSE
        )
        poll_interval = min(0.1, timeout_sec / 10)
        if ((!is.finite(poll_interval)) || (poll_interval <= 0)) {
            poll_interval = 0.05
        }
        collected = NULL
        start_time = Sys.time()
        repeat {
            collected = parallel::mccollect(child, wait=FALSE, timeout=0)
            if (!is.null(collected)) {
                break
            }
            elapsed_sec = as.numeric(difftime(Sys.time(), start_time, units='secs'))
            if (elapsed_sec >= timeout_sec) {
                break
            }
            Sys.sleep(poll_interval)
        }
        if (is.null(collected)) {
            if ((!is.null(child$pid)) && (length(child$pid)==1) && is.finite(child$pid)) {
                try(tools::pskill(child$pid), silent=TRUE)
            }
            # Try to reap the child process if available.
            try(parallel::mccollect(child, wait=FALSE, timeout=0), silent=TRUE)
            return(
                structure(
                    paste('Chronos attempt timed out after', signif(timeout_sec, digits=4), 'seconds.'),
                    class='try-error'
                )
            )
        }
        if (length(collected)==0) {
            return(structure('Chronos attempt returned no result.', class='try-error'))
        }
        return(collected[[1]])
    }

    # Fallback for non-Unix platforms.
    setTimeLimit(cpu=Inf, elapsed=timeout_sec, transient=TRUE)
    out = try(
        radte_chronos(
            phy=phy,
            lambda=chronos_lambda,
            model=chronos_model,
            calibration=calibration,
            control=chronos_control,
            quiet=TRUE
        ),
        silent=TRUE
    )
    setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE)
    return(out)
}

run_chronos_with_restarts = function(
    phy,
    calibration,
    chronos_control,
    chronos_lambda,
    chronos_model,
    context_label='chronos',
    max_restarts=3,
    seed_base=1,
    attempt_timeout_sec=Inf,
    time_budget=NULL
) {
    if (max_restarts < 1) {
        max_restarts = 1
    }
    if ((is.finite(attempt_timeout_sec)) && (attempt_timeout_sec <= 0)) {
        attempt_timeout_sec = Inf
    }
    had_seed = exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)
    previous_seed = NULL
    if (had_seed) {
        previous_seed = get('.Random.seed', envir=.GlobalEnv, inherits=FALSE)
    }
    on.exit({
        if (had_seed) {
            assign('.Random.seed', previous_seed, envir=.GlobalEnv)
        } else {
            if (exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)) {
                rm('.Random.seed', envir=.GlobalEnv)
            }
        }
    }, add=TRUE)

    out = structure('PLACEHOLDER', class='try-error')
    used_seed = NA_integer_
    for (attempt_i in seq_len(max_restarts)) {
        budget_remaining_sec = get_chronos_budget_remaining(time_budget)
        if (is.finite(budget_remaining_sec) && (budget_remaining_sec <= 0)) {
            out = structure(
                'Chronos total time budget was exhausted before attempting another retry.',
                class='try-error'
            )
            break
        }
        attempt_timeout_i = attempt_timeout_sec
        if (is.finite(budget_remaining_sec)) {
            attempt_timeout_i = min(attempt_timeout_i, budget_remaining_sec)
        }
        if (is.finite(attempt_timeout_i) && (attempt_timeout_i <= 0)) {
            out = structure(
                'Chronos attempt timeout became non-positive because total time budget was exhausted.',
                class='try-error'
            )
            break
        }
        current_seed = radte_attempt_seed(seed_base, attempt_i - 1L)
        suppressWarnings(set.seed(current_seed))
        used_seed = current_seed
        timeout_label = if (is.finite(attempt_timeout_i)) {
            format(signif(attempt_timeout_i, digits=4), scientific=FALSE)
        } else {
            'inf'
        }
        cat(
            context_label,
            ': attempt',
            attempt_i,
            'of',
            max_restarts,
            '(model=',
            chronos_model,
            ', lambda=',
            chronos_lambda,
            ', seed=',
            current_seed,
            ', timeout_sec=',
            timeout_label,
            ')\n',
            sep=''
        )
        out = run_chronos_once(
            phy=phy,
            calibration=calibration,
            chronos_control=chronos_control,
            chronos_lambda=chronos_lambda,
            chronos_model=chronos_model,
            timeout_sec=attempt_timeout_i
        )
        if (("try-error" %in% class(out)) && grepl('elapsed time limit', as.character(out), ignore.case=TRUE)) {
            out = structure(
                paste('Chronos attempt timed out after', signif(attempt_timeout_i, digits=4), 'seconds.'),
                class='try-error'
            )
        }
        if (!("try-error" %in% class(out))) {
            validation = validate_chronos_output(out)
            if (validation$valid) validation = dated_tree_validation(out, phy, calibration)
            if (isTRUE(validation$valid)) {
                return(
                    list(
                        chronos_out=out,
                        success=TRUE,
                        used_model=chronos_model,
                        used_lambda=chronos_lambda,
                        used_seed=used_seed
                    )
                )
            }
            out = structure(
                paste('Invalid chronos output:', validation$reason),
                class='try-error'
            )
        }
        cat(context_label, ': failed -> ', as.character(out), '\n', sep='')
    }
    return(
        list(
            chronos_out=out,
            success=FALSE,
            used_model=chronos_model,
            used_lambda=chronos_lambda,
            used_seed=used_seed
        )
    )
}

run_chronos_retry_pipeline = function(
    phy, calibration_table, root_num, chronos_control, chronos_lambda, chronos_model,
    alternate_strategies, calibration_label='RS', max_restarts_main=3,
    max_restarts_fallback=2, seed_cursor=1, attempt_timeout_sec=Inf, time_budget=NULL
) {
    failure = try(calibration_feasible_ranges(phy, calibration_table), silent=TRUE)
    out = list(chronos_out=structure('No calibration nodes were available.', class='try-error'),
               success=FALSE, used_model=chronos_model, used_lambda=chronos_lambda,
               used_seed=NA_integer_, seed_cursor=seed_cursor, phy=phy,
               calibration_table=calibration_table, skipped=nrow(calibration_table) == 0L)
    if (out$skipped) return(out)
    if (inherits(failure, 'try-error')) {
        out$chronos_out = failure
        return(out)
    }
    strategies = list(list(phy=phy, model=chronos_model, lambda=chronos_lambda,
                           label='', restarts=max_restarts_main))
    risk = detect_chronos_failure_risks(phy, calibration_table, root_num)
    if (any(risk$risk_flags[c('extreme_edge_ratio', 'large_edge_values')])) {
        phy = normalize_edge_length_range(phy, max_edge=100, min_edge=1e-6)
        strategies[[length(strategies) + 1L]] = list(phy=phy, model=chronos_model,
            lambda=chronos_lambda, label=' rescaled', restarts=max_restarts_main)
    }
    # Only model/lambda vary; all retries retain the original hard bounds.
    for (strategy in alternate_strategies) {
        if (strategy$model == chronos_model && strategy$lambda == chronos_lambda) next
        strategies[[length(strategies) + 1L]] = list(phy=phy, model=strategy$model,
            lambda=strategy$lambda, label=paste0(' alternate ', strategy$label),
            restarts=max_restarts_fallback)
    }
    for (strategy in strategies) {
        if (get_chronos_budget_remaining(time_budget) <= 0) break
        result = run_chronos_with_restarts(
            phy=strategy$phy, calibration=calibration_table, chronos_control=chronos_control,
            chronos_lambda=strategy$lambda, chronos_model=strategy$model,
            context_label=paste0('chronos ', calibration_label, strategy$label),
            max_restarts=strategy$restarts, seed_base=seed_cursor,
            attempt_timeout_sec=attempt_timeout_sec, time_budget=time_budget)
        seed_cursor = radte_attempt_seed(seed_cursor, strategy$restarts)
        out$chronos_out = result$chronos_out
        out$seed_cursor = seed_cursor
        if (result$success) {
            out$success = TRUE
            out$phy = strategy$phy
            out$used_model = result$used_model
            out$used_lambda = result$used_lambda
            out$used_seed = result$used_seed
            break
        }
    }
    out
}

build_dated_tree_without_chronos = function(phy, calibration_table, root_num, min_edge=NULL, require_root=TRUE) {
    tryCatch({
        if (!nrow(calibration_table)) stop('No calibration constraints were provided.')
        root_row = match(root_num, calibration_table$node)
        if (is.na(root_row) && require_root) stop('Root calibration is missing for deterministic fallback.')
        feasible = calibration_feasible_ranges(phy, calibration_table, min_edge)
        depths = ape::node.depth.edgelength(phy)
        if (any(!is.finite(depths))) stop('Could not compute finite node depths.')
        root_target = if (!require_root) {
            min(feasible$upper[[root_num]], max(max(depths), feasible$lower[[root_num]]))
        } else calibration_table$age.max[[root_row]]
        baseline = max(depths) - depths
        if (max(baseline) > 0) baseline = baseline * root_target / max(baseline)
        ages = baseline
        ages[[root_num]] = root_target
        for (i in seq_len(nrow(feasible$preorder))) {
            parent = feasible$preorder[i, 1]
            child = feasible$preorder[i, 2]
            upper = min(feasible$upper[[child]], ages[[parent]] - feasible$min_edge)
            ages[[child]] = max(feasible$lower[[child]], min(baseline[[child]], upper))
        }
        ages[seq_along(phy$tip.label)] = 0
        out = phy
        out$edge.length = ages[phy$edge[, 1]] - ages[phy$edge[, 2]]
        validate_dated_tree(out, phy, calibration_table)
        out
    }, error=function(e) structure(conditionMessage(e), class='try-error'))
}
