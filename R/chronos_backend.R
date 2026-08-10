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
        (gn_node_table[['lower_age']] >= ancestor_lower) &
        (gn_node_table[['upper_age']] >= ancestor_upper)
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
    adjusted_nodes = data.frame(
        node=integer(),
        lower_before=numeric(),
        upper_before=numeric(),
        lower_after=numeric(),
        upper_after=numeric(),
        ancestor_upper=numeric(),
        stringsAsFactors=FALSE
    )
    if (nrow(gn_node_table)==0) {
        return(list(gn_node_table=gn_node_table, adjusted_nodes=adjusted_nodes, min_delta=0))
    }
    max_constraint_age = suppressWarnings(max(gn_node_table[['upper_age']], na.rm=TRUE))
    if ((!is.finite(max_constraint_age)) || is.na(max_constraint_age)) {
        max_constraint_age = 1
    }
    if (is.null(min_delta)) {
        min_delta = max(1e-8, max_constraint_age * 1e-8)
    }
    gn_node_table2 = gn_node_table
    node_nums = as.integer(gn_node_table2[['gn_node_num']])
    total_nodes = length(gn_tree[['tip.label']]) + gn_tree[['Nnode']]
    row_by_node = integer(total_nodes)
    row_by_node[node_nums] = seq_along(node_nums)
    ancestor_lower_by_node = rep(Inf, total_nodes)
    ancestor_upper_by_node = rep(Inf, total_nodes)
    adjusted_rows = list()
    edge_preorder = ape::reorder.phylo(gn_tree, order='cladewise')[['edge']]

    for (edge_i in seq_len(nrow(edge_preorder))) {
        parent_num = edge_preorder[edge_i,1]
        gn_node_num = edge_preorder[edge_i,2]
        ancestor_lower_by_node[[gn_node_num]] = ancestor_lower_by_node[[parent_num]]
        ancestor_upper_by_node[[gn_node_num]] = ancestor_upper_by_node[[parent_num]]
        parent_row = row_by_node[[parent_num]]
        if (parent_row > 0L) {
            ancestor_lower_by_node[[gn_node_num]] = min(
                ancestor_lower_by_node[[gn_node_num]],
                gn_node_table2[["lower_age"]][[parent_row]]
            )
            ancestor_upper_by_node[[gn_node_num]] = min(
                ancestor_upper_by_node[[gn_node_num]],
                gn_node_table2[["upper_age"]][[parent_row]]
            )
        }
        row_idx = row_by_node[[gn_node_num]]
        if ((row_idx == 0L) || (gn_node_num == root_num)) {
            next
        }
        ancestor_lower = ancestor_lower_by_node[[gn_node_num]]
        ancestor_upper = ancestor_upper_by_node[[gn_node_num]]
        if ((!is.finite(ancestor_lower)) || (!is.finite(ancestor_upper))) {
            next
        }
        lower_before = gn_node_table2[row_idx,'lower_age']
        upper_before = gn_node_table2[row_idx,'upper_age']
        lower_after = lower_before
        upper_after = upper_before

        if (upper_after>=ancestor_upper) {
            if (ancestor_upper>min_delta) {
                upper_after = ancestor_upper - min_delta
            } else {
                upper_after = ancestor_upper * 0.5
            }
        }
        if (lower_after>=ancestor_lower) {
            if (ancestor_lower>min_delta) {
                lower_after = ancestor_lower - min_delta
            } else {
                lower_after = ancestor_lower * 0.5
            }
        }
        if (lower_after>=upper_after) {
            lower_after = max(0, upper_after - min_delta)
        }
        if ((upper_after-lower_after)<min_delta) {
            lower_after = max(0, upper_after - min_delta)
        }

        is_changed = (!isTRUE(all.equal(lower_before, lower_after))) || (!isTRUE(all.equal(upper_before, upper_after)))
        if (is_changed) {
            gn_node_table2[row_idx,'lower_age'] = lower_after
            gn_node_table2[row_idx,'upper_age'] = upper_after
            adjusted_rows[[length(adjusted_rows) + 1L]] = data.frame(
                node=gn_node_num,
                lower_before=lower_before,
                upper_before=upper_before,
                lower_after=lower_after,
                upper_after=upper_after,
                ancestor_upper=ancestor_upper,
                stringsAsFactors=FALSE
            )
        }
    }
    if (length(adjusted_rows) > 0L) {
        adjusted_nodes = do.call(rbind, adjusted_rows)
    }
    return(list(gn_node_table=gn_node_table2, adjusted_nodes=adjusted_nodes, min_delta=min_delta))
}

expand_narrow_calibration_ranges = function(calibration_table, root_num, min_span=NULL) {
    if (nrow(calibration_table)==0) {
        return(list(calibration_table=calibration_table, adjusted_nodes=integer(0), min_span=0))
    }
    calibration_table2 = calibration_table
    if (!('soft.bounds' %in% colnames(calibration_table2))) {
        calibration_table2$soft.bounds = NA
    }
    max_cal_age = suppressWarnings(max(calibration_table2[['age.max']], na.rm=TRUE))
    if ((!is.finite(max_cal_age)) || is.na(max_cal_age)) {
        max_cal_age = 1
    }
    if (is.null(min_span)) {
        min_span = max(1e-8, max_cal_age * 1e-8)
    }
    adjusted_nodes = c()
    for (i in seq_len(nrow(calibration_table2))) {
        node_i = calibration_table2$node[i]
        if (node_i==root_num) {
            next
        }
        age_span = calibration_table2$age.max[i] - calibration_table2$age.min[i]
        if ((!is.finite(age_span)) || is.na(age_span) || (age_span < min_span)) {
            calibration_table2$age.min[i] = max(0, calibration_table2$age.max[i] - min_span)
            calibration_table2$soft.bounds[i] = TRUE
            adjusted_nodes = c(adjusted_nodes, node_i)
        }
    }
    return(
        list(
            calibration_table=calibration_table2,
            adjusted_nodes=unique(adjusted_nodes),
            min_span=min_span
        )
    )
}

enforce_descendant_calibration_margin = function(calibration_table, phy, root_num, min_margin) {
    if (nrow(calibration_table)==0) {
        return(list(calibration_table=calibration_table, adjusted_nodes=integer(0), min_margin=min_margin))
    }
    calibration_table2 = calibration_table
    node_nums = as.integer(calibration_table2$node)
    total_nodes = length(phy[['tip.label']]) + phy[['Nnode']]
    row_by_node = integer(total_nodes)
    row_by_node[node_nums] = seq_along(node_nums)
    ancestor_lower_by_node = rep(Inf, total_nodes)
    ancestor_upper_by_node = rep(Inf, total_nodes)
    adjusted_nodes = c()
    edge_preorder = ape::reorder.phylo(phy, order='cladewise')[['edge']]
    for (edge_i in seq_len(nrow(edge_preorder))) {
        parent_num = edge_preorder[edge_i,1]
        node_i = edge_preorder[edge_i,2]
        ancestor_lower_by_node[[node_i]] = ancestor_lower_by_node[[parent_num]]
        ancestor_upper_by_node[[node_i]] = ancestor_upper_by_node[[parent_num]]
        parent_row = row_by_node[[parent_num]]
        if (parent_row > 0L) {
            ancestor_lower_by_node[[node_i]] = min(
                ancestor_lower_by_node[[node_i]], calibration_table2$age.min[[parent_row]]
            )
            ancestor_upper_by_node[[node_i]] = min(
                ancestor_upper_by_node[[node_i]], calibration_table2$age.max[[parent_row]]
            )
        }
        row_idx = row_by_node[[node_i]]
        if ((row_idx == 0L) || (node_i == root_num)) {
            next
        }
        ancestor_upper = ancestor_upper_by_node[[node_i]]
        ancestor_lower = ancestor_lower_by_node[[node_i]]
        if ((!is.finite(ancestor_lower)) || (!is.finite(ancestor_upper))) {
            next
        }

        age_min_before = calibration_table2$age.min[row_idx]
        age_max_before = calibration_table2$age.max[row_idx]
        age_max_after = min(age_max_before, max(0, ancestor_upper - min_margin))
        age_min_after = min(age_min_before, max(0, ancestor_lower - min_margin))
        if (age_min_after >= age_max_after) {
            age_min_after = max(0, age_max_after - min_margin)
        }
        if ((age_max_after - age_min_after) < min_margin) {
            age_min_after = max(0, age_max_after - min_margin)
        }

        is_changed = (!isTRUE(all.equal(age_min_before, age_min_after))) || (!isTRUE(all.equal(age_max_before, age_max_after)))
        if (is_changed) {
            calibration_table2$age.min[row_idx] = age_min_after
            calibration_table2$age.max[row_idx] = age_max_after
            adjusted_nodes = c(adjusted_nodes, node_i)
        }
    }
    return(
        list(
            calibration_table=calibration_table2,
            adjusted_nodes=unique(adjusted_nodes),
            min_margin=min_margin
        )
    )
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
            silent=TRUE
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
        current_seed = as.integer(seed_base + attempt_i - 1)
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
    phy,
    calibration_table,
    root_num,
    chronos_control,
    chronos_lambda,
    chronos_model,
    soft_attempts,
    calibration_label='RS',
    max_restarts_main=3,
    max_restarts_fallback=2,
    seed_cursor=1,
    attempt_timeout_sec=Inf,
    time_budget=NULL
) {
    if (nrow(calibration_table) == 0) {
        return(
            list(
                chronos_out=structure(
                    paste0('No calibration nodes were available for ', calibration_label, '.'),
                    class='try-error'
                ),
                success=FALSE,
                used_model=chronos_model,
                used_lambda=chronos_lambda,
                used_seed=NA_integer_,
                seed_cursor=seed_cursor,
                phy=phy,
                calibration_table=calibration_table,
                skipped=TRUE
            )
        )
    }

    chronos_out = structure('PLACEHOLDER', class='try-error')
    working_phy = phy
    working_calibration = calibration_table
    used_model = chronos_model
    used_lambda = chronos_lambda
    used_seed = NA_integer_

    risk_profile = detect_chronos_failure_risks(working_phy, working_calibration, root_num)
    if (isTRUE(risk_profile$risk_flags[['extreme_edge_ratio']])) {
        cat(
            'Potential chronos failure risk (',
            calibration_label,
            '): extreme edge-length range (max/min =',
            signif(risk_profile$edge_ratio, digits=4),
            ').\n',
            sep=''
        )
    }
    if (isTRUE(risk_profile$risk_flags[['large_edge_values']])) {
        cat(
            'Potential chronos failure risk (',
            calibration_label,
            '): large edge length (max =',
            signif(risk_profile$max_edge, digits=4),
            ').\n',
            sep=''
        )
    }
    if (isTRUE(risk_profile$risk_flags[['tight_nonroot_constraints']])) {
        cat(
            'Potential chronos failure risk (',
            calibration_label,
            '): tight non-root constraints detected at node number(s): ',
            format_limited_values(risk_profile$tight_nodes, max_items=80),
            '\n',
            sep=''
        )
    }

    cat("\nchronos, calibrated nodes: ", calibration_label, "\n", sep='')
    main_out = run_chronos_with_restarts(
        phy=working_phy,
        calibration=working_calibration,
        chronos_control=chronos_control,
        chronos_lambda=chronos_lambda,
        chronos_model=chronos_model,
        context_label=paste0('chronos ', calibration_label),
        max_restarts=max_restarts_main,
        seed_base=seed_cursor,
        attempt_timeout_sec=attempt_timeout_sec,
        time_budget=time_budget
    )
    seed_cursor = seed_cursor + max_restarts_main
    chronos_out = main_out$chronos_out
    if (main_out$success) {
        used_model = main_out$used_model
        used_lambda = main_out$used_lambda
        used_seed = main_out$used_seed
    }

    if (("try-error" %in% class(chronos_out)) &&
        (isTRUE(risk_profile$risk_flags[['extreme_edge_ratio']]) || isTRUE(risk_profile$risk_flags[['large_edge_values']]))) {
        working_phy_retry = normalize_edge_length_range(working_phy, max_edge = 100, min_edge = 1e-6)
        cat("\nchronos, calibrated nodes: ", calibration_label, " (rescaled retry)\n", sep='')
        rescaled_out = run_chronos_with_restarts(
            phy=working_phy_retry,
            calibration=working_calibration,
            chronos_control=chronos_control,
            chronos_lambda=chronos_lambda,
            chronos_model=chronos_model,
            context_label=paste0('chronos ', calibration_label, ' rescaled'),
            max_restarts=max_restarts_main,
            seed_base=seed_cursor,
            attempt_timeout_sec=attempt_timeout_sec,
            time_budget=time_budget
        )
        seed_cursor = seed_cursor + max_restarts_main
        chronos_out = rescaled_out$chronos_out
        if (rescaled_out$success) {
            working_phy = working_phy_retry
            used_model = rescaled_out$used_model
            used_lambda = rescaled_out$used_lambda
            used_seed = rescaled_out$used_seed
        }
    }

    if ("try-error" %in% class(chronos_out)) {
        expanded_calibration = expand_narrow_calibration_ranges(working_calibration, root_num)
        if (length(expanded_calibration$adjusted_nodes) > 0) {
            # Keep widened bounds for downstream soft-bound retries.
            working_calibration = expanded_calibration$calibration_table
            cat("\nchronos, calibrated nodes: ", calibration_label, " (expanded-bound retry)\n", sep='')
            expanded_out = run_chronos_with_restarts(
                phy=working_phy,
                calibration=working_calibration,
                chronos_control=chronos_control,
                chronos_lambda=chronos_lambda,
                chronos_model=chronos_model,
                context_label=paste0('chronos ', calibration_label, ' expanded'),
                max_restarts=max_restarts_main,
                seed_base=seed_cursor,
                attempt_timeout_sec=attempt_timeout_sec,
                time_budget=time_budget
            )
            seed_cursor = seed_cursor + max_restarts_main
            chronos_out = expanded_out$chronos_out
            if (expanded_out$success) {
                used_model = expanded_out$used_model
                used_lambda = expanded_out$used_lambda
                used_seed = expanded_out$used_seed
            }
        }
    }

    if ("try-error" %in% class(chronos_out)) {
        soft_calibration_table = make_soft_bounds_for_nonroot(working_calibration, root_num)
        cat("\nchronos, calibrated nodes: ", calibration_label, " (soft-bound retry)\n", sep='')
        for (sa in soft_attempts) {
            if ("try-error" %in% class(chronos_out)) {
                cat(
                    'chronos soft-bound strategy (',
                    calibration_label,
                    '): ',
                    sa[['label']],
                    '(model=',
                    sa[['model']],
                    ', lambda=',
                    sa[['lambda']],
                    ')\n',
                    sep=''
                )
                soft_out = run_chronos_with_restarts(
                    phy=working_phy,
                    calibration=soft_calibration_table,
                    chronos_control=chronos_control,
                    chronos_lambda=sa[['lambda']],
                    chronos_model=sa[['model']],
                    context_label=paste0('chronos ', calibration_label, ' soft ', sa[['label']]),
                    max_restarts=max_restarts_fallback,
                    seed_base=seed_cursor,
                    attempt_timeout_sec=attempt_timeout_sec,
                    time_budget=time_budget
                )
                seed_cursor = seed_cursor + max_restarts_fallback
                chronos_out = soft_out$chronos_out
                if (soft_out$success) {
                    working_calibration = soft_calibration_table
                    used_model = soft_out$used_model
                    used_lambda = soft_out$used_lambda
                    used_seed = soft_out$used_seed
                }
            }
        }
    }

    if ("try-error" %in% class(chronos_out)) {
        max_cal_age = suppressWarnings(max(working_calibration$age.max, na.rm=TRUE))
        if ((!is.finite(max_cal_age)) || is.na(max_cal_age)) {
            max_cal_age = 1
        }
        aggressive_min_spans = unique(c(
            max(1e-6, max_cal_age * 1e-6),
            max(1e-4, max_cal_age * 1e-4),
            max(1e-3, max_cal_age * 1e-3),
            max(1e-2, max_cal_age * 1e-2)
        ))
        aggressive_min_spans = sort(aggressive_min_spans)
        for (min_span_i in aggressive_min_spans) {
            if ("try-error" %in% class(chronos_out)) {
                marginized = enforce_descendant_calibration_margin(
                    calibration_table=working_calibration,
                    phy=working_phy,
                    root_num=root_num,
                    min_margin=min_span_i
                )
                widened = expand_narrow_calibration_ranges(
                    calibration_table=marginized$calibration_table,
                    root_num=root_num,
                    min_span=min_span_i
                )
                aggressive_calibration_table = make_soft_bounds_for_nonroot(widened$calibration_table, root_num)
                cat(
                    "\nchronos, calibrated nodes: ",
                    calibration_label,
                    " (aggressive soft retry, min_span=",
                    signif(min_span_i, digits=4),
                    ")\n",
                    sep=''
                )
                for (sa in soft_attempts) {
                    if ("try-error" %in% class(chronos_out)) {
                        cat(
                            'chronos aggressive soft strategy (',
                            calibration_label,
                            '): ',
                            sa[['label']],
                            '(model=',
                            sa[['model']],
                            ', lambda=',
                            sa[['lambda']],
                            ')\n',
                            sep=''
                        )
                        aggressive_out = run_chronos_with_restarts(
                            phy=working_phy,
                            calibration=aggressive_calibration_table,
                            chronos_control=chronos_control,
                            chronos_lambda=sa[['lambda']],
                            chronos_model=sa[['model']],
                            context_label=paste0(
                                'chronos ',
                                calibration_label,
                                ' aggressive ',
                                sa[['label']],
                                ' span',
                                format(signif(min_span_i, digits=3), scientific=TRUE)
                            ),
                            max_restarts=max_restarts_fallback,
                            seed_base=seed_cursor,
                            attempt_timeout_sec=attempt_timeout_sec,
                            time_budget=time_budget
                        )
                        seed_cursor = seed_cursor + max_restarts_fallback
                        chronos_out = aggressive_out$chronos_out
                        if (aggressive_out$success) {
                            working_calibration = aggressive_calibration_table
                            used_model = aggressive_out$used_model
                            used_lambda = aggressive_out$used_lambda
                            used_seed = aggressive_out$used_seed
                        }
                    }
                }
            }
        }
    }

    return(
        list(
            chronos_out=chronos_out,
            success=!("try-error" %in% class(chronos_out)),
            used_model=used_model,
            used_lambda=used_lambda,
            used_seed=used_seed,
            seed_cursor=seed_cursor,
            phy=working_phy,
            calibration_table=working_calibration,
            skipped=FALSE
        )
    )
}

make_soft_bounds_for_nonroot = function(calibration_table, root_num) {
    calibration2 = calibration_table
    if (!('soft.bounds' %in% colnames(calibration2))) {
        calibration2$soft.bounds = NA
    }
    calibration2$soft.bounds[calibration2$node != root_num] = TRUE
    return(calibration2)
}

build_dated_tree_without_chronos = function(phy, calibration_table, root_num, min_edge=NULL) {
    if (nrow(calibration_table)==0) {
        return(structure('No calibration constraints were provided.', class='try-error'))
    }
    max_cal_age = suppressWarnings(max(calibration_table$age.max, na.rm=TRUE))
    if ((!is.finite(max_cal_age)) || is.na(max_cal_age) || (max_cal_age <= 0)) {
        max_cal_age = 1
    }
    if (is.null(min_edge)) {
        min_edge = max(1e-8, max_cal_age * 1e-8)
    }
    root_row = which(calibration_table$node==root_num)
    if (length(root_row)==0) {
        return(structure('Root calibration is missing for deterministic fallback.', class='try-error'))
    }
    root_target = calibration_table$age.max[root_row[1]]
    if ((!is.finite(root_target)) || is.na(root_target) || (root_target <= 0)) {
        return(structure('Root calibration age is invalid for deterministic fallback.', class='try-error'))
    }
    node_count = length(phy$tip.label) + phy$Nnode
    depth_values = ape::node.depth.edgelength(phy)
    if (any(!is.finite(depth_values))) {
        return(structure('Deterministic fallback could not compute finite node depths.', class='try-error'))
    }
    baseline_ages = max(depth_values) - depth_values
    max_baseline_age = max(baseline_ages)
    if (max_baseline_age > 0) {
        baseline_ages = baseline_ages * (root_target / max_baseline_age)
    }
    node_ages = baseline_ages
    node_ages[seq_len(length(phy$tip.label))] = 0

    lower_bounds = rep(0, node_count)
    upper_bounds = rep(Inf, node_count)
    for (i in seq_len(nrow(calibration_table))) {
        node_i = as.integer(calibration_table$node[i])
        if ((node_i < 1) || (node_i > node_count)) {
            return(structure(paste('Calibration node is out of range:', node_i), class='try-error'))
        }
        lower_bounds[node_i] = calibration_table$age.min[i]
        upper_bounds[node_i] = calibration_table$age.max[i]
        node_ages[node_i] = max(lower_bounds[node_i], min(upper_bounds[node_i], node_ages[node_i]))
    }
    node_ages[root_num] = root_target
    lower_bounds[root_num] = max(lower_bounds[root_num], root_target)
    upper_bounds[root_num] = min(upper_bounds[root_num], root_target)

    node_nums = seq_len(node_count)
    parent_map = get_parent_map(phy)
    depth_map = get_node_depth_map(phy)
    process_order_up = node_nums[order(depth_map[node_nums], decreasing=FALSE)]
    process_order_down = rev(process_order_up)
    max_iter = node_count * 4
    for (iter_i in seq_len(max_iter)) {
        changed = FALSE
        for (child_node in process_order_down) {
            if (child_node==root_num) {
                next
            }
            parent_node = parent_map[[child_node]]
            if (is.na(parent_node)) {
                next
            }
            required_parent_age = node_ages[child_node] + min_edge
            parent_target_age = max(node_ages[parent_node], required_parent_age, lower_bounds[parent_node])
            parent_target_age = min(parent_target_age, upper_bounds[parent_node])
            if (parent_target_age + min_edge < required_parent_age) {
                feasible_child_age = max(0, parent_target_age - min_edge)
                lower_bounds[child_node] = min(lower_bounds[child_node], feasible_child_age)
                node_ages[child_node] = min(node_ages[child_node], feasible_child_age)
                required_parent_age = node_ages[child_node] + min_edge
                parent_target_age = max(node_ages[parent_node], required_parent_age, lower_bounds[parent_node])
                parent_target_age = min(parent_target_age, upper_bounds[parent_node])
                if (parent_target_age + min_edge < required_parent_age) {
                    return(
                        structure(
                            paste(
                                'Deterministic fallback cannot satisfy parent bound for child node',
                                child_node
                            ),
                            class='try-error'
                        )
                    )
                }
            }
            if (!isTRUE(all.equal(parent_target_age, node_ages[parent_node]))) {
                node_ages[parent_node] = parent_target_age
                changed = TRUE
            }
        }
        for (parent_node in process_order_up) {
            child_rows = which(phy$edge[,1]==parent_node)
            if (length(child_rows)==0) {
                next
            }
            for (edge_idx in child_rows) {
                child_node = phy$edge[edge_idx,2]
                max_child_age = max(0, node_ages[parent_node] - min_edge)
                child_target_age = min(node_ages[child_node], max_child_age, upper_bounds[child_node])
                child_target_age = max(child_target_age, lower_bounds[child_node])
                if (child_target_age > (max_child_age + min_edge)) {
                    lower_bounds[child_node] = min(lower_bounds[child_node], max_child_age)
                    child_target_age = max(0, max_child_age)
                }
                if (!isTRUE(all.equal(child_target_age, node_ages[child_node]))) {
                    node_ages[child_node] = child_target_age
                    changed = TRUE
                }
            }
        }
        node_ages[seq_len(length(phy$tip.label))] = 0
        if (!changed) {
            break
        }
    }

    constrained_nodes = as.integer(calibration_table$node)
    for (node_i in constrained_nodes) {
        if (node_ages[node_i] < (lower_bounds[node_i] - min_edge)) {
            return(
                structure(
                    paste('Deterministic fallback violated lower bound at node', node_i),
                    class='try-error'
                )
            )
        }
        if (node_ages[node_i] > (upper_bounds[node_i] + min_edge)) {
            return(
                structure(
                    paste('Deterministic fallback violated upper bound at node', node_i),
                    class='try-error'
                )
            )
        }
    }

    edge_lengths = node_ages[phy$edge[,1]] - node_ages[phy$edge[,2]]
    if (any(!is.finite(edge_lengths))) {
        return(structure('Deterministic fallback generated non-finite edge lengths.', class='try-error'))
    }
    edge_lengths[edge_lengths < min_edge] = min_edge
    phy2 = phy
    phy2$edge.length = edge_lengths
    return(phy2)
}
