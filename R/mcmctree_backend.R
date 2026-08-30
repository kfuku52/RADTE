format_mcmctree_number = function(value) {
    value_num = suppressWarnings(as.numeric(value))
    if (is.na(value_num) || (!is.finite(value_num))) {
        stop('MCMCTree numeric value is invalid: ', value)
    }
    txt = formatC(value_num, digits=15, format='fg', flag='#')
    txt = sub('\\.$', '', txt)
    return(txt)
}

assert_no_input_output_collision = function(inputs, outputs) {
    inputs = inputs[file.exists(inputs)]
    resolved_inputs = vapply(inputs, normalizePath, character(1), mustWork=TRUE)
    existing = outputs[file.exists(outputs)]
    if (length(existing)) {
        resolved_outputs = vapply(existing, normalizePath, character(1), mustWork=TRUE)
        collisions = existing[resolved_outputs %in% resolved_inputs]
        if (length(collisions)) {
            stop('Input/output path collision; refusing to overwrite input: ', paste(collisions, collapse=', '))
        }
    }
    invisible(TRUE)
}

stage_external_input_file = function(source, staged_file) {
    if (!is_nonempty_scalar_string(source) || !file.exists(source) || dir.exists(source)) {
        stop('External input file was not found: ', source)
    }
    source_abs = normalizePath(source, mustWork=TRUE)
    assert_no_input_output_collision(source_abs, staged_file)
    # Copy a snapshot: the external program must never write through a link to
    # the original alignment. Atomic replacement also handles dangling links.
    atomic_copy_file(source_abs, staged_file)
    invisible(staged_file)
}

make_mcmctree_calibration_text = function(age_min, age_max, exact_ratio=1e-6, exact_min=1e-8) {
    age_min_num = suppressWarnings(as.numeric(age_min))
    age_max_num = suppressWarnings(as.numeric(age_max))
    has_min = !is.na(age_min_num) && is.finite(age_min_num)
    has_max = !is.na(age_max_num) && is.finite(age_max_num)
    if ((!has_min) && (!has_max)) {
        return(NA_character_)
    }
    if (has_min && has_max) {
        if (age_max_num < age_min_num) {
            stop('MCMCTree calibration upper bound is younger than lower bound.')
        }
        if (age_min_num == age_max_num) {
            eps = max(exact_min, abs(age_max_num) * exact_ratio)
            age_min_num = max(0, age_min_num - eps)
            age_max_num = age_max_num + eps
        }
        return(
            paste0(
                'B{',
                format_mcmctree_number(age_min_num),
                ', ',
                format_mcmctree_number(age_max_num),
                '}'
            )
        )
    }
    if (has_min) {
        return(paste0('L(', format_mcmctree_number(age_min_num), ', 0.1, 1, 1e-300)'))
    }
    return(paste0('U(', format_mcmctree_number(age_max_num), ', 1e-300)'))
}

make_mcmctree_root_calibration_text = function(age_min, age_max, exact_ratio=1e-6, exact_min=1e-8) {
    age_min_num = suppressWarnings(as.numeric(age_min))
    age_max_num = suppressWarnings(as.numeric(age_max))
    has_min = !is.na(age_min_num) && is.finite(age_min_num)
    has_max = !is.na(age_max_num) && is.finite(age_max_num)
    if ((!has_min) && (!has_max)) {
        return(NA_character_)
    }
    if (has_min && has_max) {
        if (age_max_num < age_min_num) {
            stop('MCMCTree root calibration upper bound is younger than lower bound.')
        }
        if (age_min_num == age_max_num) {
            eps = max(exact_min, abs(age_max_num) * exact_ratio)
            age_min_num = max(0, age_min_num - eps)
            age_max_num = age_max_num + eps
        }
        return(
            paste0(
                '>',
                format_mcmctree_number(age_min_num),
                '<',
                format_mcmctree_number(age_max_num)
            )
        )
    }
    if (has_min) {
        return(paste0('>', format_mcmctree_number(age_min_num)))
    }
    return(paste0('<', format_mcmctree_number(age_max_num)))
}

validate_mcmctree_calibration_constraints = function(phy, calibration_table) {
    required_cols = c('node', 'age.min', 'age.max')
    if (any(!required_cols %in% colnames(calibration_table))) {
        stop(
            'MCMCTree calibration table is missing required column(s): ',
            paste(setdiff(required_cols, colnames(calibration_table)), collapse=', ')
        )
    }
    if (nrow(calibration_table) == 0) {
        return(invisible(TRUE))
    }

    node_nums = suppressWarnings(as.integer(calibration_table[['node']]))
    age_min = suppressWarnings(as.numeric(calibration_table[['age.min']]))
    age_max = suppressWarnings(as.numeric(calibration_table[['age.max']]))
    total_nodes = ape::Ntip(phy) + phy$Nnode
    invalid_node = is.na(node_nums) | (node_nums <= ape::Ntip(phy)) | (node_nums > total_nodes)
    if (any(invalid_node)) {
        stop(
            'MCMCTree calibration table contains invalid internal node number(s): ',
            paste(unique(node_nums[invalid_node]), collapse=', ')
        )
    }
    if (any(duplicated(node_nums))) {
        stop(
            'MCMCTree calibration table contains duplicated node(s): ',
            paste(unique(node_nums[duplicated(node_nums)]), collapse=', ')
        )
    }
    invalid_bounds = is.na(age_min) | (!is.finite(age_min)) | is.na(age_max) | (!is.finite(age_max))
    if (any(invalid_bounds)) {
        stop(
            'MCMCTree calibration table contains missing or non-finite bound(s) for node(s): ',
            paste(node_nums[invalid_bounds], collapse=', ')
        )
    }
    if (any(age_min < 0) || any(age_max < 0)) {
        stop('MCMCTree calibration table contains negative age bound(s).')
    }
    reversed_bounds = age_max < age_min
    if (any(reversed_bounds)) {
        stop(
            'MCMCTree calibration table contains age.max younger than age.min for node(s): ',
            paste(node_nums[reversed_bounds], collapse=', ')
        )
    }

    root_num = get_root_num(phy)
    ancestor_context = get_ancestor_constraint_minima(phy, node_nums, age_max)
    ancestor_upper = ancestor_context[['minimum']][node_nums]
    invalid_temporal = (node_nums != root_num) & is.finite(ancestor_upper) &
        (age_min >= ancestor_upper)
    if (any(invalid_temporal)) {
        i = which(invalid_temporal)[[1]]
        child_node = node_nums[[i]]
        ancestor_node = ancestor_context[['source']][[child_node]]
        ancestor_i = match(ancestor_node, node_nums)
        child_name = get_node_name_by_num(phy, child_node)
        ancestor_name = get_node_name_by_num(phy, ancestor_node)
        stop(
            'MCMCTree calibration bounds are temporally infeasible: descendant ',
            child_name,
            ' (node ',
            child_node,
            ', age.min=',
            age_min[[i]],
            ') is not younger than ancestor ',
            ancestor_name,
            ' (node ',
            ancestor_node,
            ', age.max=',
            age_max[[ancestor_i]],
            '). Bounds are preserved; revise the input constraints.'
        )
    }
    return(invisible(TRUE))
}

empty_mcmctree_mirror_table = function() {
    return(
        data.frame(
            shared_speciation_group=character(),
            node=integer(),
            gn_node=character(),
            mirror_label=character(),
            mirror_role=character(),
            prior_emitted=logical(),
            stringsAsFactors=FALSE
        )
    )
}

build_mcmctree_annotation_map = function(phy, gn_node_table, calibration_table, root_num) {
    total_nodes = ape::Ntip(phy) + phy$Nnode
    annotation_map = stats::setNames(as.list(rep('', total_nodes)), as.character(seq_len(total_nodes)))
    validate_mcmctree_calibration_constraints(phy, calibration_table)
    if (nrow(gn_node_table) == 0) {
        return(
            list(
                annotation_map=annotation_map,
                duplication_flag=0L,
                mirror_table=empty_mcmctree_mirror_table()
            )
        )
    }

    calibration_text_by_node = rep(NA_character_, total_nodes)
    if (nrow(calibration_table) > 0) {
        for (i in seq_len(nrow(calibration_table))) {
            node_i = as.integer(calibration_table$node[i])
            if (node_i == root_num) {
                calibration_text_by_node[[node_i]] = make_mcmctree_root_calibration_text(
                    age_min=calibration_table$age.min[i],
                    age_max=calibration_table$age.max[i]
                )
            } else {
                calibration_text_by_node[[node_i]] = make_mcmctree_calibration_text(
                    age_min=calibration_table$age.min[i],
                    age_max=calibration_table$age.max[i]
                )
            }
        }
    }

    mirror_label_by_node = rep(NA_character_, total_nodes)
    mirror_group = rep(NA_character_, nrow(gn_node_table))
    if ('shared_speciation_group' %in% colnames(gn_node_table)) {
        mirror_group = as.character(gn_node_table[['shared_speciation_group']])
    } else if ('constraint_sp_node' %in% colnames(gn_node_table)) {
        mirror_group = as.character(gn_node_table[['constraint_sp_node']])
    } else if (all(c('lower_sp_node', 'upper_sp_node') %in% colnames(gn_node_table))) {
        is_same_sp_node = (!is.na(gn_node_table[['lower_sp_node']])) &
            (!is.na(gn_node_table[['upper_sp_node']])) &
            (gn_node_table[['lower_sp_node']] == gn_node_table[['upper_sp_node']])
        mirror_group[is_same_sp_node] = as.character(gn_node_table[['lower_sp_node']][is_same_sp_node])
    }
    spec_rows = gn_node_table[
        (gn_node_table$gn_node_num > ape::Ntip(phy)) &
        grepl('^S', gn_node_table$event) &
        !is.na(mirror_group) &
        (mirror_group != ''),
        ,
        drop=FALSE
    ]
    if (nrow(spec_rows) > 0) {
        spec_rows[['mirror_group']] = mirror_group[match(spec_rows[['gn_node_num']], gn_node_table[['gn_node_num']])]
    }
    mirror_table = empty_mcmctree_mirror_table()
    if (nrow(spec_rows) > 0) {
        parent_map = get_parent_map(phy)
        spec_groups = split(spec_rows, spec_rows$mirror_group)
        label_counter = 1L
        for (group_name in names(spec_groups)) {
            group_df = spec_groups[[group_name]]
            if (nrow(group_df) < 2) {
                stop(
                    'MCMCTree shared speciation group ',
                    group_name,
                    ' has fewer than two member nodes.'
                )
            }
            node_nums = sort(unique(as.integer(group_df$gn_node_num)))
            for (node_i in node_nums) {
                other_nodes = setdiff(node_nums, node_i)
                ancestral_members = intersect(
                    get_ancestor_num(phy, node_i, parent_map, root_num),
                    other_nodes
                )
                if (length(ancestral_members) > 0) {
                    stop(
                        'MCMCTree shared speciation group ',
                        group_name,
                        ' contains an ancestor-descendant pair: ',
                        get_node_name_by_num(phy, node_i),
                        ' (node ',
                        node_i,
                        ') and ',
                        get_node_name_by_num(phy, ancestral_members[[1]]),
                        ' (node ',
                        ancestral_members[[1]],
                        '). Shared ages would create a zero-length path.'
                    )
                }
            }
            group_rows = calibration_table[calibration_table$node %in% node_nums, , drop=FALSE]
            has_complete_group = (nrow(group_rows) == length(node_nums))
            if (!has_complete_group) {
                missing_nodes = setdiff(node_nums, as.integer(group_rows$node))
                missing_names = vapply(
                    missing_nodes,
                    function(node_i) paste0(get_node_name_by_num(phy, node_i), ' (node ', node_i, ')'),
                    character(1)
                )
                stop(
                    'MCMCTree shared speciation group ',
                    group_name,
                    ' is missing calibration row(s) for: ',
                    paste(missing_names, collapse=', '),
                    '.'
                )
            }
            group_rows = group_rows[match(node_nums, as.integer(group_rows$node)), , drop=FALSE]
            group_min = as.numeric(group_rows$age.min)
            group_max = as.numeric(group_rows$age.max)
            bounds_scale = max(1, abs(c(group_min, group_max)))
            bounds_tolerance = bounds_scale * 1e-10
            is_same_bounds = (
                all(abs(group_min - group_min[[1]]) <= bounds_tolerance) &&
                all(abs(group_max - group_max[[1]]) <= bounds_tolerance)
            )
            if (!is_same_bounds) {
                bounds_text = vapply(
                    seq_along(node_nums),
                    function(i) {
                        paste0(
                            get_node_name_by_num(phy, node_nums[[i]]),
                            ' (node ',
                            node_nums[[i]],
                            ')=[',
                            format_mcmctree_number(group_min[[i]]),
                            ',',
                            format_mcmctree_number(group_max[[i]]),
                            ']'
                        )
                    },
                    character(1)
                )
                stop(
                    'MCMCTree shared speciation group ',
                    group_name,
                    ' has inconsistent calibration bounds: ',
                    paste(bounds_text, collapse='; '),
                    '.'
                )
            }
            label_text = paste0('#', label_counter)
            mirror_label_by_node[node_nums] = label_text
            keep_node = min(node_nums)
            drop_nodes = setdiff(node_nums, keep_node)
            if (length(drop_nodes) > 0) {
                calibration_text_by_node[drop_nodes] = NA_character_
            }
            group_node_names = vapply(
                node_nums,
                function(node_i) get_node_name_by_num(phy, node_i),
                character(1)
            )
            mirror_table = rbind(
                mirror_table,
                data.frame(
                    shared_speciation_group=rep(group_name, length(node_nums)),
                    node=node_nums,
                    gn_node=group_node_names,
                    mirror_label=rep(label_text, length(node_nums)),
                    mirror_role=ifelse(node_nums == keep_node, 'driver', 'mirror'),
                    prior_emitted=(node_nums == keep_node),
                    stringsAsFactors=FALSE
                )
            )
            label_counter = label_counter + 1L
        }
    }

    internal_nodes = (ape::Ntip(phy) + 1):total_nodes
    for (node_i in internal_nodes) {
        label_i = mirror_label_by_node[[node_i]]
        cal_i = calibration_text_by_node[[node_i]]
        if (node_i == root_num) {
            pieces = c(label_i, cal_i)
            pieces = pieces[!is.na(pieces) & (nchar(pieces) > 0)]
            if (length(pieces) > 0) {
                annotation_map[[as.character(node_i)]] = paste0(' ', paste(pieces, collapse=' '))
            }
            next
        }
        if (!is.na(cal_i) && (nchar(cal_i) > 0)) {
            pieces = c(label_i, cal_i)
            pieces = pieces[!is.na(pieces) & (nchar(pieces) > 0)]
            annotation_map[[as.character(node_i)]] = paste0(' [', paste(pieces, collapse=' '), ']')
            next
        }
        if (!is.na(label_i) && (nchar(label_i) > 0)) {
            annotation_map[[as.character(node_i)]] = paste0(' ', label_i)
        }
    }

    duplication_flag = as.integer(any(!is.na(mirror_label_by_node) & (nchar(mirror_label_by_node) > 0)))
    return(
        list(
            annotation_map=annotation_map,
            duplication_flag=duplication_flag,
            mirror_table=mirror_table
        )
    )
}

build_mcmctree_tree_text = function(phy, gn_node_table, calibration_table, root_num) {
    annotation_info = build_mcmctree_annotation_map(phy, gn_node_table, calibration_table, root_num)
    annotation_map = annotation_info$annotation_map

    total = length(phy$tip.label) + phy$Nnode
    annotations = unname(annotation_map[as.character(seq_len(total))])
    children = split(phy$edge[, 2], factor(phy$edge[, 1], levels=seq_len(total)))
    stack = integer(total * 3L)
    stack[[1]] = root_num
    top = 1L
    tokens = character(total * 3L)
    count = 0L
    while (top > 0L) {
        node = stack[[top]]
        top = top - 1L
        count = count + 1L
        if (node == 0L) {
            tokens[[count]] = ','
        } else if (node < 0L) {
            annotation = annotations[[-node]]
            if (is.null(annotation) || is.na(annotation)) annotation = ''
            tokens[[count]] = paste0(')', annotation)
        } else if (node <= length(phy$tip.label)) {
            tokens[[count]] = phy$tip.label[[node]]
        } else {
            tokens[[count]] = '('
            top = top + 1L
            stack[[top]] = -node
            child_nodes = children[[node]]
            for (i in rev(seq_along(child_nodes))) {
                top = top + 1L
                stack[[top]] = child_nodes[[i]]
                if (i > 1L) {
                    top = top + 1L
                    stack[[top]] = 0L
                }
            }
        }
    }
    list(tree_text=paste0(paste0(tokens[seq_len(count)], collapse=''), ';'),
         duplication_flag=annotation_info$duplication_flag,
         mirror_table=annotation_info$mirror_table)

}

write_mcmctree_tree_file = function(file, phy, gn_node_table, calibration_table, root_num) {
    tree_info = build_mcmctree_tree_text(phy, gn_node_table, calibration_table, root_num)
    writeLines(
        c(
            paste(ape::Ntip(phy), 1),
            tree_info$tree_text
        ),
        con=file
    )
    return(tree_info)
}

read_mcmctree_posterior_tree = function(file) {
    if (!file.exists(file)) {
        stop('MCMCTree output file was not found: ', file)
    }
    lines = readLines(file, warn=FALSE)
    section_idx = grep('Species tree for FigTree', lines, fixed=TRUE)
    if (length(section_idx) == 0) {
        stop('Could not locate the FigTree tree section in the MCMCTree output.')
    }
    trees = character(0)
    current_tree = ''
    for (i in seq.int(section_idx[[1]] + 1, length(lines))) {
        line_i = trimws(lines[[i]])
        if (nchar(line_i) == 0) {
            next
        }
        current_tree = paste0(current_tree, line_i)
        if (grepl(';$', line_i)) {
            trees = c(trees, current_tree)
            current_tree = ''
            if (length(trees) >= 2) {
                break
            }
        }
    }
    if (length(trees) < 2) {
        stop('Could not extract the posterior mean time tree from the MCMCTree output.')
    }
    return(trees[[2]])
}

empty_shared_speciation_age_summary = function() {
    return(
        data.frame(
            species_node=character(),
            gn_nodes=character(),
            member_count=integer(),
            posterior_mean=numeric(),
            posterior_min=numeric(),
            posterior_max=numeric(),
            max_member_age_diff=numeric(),
            tolerance=numeric(),
            stringsAsFactors=FALSE
        )
    )
}

summarize_shared_speciation_ages = function(
    phy,
    gn_node_table,
    tolerance=NULL,
    fail_on_mismatch=TRUE
) {
    required_cols = c('gn_node', 'gn_node_num', 'event', 'shared_speciation_group')
    if ((nrow(gn_node_table) == 0) || any(!required_cols %in% colnames(gn_node_table))) {
        return(empty_shared_speciation_age_summary())
    }
    is_shared = (
        grepl('^S', gn_node_table[['event']]) &
        !is.na(gn_node_table[['shared_speciation_group']]) &
        (gn_node_table[['shared_speciation_group']] != '')
    )
    is_shared[is.na(is_shared)] = FALSE
    shared_rows = gn_node_table[is_shared, , drop=FALSE]
    if (nrow(shared_rows) == 0) {
        return(empty_shared_speciation_age_summary())
    }

    depth_values = ape::node.depth.edgelength(phy)
    if (length(depth_values) != (ape::Ntip(phy) + phy$Nnode) || any(!is.finite(depth_values))) {
        stop('Could not compute finite node ages for MCMCTree shared speciation QA.')
    }
    tree_height = max(depth_values[seq_len(ape::Ntip(phy))])
    if (is.null(tolerance)) {
        tolerance = max(1e-6, abs(tree_height) * 1e-6)
    }
    tolerance = suppressWarnings(as.numeric(tolerance))
    if (is.na(tolerance) || (!is.finite(tolerance)) || (tolerance < 0)) {
        stop('Shared speciation posterior age tolerance should be a non-negative finite number.')
    }

    summary_table = empty_shared_speciation_age_summary()
    shared_groups = split(shared_rows, shared_rows$shared_speciation_group)
    for (group_name in names(shared_groups)) {
        group_df = shared_groups[[group_name]]
        if (nrow(group_df) < 2) {
            stop('Shared speciation posterior QA group ', group_name, ' has fewer than two member nodes.')
        }
        gn_nodes = as.character(group_df[['gn_node']])
        posterior_node_nums = vapply(
            gn_nodes,
            function(node_name) {
                matched_nodes = get_node_num_by_name(phy, node_name)
                if (length(matched_nodes) != 1) {
                    stop(
                        'Could not uniquely map shared speciation node ',
                        node_name,
                        ' into the MCMCTree posterior tree.'
                    )
                }
                return(as.integer(matched_nodes[[1]]))
            },
            integer(1)
        )
        member_ages = tree_height - depth_values[posterior_node_nums]
        age_min = min(member_ages)
        age_max = max(member_ages)
        max_diff = age_max - age_min
        summary_table = rbind(
            summary_table,
            data.frame(
                species_node=group_name,
                gn_nodes=paste(gn_nodes, collapse=','),
                member_count=length(gn_nodes),
                posterior_mean=mean(member_ages),
                posterior_min=age_min,
                posterior_max=age_max,
                max_member_age_diff=max_diff,
                tolerance=tolerance,
                stringsAsFactors=FALSE
            )
        )
        if (isTRUE(fail_on_mismatch) && (max_diff > tolerance)) {
            age_text = paste0(gn_nodes, '=', format(member_ages, digits=15), collapse=', ')
            stop(
                'MCMCTree shared speciation posterior ages differ for group ',
                group_name,
                ': ',
                age_text,
                '. Maximum difference ',
                format(max_diff, digits=15),
                ' exceeds tolerance ',
                format(tolerance, digits=15),
                '.'
            )
        }
    }
    return(summary_table)
}

annotate_calibration_output = function(
    calibration_table,
    gn_node_table,
    mirror_table=NULL,
    emitted_nodes=NULL
) {
    calibration_out = calibration_table
    calibration_out[['.input_order']] = seq_len(nrow(calibration_out))
    metadata_cols = intersect(
        c('gn_node_num', 'gn_node', 'event', 'shared_speciation_group'),
        colnames(gn_node_table)
    )
    node_metadata = gn_node_table[, metadata_cols, drop=FALSE]
    calibration_out = merge(
        calibration_out,
        node_metadata,
        by.x='node',
        by.y='gn_node_num',
        all.x=TRUE,
        sort=FALSE
    )
    calibration_out = calibration_out[order(calibration_out[['.input_order']]), , drop=FALSE]
    calibration_out[['.input_order']] = NULL
    if (!('shared_speciation_group' %in% colnames(calibration_out))) {
        calibration_out[['shared_speciation_group']] = NA_character_
    }
    calibration_out[['mirror_role']] = NA_character_
    calibration_out[['prior_emitted']] = NA

    if (!is.null(mirror_table)) {
        if (is.null(emitted_nodes)) {
            emitted_nodes = as.integer(calibration_out[['node']])
        }
        emitted_nodes = as.integer(emitted_nodes)
        calibration_out[['prior_emitted']] = as.integer(calibration_out[['node']]) %in% emitted_nodes
        if (nrow(mirror_table) > 0) {
            mirror_idx = match(as.integer(calibration_out[['node']]), as.integer(mirror_table[['node']]))
            has_mirror = !is.na(mirror_idx)
            calibration_out[has_mirror, 'shared_speciation_group'] = mirror_table[
                mirror_idx[has_mirror],
                'shared_speciation_group'
            ]
            calibration_out[has_mirror, 'mirror_role'] = mirror_table[
                mirror_idx[has_mirror],
                'mirror_role'
            ]
            calibration_out[has_mirror, 'prior_emitted'] = (
                mirror_table[mirror_idx[has_mirror], 'prior_emitted'] &
                (as.integer(calibration_out[has_mirror, 'node']) %in% emitted_nodes)
            )
        }
    }
    return(calibration_out)
}

write_mcmctree_control_file = function(
    file,
    seqfile_name,
    treefile_name,
    usedata=1,
    seqtype=0,
    clock=2,
    model=0,
    burnin=2000,
    sampfreq=10,
    nsample=20000,
    ncatG=5,
    seed=1L,
    duplication_flag=0
) {
    lines = c(
        paste0('seed = ', parse_integer_arg(seed, '--seed', min_value=1L)),
        paste0('seqfile = ', seqfile_name),
        paste0('treefile = ', treefile_name),
        'outfile = out.txt',
        'mcmcfile = mcmc.txt',
        'ndata = 1',
        paste0('seqtype = ', seqtype),
        paste0('usedata = ', usedata),
        paste0('clock = ', clock),
        paste0('model = ', model),
        'alpha = 0',
        paste0('ncatG = ', ncatG),
        'cleandata = 0',
        'BDparas = 1 1 0.1 M',
        'rgene_gamma = 2 20 1',
        'sigma2_gamma = 1 10 1',
        'kappa_gamma = 6 2',
        'alpha_gamma = 1 1',
        paste0('burnin = ', burnin),
        paste0('sampfreq = ', sampfreq),
        paste0('nsample = ', nsample),
        paste0('duplication = ', duplication_flag),
        'print = 1',
        'finetune = 1: .1 .1 .1 .1 .1 .1'
    )
    writeLines(lines, con=file)
    return(invisible(file))
}

run_mcmctree_backend = function(
    phy,
    gn_node_table,
    calibration_table,
    root_num,
    seqfile,
    bin='mcmctree',
    workdir='radte_mcmctree_run',
    usedata=1,
    seqtype=0,
    clock=2,
    model=0,
    burnin=2000,
    sampfreq=10,
    nsample=20000,
    ncatG=5,
    seed=1L,
    timeout_sec=Inf,
    protected_inputs=character()
) {
    if (!is_nonempty_scalar_string(seqfile)) {
        stop('--mcmctree_seqfile is required when --dating_backend=mcmctree.')
    }
    if (!file.exists(seqfile)) {
        stop('MCMCTree seqfile was not found: ', seqfile)
    }
    bin_path = as.character(bin)
    if (!file.exists(bin_path)) {
        bin_path = Sys.which(bin_path)
    }
    if ((!is_nonempty_scalar_string(bin_path)) || (!file.exists(bin_path)) || dir.exists(bin_path)) {
        stop('MCMCTree executable was not found: ', bin)
    }

    # Keep the final link name: PAML's Debian wrapper dispatches by basename($0).
    bin_path = file.path(normalizePath(dirname(bin_path), mustWork=TRUE), basename(bin_path))
    if (file.access(bin_path, 1L) != 0L) stop('MCMCTree executable is not executable: ', bin_path)
    seqfile = normalizePath(seqfile, mustWork=TRUE)
    tree_info = build_mcmctree_tree_text(phy, gn_node_table, calibration_table, root_num)

    dir.create(workdir, recursive=TRUE, showWarnings=FALSE)
    workdir_abs = normalizePath(workdir, mustWork=TRUE)
    known_outputs = c(
        'seqfile.phy', 'input.trees', 'mcmctree.ctl', 'out.txt', 'mcmc.txt',
        'FigTree.tre', 'mcmctree.stdout.log', 'mcmctree.stderr.log'
    )
    known_output_paths = file.path(workdir_abs, known_outputs)
    assert_no_input_output_collision(c(seqfile, bin_path, protected_inputs), known_output_paths)
    if (any(dir.exists(known_output_paths))) stop('MCMCTree output path is an existing directory.')
    existing_outputs = known_output_paths[file.exists(known_output_paths) | (!is.na(Sys.readlink(known_output_paths)) & nzchar(Sys.readlink(known_output_paths)))]
    staging_dir = tempfile(pattern='.mcmctree-inputs-', tmpdir=dirname(workdir_abs))
    dir.create(staging_dir)
    on.exit(unlink(staging_dir, recursive=TRUE), add=TRUE)
    seqfile_staged = file.path(staging_dir, 'seqfile.phy')
    treefile_staged = file.path(staging_dir, 'input.trees')
    ctlfile_staged = file.path(staging_dir, 'mcmctree.ctl')

    stage_external_input_file(seqfile, seqfile_staged)
    writeLines(c(paste(ape::Ntip(phy), 1), tree_info$tree_text), treefile_staged)
    write_mcmctree_control_file(
        file=ctlfile_staged,
        seqfile_name=basename(seqfile_staged),
        treefile_name=basename(treefile_staged),
        usedata=usedata,
        seqtype=seqtype,
        clock=clock,
        model=model,
        burnin=burnin,
        sampfreq=sampfreq,
        nsample=nsample,
        ncatG=ncatG,
        seed=seed,
        duplication_flag=tree_info$duplication_flag
    )

    if (length(existing_outputs) > 0) {
        unlink(existing_outputs)
    }
    for (input in c(seqfile_staged, treefile_staged, ctlfile_staged)) {
        if (!file.rename(input, file.path(workdir_abs, basename(input)))) stop('Could not stage MCMCTree input: ', input)
    }

    old_wd = getwd()
    on.exit(setwd(old_wd), add=TRUE)
    setwd(workdir_abs)
    system_timeout = if (is.finite(timeout_sec)) max(1L, as.integer(ceiling(timeout_sec))) else 0L
    exit_code = suppressWarnings(
        system2(
            command=bin_path,
            args=basename(ctlfile_staged),
            stdout='mcmctree.stdout.log',
            stderr='mcmctree.stderr.log',
            timeout=system_timeout
        )
    )
    if (exit_code != 0) {
        stderr_tail = character(0)
        stdout_tail = character(0)
        if (file.exists('mcmctree.stderr.log')) {
            stderr_lines = readLines('mcmctree.stderr.log', warn=FALSE)
            stderr_tail = utils::tail(stderr_lines, 20)
        }
        if (file.exists('mcmctree.stdout.log')) {
            stdout_lines = readLines('mcmctree.stdout.log', warn=FALSE)
            stdout_tail = utils::tail(stdout_lines, 20)
        }
        timeout_note = if (is.finite(timeout_sec) && exit_code == 124L) {
            paste0(' after reaching the ', timeout_sec, '-second timeout')
        } else {
            ''
        }
        stop(
            'MCMCTree run failed with exit code ',
            exit_code,
            timeout_note,
            '. stderr tail: ',
            paste(stderr_tail, collapse=' | '),
            '. stdout tail: ',
            paste(stdout_tail, collapse=' | ')
        )
    }

    out_file = file.path(workdir_abs, 'out.txt')
    posterior_tree_text = read_mcmctree_posterior_tree(out_file)
    posterior_tree = ape::read.tree(text=posterior_tree_text)
    if (is.null(posterior_tree$node.label) || (length(posterior_tree$node.label) != posterior_tree$Nnode)) {
        posterior_tree$node.label = paste0('mcmctree_node_', seq_len(posterior_tree$Nnode))
    }
    validate_dated_tree(posterior_tree, phy, tolerance=1e-5)
    posterior_tree = transfer_dated_ages(phy, posterior_tree, tolerance=1e-5)
    shared_speciation_ages = summarize_shared_speciation_ages(
        phy=posterior_tree,
        gn_node_table=gn_node_table,
        fail_on_mismatch=TRUE
    )
    stdout_lines = if (file.exists('mcmctree.stdout.log')) readLines('mcmctree.stdout.log', warn=FALSE) else character()
    mcmctree_banner = paste(utils::head(stdout_lines[nzchar(trimws(stdout_lines))], 3), collapse=' | ')
    return(
        list(
            tree=posterior_tree,
            workdir=workdir_abs,
            duplication_flag=tree_info$duplication_flag,
            mirror_table=tree_info$mirror_table,
            shared_speciation_ages=shared_speciation_ages,
            executable=bin_path,
            banner=mcmctree_banner,
            seed=as.integer(seed)
        )
    )
}
