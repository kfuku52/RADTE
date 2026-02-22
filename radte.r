#!/usr/bin/env Rscript

radte_version = '0.3.5'

#devtools::install_github(repo="cran/ape", ref="master")

library(ape)

cat(paste(version[['version.string']], '\n'))
cat(paste('ape version:', packageVersion('ape'), '\n'))
cat(paste('RADTE version:', radte_version, '\n'))
cat(paste('RADTE command:', paste(commandArgs(trailingOnly=FALSE), collapse=' '), '\n'))
cat(paste('RADTE bug report:', 'https://github.com/kfuku52/RADTE/issues', '\n'))


# copied from rkftools https://github.com/kfuku52/rkftools
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
        param = split[[i]][1]
        if (nchar(param)==0) {
            stop('Argument key is empty in --key=value format: ', args[[i]])
        }
        if (param %in% names(parsed)) {
            stop('Argument key is duplicated: --', param)
        }
        value = paste(split[[i]][2:length(split[[i]])], collapse='=')
        if (!is.na(suppressWarnings(as.numeric(value)))) {
            value = as.numeric(value)
        }
        parsed[[param]] = value
    }
    if (print) {
        for (name in names(parsed)) {
            cat(name, '=', parsed[[name]], '\n')
        }
        cat('\n')
    }
    return(parsed)
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

ensure_root_event_tag = function(event_values) {
    if (length(event_values)==0) {
        return(event_values)
    }
    event_values = as.character(event_values)
    has_root_tag = grepl('\\(R\\)$', event_values) | (event_values=='R')
    event_values[!has_root_tag] = paste0(event_values[!has_root_tag], '(R)')
    return(event_values)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_node_name_by_num = function(phy, node_num) {
    if (length(node_num)==0) {
        return(character(0))
    }
    node_num = suppressWarnings(as.integer(node_num))
    node_num = node_num[!is.na(node_num)]
    if (length(node_num)==0) {
        return(character(0))
    }
    node_names = c(phy[['tip.label']], phy$node.label)
    node_nums = seq_along(node_names)
    node_name = unlist(
        lapply(node_num, function(x) {
            node_names[node_nums == x]
        }),
        use.names=FALSE
    )
    if (is.null(node_name)) {
        return(character(0))
    }
    return(node_name)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_node_num_by_name = function(phy, node_name) {
    if (length(node_name)==0) {
        return(integer(0))
    }
    node_name = as.character(node_name)
    node_name = node_name[!is.na(node_name)]
    if (length(node_name)==0) {
        return(integer(0))
    }
    node_names = c(phy[['tip.label']], phy$node.label)
    node_nums = seq_along(node_names)
    node_num = unlist(
        lapply(node_name, function(x) {
            node_nums[node_names == x]
        }),
        use.names=FALSE
    )
    if (is.null(node_num)) {
        return(integer(0))
    }
    return(node_num)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_root_num = function(phy) {
    root_num = setdiff(phy[['edge']][,1], phy[['edge']][,2])
    return(root_num)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_parent_num = function(phy, node_num) {
    parent_num = phy[['edge']][(phy[['edge']][,2]==node_num),1]
    return(parent_num)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_sister_num = function(phy, node_num) {
    parent_num = phy[['edge']][(phy[['edge']][,2]==node_num),1]
    sibling_num = phy[['edge']][(phy[['edge']][,1]==parent_num),2]
    sister_num = sibling_num[sibling_num!=node_num]
    return(sister_num)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_ancestor_num = function(phy, node_num) {
    ancestor_num = c()
    root_num = get_root_num(phy)
    current_node_num = node_num
    for (i in seq_len(phy[['Nnode']])) {
        if (current_node_num==root_num) {
            break
        }
        parent_num = get_parent_num(phy, current_node_num)
        ancestor_num = c(ancestor_num, parent_num)
        current_node_num = parent_num
    }
    return(ancestor_num)
}

# copied from rkftools https://github.com/kfuku52/rkftools
read_notung_parsable = function(file, mode='D') {
    cols = c('event', 'gn_node', 'lower_sp_node', 'upper_sp_node')
    if (mode=='D') {
        event_lines = readLines(file, warn=FALSE)
        dup_positions = grep("^\\s*#D\\b", event_lines, perl=TRUE, ignore.case=TRUE)
        if (length(dup_positions)>0) {
            dup_lines = trimws(event_lines[dup_positions])
            dup_items = lapply(dup_lines, function(line) {
                strsplit(line, "\\s+")[[1]]
            })
                malformed_idx = sapply(dup_items, function(items) {
                    (length(items) > 0) && (toupper(items[[1]]) == '#D') && (length(items) < 4)
                })
            if (any(malformed_idx)) {
                stop('Malformed #D line(s) were found in the NOTUNG parsable file.')
            }
                dup_items = Filter(function(items) {
                    (length(items) >= 4) &&
                    (toupper(items[[1]]) == '#D') &&
                    !(
                        (toupper(items[[2]]) == 'DUPLICATION') ||
                        ((toupper(items[[2]]) == 'GENE') && (toupper(items[[3]]) == 'NODE'))
                )
            }, dup_items)
            if (length(dup_items) > 0) {
                normalize_sp_node = function(x) {
                    x2 = trimws(x)
                    if (toupper(x2) %in% c('N/A', 'NA', 'NONE', 'NULL', '-', 'NIL')) {
                        return(NA_character_)
                    }
                    return(x2)
                }
                df = data.frame(
                    event = rep('D', length(dup_items)),
                    gn_node = sapply(dup_items, function(items) items[[2]]),
                    lower_sp_node = sapply(dup_items, function(items) normalize_sp_node(items[[3]])),
                    upper_sp_node = sapply(dup_items, function(items) normalize_sp_node(items[[4]])),
                    stringsAsFactors = FALSE
                )
            } else {
                df = data.frame(matrix(NA, 0, length(cols)))
                colnames(df) = cols
            }
        } else {
            df = data.frame(matrix(NA, 0, length(cols)))
            colnames(df) = cols
        }
    } else {
        cat('mode', mode, 'is not supported.')
        df = data.frame(matrix(NA, 0, length(cols)))
        colnames(df) = cols
    }
    return(df)
}

# copied from rkftools https://github.com/kfuku52/rkftools
pad_short_edges = function(tree, threshold=1e-6, external_only=FALSE) {
    stopifnot(ape::is.binary(tree))
    edge_idx = 1:nrow(tree$edge)
    is_target_edge = TRUE
    if (external_only) {
        is_target_edge = is_target_edge & (tree$edge[,2]<=length(tree$tip.label))
    }
    edge_lengths = tree[['edge.length']][is_target_edge]
    min_eel = min(edge_lengths)
    cat('Minimum edge length:', min_eel, '\n')
    is_short_eel = (is_target_edge)&(tree$edge.length<threshold)
    num_short_eel = sum(is_short_eel)
    cat('Number of short edges ( length <', threshold, '):', num_short_eel, '\n')
    if (num_short_eel>0) {
        short_eel_idx = edge_idx[is_short_eel]
        for (i in short_eel_idx) {
            if (tree$edge.length[i]<threshold) {
                shift_value = threshold - tree$edge.length[i]
                sister_node_num = get_sister_num(tree, tree$edge[i,2])
                sister_edge_idx = edge_idx[tree$edge[,2]==sister_node_num]
                root_num = get_root_num(tree)
                flag = TRUE
                flag_root = FALSE
                current_idx = i
                while (flag==TRUE) {
                    parent_node_num = tree$edge[current_idx,1]
                    parent_edge_idx = edge_idx[tree$edge[,2]==parent_node_num]
                    parent_edge_length = tree$edge.length[parent_edge_idx]
                    if (parent_node_num==root_num) {
                        flag = FALSE
                        flag_root = TRUE
                    } else if (parent_edge_length>=threshold+shift_value) {
                        flag = FALSE
                    } else {
                        current_idx = edge_idx[tree$edge[,2]==parent_node_num]
                    }
                }

                tree$edge.length[i] = tree$edge.length[i] +shift_value
                tree$edge.length[sister_edge_idx] = tree$edge.length[sister_edge_idx] + shift_value
                if (flag_root) {
                    cat('Adding branch length to subroot edges,', i, 'and', sister_edge_idx, '\n')
                } else {
                    cat('Transfering branch length from edge', parent_edge_idx, 'to', i, 'and', sister_edge_idx, '\n')
                    tree$edge.length[parent_edge_idx] = tree$edge.length[parent_edge_idx] - shift_value
                }
            }
        }
    }
    return(tree)
}

# copied from rkftools https://github.com/kfuku52/rkftools
force_ultrametric = function(tree, stop_if_larger_change=0.01) {
    if (ape::is.ultrametric(tree)) {
        cat('The tree is ultrametric.\n')
    } else {
        cat('The tree is not ultrametric. Adjusting the branch length.\n')
        edge_length_before = tree[['edge.length']]
        tree = ape::chronoMPL(tree)
        edge_length_after = tree[['edge.length']]
        sum_adjustment = sum(abs(edge_length_after-edge_length_before))
        cat('Total branch length difference between before- and after-adjustment:', sum_adjustment, '\n')
        stopifnot(sum_adjustment<(sum(tree[['edge.length']]) * stop_if_larger_change))
    }
    return(tree)
}

# copied from rkftools https://github.com/kfuku52/rkftools
contains_polytomy = function(phy) {
    # RADTE expects rooted, fully bifurcated trees.
    # Internal nodes must have exactly two children.
    child_counts = table(phy[['edge']][,1])
    return(any(child_counts != 2))
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_species_names = function(phy, sep='_') {
    split_names = strsplit(phy[['tip.label']], sep)
    species_names = c()
    for (sn in split_names) {
        species_names = c(species_names, paste0(sn[1], sep, sn[2]))
    }
    return(species_names)
}

# copied from rkftools https://github.com/kfuku52/rkftools
leaf2species = function(leaf_names) {
    split = strsplit(leaf_names, '_')
    species_names = c()
    for (i in seq_along(split)) {
        if (length(split[[i]])>=3) {
            species_names = c(
                species_names,
                paste(split[[i]][[1]], split[[i]][[2]])
            )
        } else {
            warning('leaf name could not be interpreted as genus_species_gene: ', split[[i]], '\n')
        }
    }
    return(species_names)
}

validate_species_tree_labels = function(sp_tree) {
    node_labels = sp_tree[['node.label']]
    if (is.null(node_labels)) {
        node_labels = character(0)
    }
    if (length(node_labels) != sp_tree[['Nnode']]) {
        stop('Input species tree contains non-labeled node(s).')
    }
    if (any(is.na(node_labels) | node_labels=='')) {
        stop('Input species tree contains non-labeled node(s).')
    }
    if (any(duplicated(node_labels))) {
        duplicated_labels = unique(node_labels[duplicated(node_labels)])
        stop(
            'Input species tree contains duplicated internal node label(s): ',
            paste(duplicated_labels, collapse=', ')
        )
    }
    all_node_labels = c(sp_tree[['tip.label']], node_labels)
    if (any(duplicated(all_node_labels))) {
        duplicated_labels = unique(all_node_labels[duplicated(all_node_labels)])
        stop(
            'Input species tree contains duplicated tip/internal node label(s): ',
            paste(duplicated_labels, collapse=', ')
        )
    }
    sp_tree[['node.label']] = node_labels
    return(sp_tree)
}

validate_gene_tip_labels = function(tip_labels) {
    split_labels = strsplit(tip_labels, '_')
    is_invalid = sapply(split_labels, function(items) {
        if (length(items) < 3) {
            return(TRUE)
        }
        genus = items[[1]]
        species = items[[2]]
        gene_id = paste(items[3:length(items)], collapse='_')
        (nchar(genus)==0) || (nchar(species)==0) || (nchar(gene_id)==0)
    })
    if (any(is_invalid)) {
        invalid_labels = tip_labels[is_invalid]
        max_show = min(5, length(invalid_labels))
        shown_labels = paste(invalid_labels[1:max_show], collapse=', ')
        extra_suffix = ''
        if (length(invalid_labels) > max_show) {
            extra_suffix = paste0(' ... (', length(invalid_labels)-max_show, ' more)')
        }
        stop(
            'Input gene tree tip label(s) must follow GENUS_SPECIES_GENEID format. Invalid label(s): ',
            shown_labels,
            extra_suffix
        )
    }
}

validate_gene_species_coverage = function(gn_tree, sp_tree) {
    gn_species = unique(get_species_names(gn_tree))
    missing_species = setdiff(gn_species, sp_tree[['tip.label']])
    if (length(missing_species) > 0) {
        stop(
            'Species in the gene tree were not found in the species tree: ',
            paste(missing_species, collapse=', ')
        )
    }
}

validate_gene_tree_labels = function(gn_tree) {
    node_labels = gn_tree[['node.label']]
    if (is.null(node_labels)) {
        node_labels = character(0)
    }
    if (length(node_labels) != gn_tree[['Nnode']]) {
        stop('Input gene tree contains non-labeled internal node(s).')
    }
    if (any(is.na(node_labels) | node_labels=='')) {
        stop('Input gene tree contains non-labeled internal node(s).')
    }
    if (any(duplicated(node_labels))) {
        duplicated_labels = unique(node_labels[duplicated(node_labels)])
        stop(
            'Input gene tree contains duplicated internal node label(s): ',
            paste(duplicated_labels, collapse=', ')
        )
    }
    all_node_labels = c(gn_tree[['tip.label']], node_labels)
    if (any(duplicated(all_node_labels))) {
        duplicated_labels = unique(all_node_labels[duplicated(all_node_labels)])
        stop(
            'Input gene tree contains duplicated tip/internal node label(s): ',
            paste(duplicated_labels, collapse=', ')
        )
    }
    gn_tree[['node.label']] = node_labels
    return(gn_tree)
}

validate_gn_node_table = function(gn_node_table) {
    if (nrow(gn_node_table) == 0) {
        return(invisible(NULL))
    }
    if (any(is.na(gn_node_table[['gn_node']]) | gn_node_table[['gn_node']]=='')) {
        stop('Gene node table contains missing gn_node value(s).')
    }
    if (any(is.na(gn_node_table[['gn_node_num']]))) {
        stop('Gene node table contains missing gn_node_num value(s).')
    }
    if (any(duplicated(gn_node_table[['gn_node']]))) {
        duplicated_names = unique(gn_node_table[['gn_node']][duplicated(gn_node_table[['gn_node']])])
        stop('Gene node table contains duplicated gn_node value(s): ', paste(duplicated_names, collapse=', '))
    }
    if (any(duplicated(gn_node_table[['gn_node_num']]))) {
        duplicated_nums = unique(gn_node_table[['gn_node_num']][duplicated(gn_node_table[['gn_node_num']])])
        stop('Gene node table contains duplicated gn_node_num value(s): ', paste(duplicated_nums, collapse=', '))
    }
}

validate_duplication_nodes_internal = function(gn_node_table, gn_tree) {
    if (nrow(gn_node_table) == 0) {
        return(invisible(NULL))
    }
    is_dup = grepl('^D', gn_node_table[['event']])
    if (!any(is_dup)) {
        return(invisible(NULL))
    }
    dup_rows = gn_node_table[is_dup, , drop=FALSE]
    is_tip_dup = dup_rows[['gn_node_num']] <= ape::Ntip(gn_tree)
    if (any(is_tip_dup)) {
        invalid_nodes = unique(dup_rows[is_tip_dup, 'gn_node'])
        stop(
            'Gene node table contains duplication event(s) mapped to tip node(s): ',
            paste(invalid_nodes, collapse=', '),
            '. Duplication annotations must refer to internal gene-tree nodes.'
        )
    }
}

validate_gn_species_nodes = function(gn_node_table, sp_tree) {
    valid_sp_nodes = c(sp_tree[['tip.label']], sp_tree[['node.label']])
    lower_nodes = unique(gn_node_table[['lower_sp_node']])
    lower_nodes = lower_nodes[!is.na(lower_nodes) & (lower_nodes != '')]
    missing_lower = setdiff(lower_nodes, valid_sp_nodes)
    if (length(missing_lower) > 0) {
        stop(
            'Gene node table contains lower_sp_node value(s) not found in the species tree: ',
            paste(missing_lower, collapse=', ')
        )
    }
    upper_nodes = unique(gn_node_table[['upper_sp_node']])
    upper_nodes = upper_nodes[!is.na(upper_nodes) & (upper_nodes != '')]
    missing_upper = setdiff(upper_nodes, valid_sp_nodes)
    if (length(missing_upper) > 0) {
        stop(
            'Gene node table contains upper_sp_node value(s) not found in the species tree: ',
            paste(missing_upper, collapse=', ')
        )
    }
}

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

# copied from rkftools https://github.com/kfuku52/rkftools
transfer_node_labels = function(phy_from, phy_to) {
    for (t in seq_along(phy_to$node.label)) {
        to_node = phy_to$node.label[t]
        to_clade = extract.clade(phy=phy_to, node=to_node, root.edge = 0, interactive = FALSE)
        to_leaves = to_clade$tip.label
        for (f in seq_along(phy_from$node.label)) {
            from_node = phy_from$node.label[f]
            from_clade = extract.clade(phy=phy_from, node=from_node, root.edge = 0, interactive = FALSE)
            from_leaves = from_clade$tip.label
            if (setequal(to_leaves, from_leaves)) {
                phy_to$node.label[t] = from_node
                break
            }
        }
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
    for (gn_node_num in gn_node_table[['gn_node_num']]) {
        if (gn_node_num==root_num) {
            next
        }
        ancestor_nums = get_ancestor_num(gn_tree, gn_node_num)
        ancestor_rows = gn_node_table[['gn_node_num']] %in% ancestor_nums
        if (!any(ancestor_rows)) {
            next
        }
        child_lower = gn_node_table[(gn_node_table[['gn_node_num']]==gn_node_num),'lower_age']
        child_upper = gn_node_table[(gn_node_table[['gn_node_num']]==gn_node_num),'upper_age']
        ancestor_lower = min(gn_node_table[ancestor_rows,'lower_age'])
        ancestor_upper = min(gn_node_table[ancestor_rows,'upper_age'])
        is_conflict = isTRUE((child_lower>=ancestor_lower) & (child_upper>=ancestor_upper))
        if (is_conflict) {
            conflicts = rbind(
                conflicts,
                data.frame(
                    node=gn_node_num,
                    child_lower=child_lower,
                    child_upper=child_upper,
                    ancestor_lower=ancestor_lower,
                    ancestor_upper=ancestor_upper,
                    stringsAsFactors=FALSE
                )
            )
        }
    }
    return(conflicts)
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
    node_nums = gn_node_table2[['gn_node_num']]
    ancestor_count = sapply(node_nums, function(node_num) {
        length(get_ancestor_num(gn_tree, node_num))
    })
    process_order = node_nums[order(ancestor_count, decreasing=FALSE)]

    for (gn_node_num in process_order) {
        if (gn_node_num==root_num) {
            next
        }
        row_idx = which(gn_node_table2[['gn_node_num']]==gn_node_num)
        ancestor_nums = get_ancestor_num(gn_tree, gn_node_num)
        ancestor_rows = gn_node_table2[['gn_node_num']] %in% ancestor_nums
        if (!any(ancestor_rows)) {
            next
        }
        ancestor_lower = min(gn_node_table2[ancestor_rows,'lower_age'])
        ancestor_upper = min(gn_node_table2[ancestor_rows,'upper_age'])
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
            adjusted_nodes = rbind(
                adjusted_nodes,
                data.frame(
                    node=gn_node_num,
                    lower_before=lower_before,
                    upper_before=upper_before,
                    lower_after=lower_after,
                    upper_after=upper_after,
                    ancestor_upper=ancestor_upper,
                    stringsAsFactors=FALSE
                )
            )
        }
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
    node_nums = calibration_table2$node
    ancestor_count = sapply(node_nums, function(node_num) {
        length(get_ancestor_num(phy, node_num))
    })
    process_order = node_nums[order(ancestor_count, decreasing=FALSE)]
    adjusted_nodes = c()
    for (node_i in process_order) {
        if (node_i==root_num) {
            next
        }
        row_idx = which(calibration_table2$node==node_i)
        ancestor_nodes = intersect(get_ancestor_num(phy, node_i), calibration_table2$node)
        if (length(ancestor_nodes)==0) {
            next
        }
        ancestor_rows = calibration_table2$node %in% ancestor_nodes
        ancestor_upper = min(calibration_table2$age.max[ancestor_rows])
        ancestor_lower = min(calibration_table2$age.min[ancestor_rows])

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
                chronos(
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
                chronos(
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
        chronos(
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
    depth_values = node.depth.edgelength(phy)
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
    ancestor_count = sapply(node_nums, function(node_num) {
        length(get_ancestor_num(phy, node_num))
    })
    process_order_up = node_nums[order(ancestor_count, decreasing=FALSE)]
    process_order_down = rev(process_order_up)
    max_iter = node_count * 4
    for (iter_i in seq_len(max_iter)) {
        changed = FALSE
        for (child_node in process_order_down) {
            if (child_node==root_num) {
                next
            }
            parent_node = get_parent_num(phy, child_node)
            if (length(parent_node)==0) {
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

save_tree_pdf = function(phy, file, show.age=FALSE, edge_colors=list()) {
    phy = ape::ladderize(phy)
    if (show.age) {
        root_depth = max(node.depth.edgelength(phy))
        node_ages = abs(node.depth.edgelength(phy) - root_depth)
        int_node_ages = node_ages[(length(phy$tip.label)+1):length(node_ages)]
        phy$node.label = paste(phy$node.label, as.character(round(int_node_ages, digits=1)))
    }
    ec2 = rep('black', nrow(phy[['edge']]))
    node_colors = 'black'
    if (length(edge_colors)!=0) {
        for (col in names(edge_colors)) {
            ec2[(phy[['edge']][,2]%in%edge_colors[[col]])] = col
        }
    }
    if (length(edge_colors)!=0) { # Should not be merged to the previous if block
        is_node = (phy[['edge']][,2]>length(phy[['tip.label']]))
        node_order = order(phy[['edge']][,2][is_node])
        node_colors = ec2[is_node][node_order]
        root_num = get_root_num(phy)
        for (col in names(edge_colors)) {
            if (root_num %in% edge_colors[[col]]) {
                node_colors = c(col, node_colors) # Adding root
                break
            }
        }
    }
    pdf(file, height=max(3, length(phy$tip.label)/5+1), width=7.2)
    plot(phy, show.node.label=FALSE, show.tip.label=TRUE, cex=0.5, label.offset=0, 
         edge.color=ec2, root.edge=TRUE)
    nodelabels(text=phy[['node.label']], col=node_colors, bg='white', cex=0.5)
    invisible(dev.off())
}


cat('arguments:\n')
args = commandArgs(trailingOnly=TRUE)
args = get_parsed_args(args, print=TRUE)

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

required_args = c('species_tree', 'max_age', 'chronos_lambda', 'chronos_model')
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
chronos_attempt_timeout_sec = Inf
chronos_total_timeout_sec = Inf
if (!allow_constraint_drop) {
    # In no-drop mode, avoid hanging forever and fall back deterministically.
    chronos_attempt_timeout_sec = 60
    chronos_total_timeout_sec = 300
}
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

cat('\nStart: species tree processing', '\n')
tree_text0 = scan(sp_file, what=character(), sep="\n", blank.lines.skip=FALSE)
tree_text1 = gsub("'([0-9]+)'", "PLACEHOLDER\\1", tree_text0)
sp_tree = read.tree(text=tree_text1)
if (!is.null(sp_tree[['node.label']])) {
    sp_tree[['node.label']] = sub('^PLACEHOLDER', '', sp_tree[['node.label']])
}
sp_tree = validate_species_tree_labels(sp_tree)
if (contains_polytomy(sp_tree)) {
    stop('Input species tree contains polytomy. A completely bifurcated tree is expected as input.')
}
validate_tree_edge_lengths(sp_tree, 'species tree')
if (length(args[['pad_short_edge']])) {
    sp_tree = pad_short_edges(sp_tree, threshold=args[['pad_short_edge']], external_only=FALSE)
}
sp_tree = force_ultrametric(sp_tree, stop_if_larger_change=0.01)
root_depth = max(node.depth.edgelength(sp_tree))
sp_node_ages = abs(node.depth.edgelength(sp_tree) - root_depth)
sp_node_names = c(sp_tree[['tip.label']], sp_tree[['node.label']])
sp_node_table = data.frame(node=sp_node_names, age=sp_node_ages, spp=NA, stringsAsFactors=FALSE)
for (sp_sub in ape::subtrees(sp_tree)) {
    subroot_node = sp_sub[['node.label']][1]
    sp_node_table[(sp_node_table$node==subroot_node),'spp'] = paste(sp_sub[['tip.label']], collapse=',')
}
max_tip_age = max(sp_node_table[is.na(sp_node_table[['spp']]),'age'])
if (max_tip_age!=0) {
    cat(paste0('Nonzero tip age(s) were detected (max=', max_tip_age, '). Coercing to 0.\n'))
    sp_node_table[is.na(sp_node_table[['spp']]),'age'] = 0
}
cat('End: species tree processing', '\n\n')

cat('Start: gene tree processing', '\n')

read_generax_nhx = function(generax_file) {
    treetext = readLines(generax_file, warn=FALSE)
    count_matches = function(pattern, text_vec) {
        matches = gregexpr(pattern, text_vec, perl=TRUE)
        return(sum(vapply(matches, function(m) {
            if ((length(m)==1) && (m[[1]]==-1)) {
                return(0L)
            }
            return(as.integer(length(m)))
        }, integer(1))))
    }
    num_open = count_matches('\\(', treetext)
    num_close = count_matches('\\)', treetext)
    if ((num_open - num_close)==-1) {
        cat('Number of parentheses in the .nhx is not consistent. Trying to fix.')
        treetext <- gsub("\\)\\s*;", ";", treetext, perl=TRUE)
    }
    tmp_nhx_file = tempfile(pattern='radte_nhx_', fileext='.txt')
    write(treetext, tmp_nhx_file)
    on.exit(file.remove(tmp_nhx_file), add=TRUE)
    nhxtree = treeio::read.nhx(tmp_nhx_file)
    return(nhxtree)
}

if (mode=='generax') {
    cat('Reading GeneRax tree.\n')
    nhxtree = read_generax_nhx(generax_file)

    gn_tree = nhxtree@phylo
    if (contains_polytomy(gn_tree)) {
        stop('Input tree contains polytomy. A completely bifurcated tree is expected as input.')
    }
    gn_tree = validate_gene_tree_labels(gn_tree)
    validate_gene_tip_labels(gn_tree[['tip.label']])
    validate_gene_species_coverage(gn_tree, sp_tree)
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
    gn_node_table[(gn_node_table[['event']]=='D'),'upper_sp_node'] = NA
    for (sp_node in unique(gn_node_table[['lower_sp_node']])) {
        node_num = get_node_num_by_name(sp_tree, sp_node)
        parent_num = get_parent_num(sp_tree, node_num)
        parent_name = get_node_name_by_num(sp_tree, parent_num)
        if (identical(parent_name, character(0))) {
            parent_name = NA
        }
        conditions = (gn_node_table[['lower_sp_node']]==sp_node)
        conditions = conditions & (gn_node_table[['event']]=='D')
        gn_node_table[conditions,'upper_sp_node'] = parent_name
    }
    gn_node_table[,'lower_age'] = NA
    gn_node_table[,'upper_age'] = NA
    for (sp_node in sp_node_table[['node']]) {
        node_age = as.numeric(sp_node_table[(sp_node_table[['node']]==sp_node),'age'])
        conditions = (gn_node_table[['lower_sp_node']]==sp_node)
        conditions = conditions & (gn_node_table[['upper_sp_node']]==sp_node)
        conditions[is.na(conditions)] = FALSE
        gn_node_table[conditions,'lower_age'] = node_age
        gn_node_table[conditions,'upper_age'] = node_age
    }
    gn_node_table[,'gn_node_num'] = get_node_num_by_name(gn_tree, gn_node_table[['gn_node']])
    gn_node_table = data.frame(gn_node_table[,cols], stringsAsFactors=FALSE)
    validate_gn_species_nodes(gn_node_table, sp_tree)
} else if (mode=='notung') {
    cat('Reading NOTUNG tree.\n')
    gn_tree = read.tree(gn_file)
    gn_tree[['node.label']] = gsub("\\'", "",gn_tree[['node.label']])
    if (contains_polytomy(gn_tree)) {
        stop('Input tree contains polytomy. A completely bifurcated tree is expected as input.')
    }
    gn_tree = validate_gene_tree_labels(gn_tree)
    validate_gene_tip_labels(gn_tree[['tip.label']])
    validate_gene_species_coverage(gn_tree, sp_tree)
    gn_tree = pad_branch_length(gn_tree, pad_size=args[['pad_short_edge']])

    gn_node_table = read_notung_parsable(file=parsable_file, mode='D')
    gn_node_table = merge(gn_node_table, data.frame(lower_age=NA, upper_age=NA, spp=NA), all=TRUE)
    check_gn_node_name_uniqueness(gn_node_table, gn_tree)
    if (nrow(gn_node_table) > 0) {
        gn_node_nums = sapply(gn_node_table[,'gn_node'], function(x){get_node_num_by_name(gn_tree, x)})
        gn_node_table$gn_node_num = gn_node_nums
        for (i in 1:nrow(gn_node_table)) {
            lower_sp_node_i = gn_node_table$lower_sp_node[i]
            if ((!is.na(lower_sp_node_i)) && any(sp_node_table$node==lower_sp_node_i)) {
                gn_node_table$lower_age[i] = sp_node_table$age[sp_node_table$node==lower_sp_node_i]
            }
            upper_sp_node_i = gn_node_table$upper_sp_node[i]
            if ((!is.na(upper_sp_node_i)) && any(sp_node_table$node==upper_sp_node_i)) {
                gn_node_table$upper_age[i] = sp_node_table$age[sp_node_table$node==upper_sp_node_i]
            }
        }
    } else {
        gn_node_table = gn_node_table[0,]
    }
    
    for (gn_sub in ape::subtrees(gn_tree)) {
        root_node = gn_sub$node.label[1]
        if (! root_node %in% gn_node_table$gn_node) {
            root_num = get_node_num_by_name(gn_tree, root_node)
            node_spp = unique(leaf2species(gn_sub[['tip.label']]))
            node_spp = sub(' ', '_', node_spp)
            sp_node_sets = strsplit(sp_node_table$spp, ',', fixed=TRUE)
            is_spnode_species = sapply(sp_node_sets, function(spp_set) {
                all(node_spp %in% spp_set)
            })
            if (!any(is_spnode_species)) {
                stop(
                    'Could not map gene subtree species to a species-tree clade: ',
                    paste(node_spp, collapse=', ')
                )
            }
            node_age = min(sp_node_table[is_spnode_species,'age'])
            is_min = (sp_node_table[,'age']==node_age)
            sp_node = sp_node_table[is_min&is_spnode_species,'node']
            if (length(sp_node) != 1) {
                stop(
                    'Ambiguous species-tree clade mapping for gene subtree species: ',
                    paste(node_spp, collapse=', ')
                )
            }
            ind = nrow(gn_node_table)+1
            gn_node_table[ind,'event'] = "S"
            gn_node_table[ind,'gn_node'] = root_node
            gn_node_table[ind,'gn_node_num'] = root_num
            gn_node_table[ind,'lower_sp_node'] = sp_node
            gn_node_table[ind,'upper_sp_node'] = sp_node
            gn_node_table[ind,'lower_age'] = node_age
            gn_node_table[ind,'upper_age'] = node_age
            gn_node_table[ind,'spp'] = paste(node_spp, collapse='|')
        }
    }
}

cat('End: gene tree processing', '\n\n')
validate_gn_node_table(gn_node_table)
validate_duplication_nodes_internal(gn_node_table, gn_tree)

# Calibration node check
if ((sum(gn_node_table[['event']]=="D") > 0)&(any(is.na(gn_node_table[['upper_age']])))) {
    gn_spp = unique(get_species_names(gn_tree))
    num_sp = length(gn_spp)
    cat('# species in the gene tree:', num_sp, '\n')
    cat('Species in the gene tree:', paste(gn_spp, collapse=', '), '\n')
    num_sp_gntree = max(1, ape::getMRCA(sp_tree, gn_spp))
    if (num_sp_gntree==get_root_num(sp_tree)) {
        divtime_max = max_age
        divtime_min = max(ape::node.depth.edgelength(sp_tree))
    } else {
        if (length(gn_spp)==1) {
            num_mrca = get_node_num_by_name(sp_tree, gn_spp)
        } else {
            num_mrca = ape::getMRCA(sp_tree, gn_spp)
        }
        num_parent = sp_tree$edge[,1][sp_tree$edge[,2]==num_mrca]
        label_mrca = get_node_name_by_num(phy=sp_tree, node_num=num_mrca)
        label_parent = get_node_name_by_num(phy=sp_tree, node_num=num_parent)
        divtime_max = sp_node_table[sp_node_table$node==label_parent,'age']
        divtime_min = sp_node_table[sp_node_table$node==label_mrca,'age']
        cat('Species in the MRCA species tree clade:', paste(sp_node_table[sp_node_table$node==label_mrca,'spp'], collapse=', '), '\n')
        cat('Species in the parent species tree clade:', paste(sp_node_table[sp_node_table$node==label_parent,'spp'], collapse=', '), '\n')
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
cat('\n')
gn_node_table_calibration = gn_node_table[(gn_node_table[,'gn_node_num']>ape::Ntip(gn_tree)),]
num_constrained_speciation = sum(grepl('^S', gn_node_table_calibration[,'event']))
num_constrained_duplication = sum(grepl('^D', gn_node_table_calibration[,'event']))
cat('Number of constrained speciation nodes:', num_constrained_speciation, '\n')
cat('Number of constrained duplication nodes:', num_constrained_duplication, '\n')

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

# chronos
chronos_out = NULL
chronos_control = chronos.control()
chronos_control$iter.max = 100000
chronos_control$eval.max = 100000
chronos_control$dual.iter.max = 200

has_duplication_event = any(grepl('^D', gn_node_table$event))
if (!has_duplication_event) {
    # Gene tree without duplication nodes
    calibrated_node = "allS"
    cat("Constrained nodes:", calibrated_node, '\n')
    cat("All nodes are speciation nodes. Transferring node ages from species tree without age inference by chronos.", '\n')
    dup_constraint = NA
    gn_spp = unique(get_species_names(gn_tree))
    drop_spp = sp_tree$tip.label[! sp_tree$tip.label %in% gn_spp]
    if (length(drop_spp) > 0) {
        chronos_out = drop.tip(phy=sp_tree, tip=drop_spp, trim.internal = TRUE)
    } else {
        chronos_out = sp_tree
    }
    gn_tip_index = c()
    for (sp in chronos_out$tip.label) {
        tip_matches = which(startsWith(gn_tree$tip.label, paste0(sp, '_')) | (gn_tree$tip.label==sp))
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
    # Gene tree with duplication nodes
    validate_tree_edge_lengths(gn_tree, 'gene tree for chronos')
    # Normalize edge lengths to prevent numerical overflow in chronos
    gn_tree = normalize_edge_length_range(gn_tree, max_edge = 1000, min_edge = 1e-8)
    chronos_out = structure('PLACEHOLDER', class='try-error')
    calibrated_node = 'RS'
    current_calibration_table = calibration_tables[['RS']]

    chronos_model_used = chronos_model
    chronos_lambda_used = chronos_lambda
    chronos_seed_used = NA_integer_
    seed_cursor = 1
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

if ("try-error" %in% class(chronos_out)) {
    stop('All attempts for divergence time estimation were failed. Exiting.')
} else {
    cat('Writing output files.\n')
    chronos_out2 = chronos_out
    num_neg = 1
    counter = 1
    if (length(args[['pad_short_edge']])) {
        while ((num_neg>0)&(counter<100)) {
            cat(paste0('Branch length padding round ', counter, ' started.\n'))
            chronos_out2 = pad_short_edges(chronos_out2, threshold=args[['pad_short_edge']], external_only=FALSE)
            chronos_out2 = force_ultrametric(chronos_out2, stop_if_larger_change=0.01)
            num_neg = sum(chronos_out2[['edge.length']]<0)
            cat(num_neg, 'negative value(s) were detected in estimated branch length.\n\n')
            counter = counter + 1
        }
    }

    write(calibrated_node, file='radte_calibrated_nodes.txt')

    write.tree(chronos_out2, file="radte_gene_tree_output.nwk")
    current_calibration_table = merge(current_calibration_table, gn_node_table[,c('gn_node_num','event')], by.x='node', by.y='gn_node_num', all.x=TRUE)
    current_calibration_table[current_calibration_table$node==calibration_table_R$node,'event'] = 'R'
    write.table(current_calibration_table, file='radte_calibration_used.tsv', sep='\t', quote=FALSE, row.names=FALSE)

    calibration_table = merge(calibration_table, gn_node_table[,c('gn_node_num','event')], by.x='node', by.y='gn_node_num', all.x=TRUE)
    write.table(calibration_table, file='radte_calibration_all.tsv', sep='\t', quote=FALSE, row.names=FALSE)

    gn_node_table$spp = NULL
    write.table(gn_node_table, file='radte_gene_tree.tsv', sep='\t', quote=FALSE, row.names=FALSE)

    sp_node_table$spp = NULL
    write.table(sp_node_table, file='radte_species_tree.tsv', sep='\t', quote=FALSE, row.names=FALSE)
    
    node_nums = (length(chronos_out2[['tip.label']])+1):max(chronos_out2[['edge']])
    noncalibrated_nodes = node_nums[!node_nums %in% current_calibration_table[['node']]]
    ec = list('red'=noncalibrated_nodes, 'blue'=current_calibration_table[['node']])
	    save_tree_pdf(phy=gn_tree, file="radte_gene_tree_input.pdf", show.age=FALSE, edge_colors=ec)
	    save_tree_pdf(phy=chronos_out2, file="radte_gene_tree_output.pdf", show.age=TRUE, edge_colors=ec)
	    save_tree_pdf(phy=sp_tree, file="radte_species_tree.pdf", show.age=TRUE)

	    if (exists('chronos_model_used')) {
	        cat('chronos model used:', chronos_model_used, '\n')
	    }
	    if (exists('chronos_lambda_used')) {
	        cat('chronos lambda used:', chronos_lambda_used, '\n')
	    }
	    if (exists('chronos_seed_used') && (!is.na(chronos_seed_used))) {
	        cat('chronos seed used:', chronos_seed_used, '\n')
	    }
	    cat('Calibrated nodes:', calibrated_node, '\n')
	    cat('Tree height:', max(ape::node.depth.edgelength(sp_tree)), '\n')
    is_max_age = (calibration_table[,'age.max']==max_age)
    num_spnode_used_for_constraint = nrow(unique(calibration_table[!is_max_age,c('age.min','age.max')]))
    cat('Number of species tree node used for the gene tree constraint:', num_spnode_used_for_constraint, '\n')    
    cat('Completed: RADTE divergence time estimation', '\n')
}
