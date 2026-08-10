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
    is_valid = (node_num >= 1L) & (node_num <= length(node_names))
    node_num = node_num[is_valid]
    if (length(node_num)==0) {
        return(character(0))
    }
    return(unname(node_names[node_num]))
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
    if (anyDuplicated(node_names)==0) {
        node_num = match(node_name, node_names)
        node_num = node_num[!is.na(node_num)]
    } else {
        # Preserve the historical behavior for invalid trees with duplicated
        # labels: return every matching node so validation can report it.
        node_nums = seq_along(node_names)
        node_num = unlist(
            lapply(node_name, function(x) node_nums[node_names == x]),
            use.names=FALSE
        )
    }
    if (length(node_num)==0) {
        return(integer(0))
    }
    return(as.integer(node_num))
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_root_num = function(phy) {
    root_num = setdiff(phy[['edge']][,1], phy[['edge']][,2])
    return(root_num)
}

get_parent_map = function(phy) {
    num_nodes = length(phy[['tip.label']]) + phy[['Nnode']]
    parent_map = rep(NA_integer_, num_nodes)
    parent_map[phy[['edge']][,2]] = phy[['edge']][,1]
    return(parent_map)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_parent_num = function(phy, node_num) {
    if (length(node_num)==0) {
        return(integer(0))
    }
    node_num = suppressWarnings(as.integer(node_num))
    node_num = node_num[!is.na(node_num)]
    parent_map = get_parent_map(phy)
    is_valid = (node_num >= 1L) & (node_num <= length(parent_map))
    parent_num = parent_map[node_num[is_valid]]
    return(as.integer(parent_num[!is.na(parent_num)]))
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_sister_num = function(phy, node_num) {
    parent_num = get_parent_num(phy, node_num)
    if (length(parent_num)==0) {
        return(integer(0))
    }
    sibling_num = phy[['edge']][(phy[['edge']][,1]==parent_num),2]
    sister_num = sibling_num[sibling_num!=node_num]
    return(as.integer(sister_num))
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_ancestor_num = function(phy, node_num, parent_map=NULL, root_num=NULL) {
    ancestor_num = integer(0)
    if (length(node_num)!=1 || is.na(node_num)) {
        return(ancestor_num)
    }
    if (is.null(parent_map)) {
        parent_map = get_parent_map(phy)
    }
    if (is.null(root_num)) {
        root_num = get_root_num(phy)
    }
    current_node_num = node_num
    for (i in seq_len(phy[['Nnode']])) {
        if (current_node_num==root_num) {
            break
        }
        parent_num = parent_map[[current_node_num]]
        if (is.na(parent_num)) {
            break
        }
        ancestor_num = c(ancestor_num, parent_num)
        current_node_num = parent_num
    }
    return(as.integer(ancestor_num))
}

get_node_depth_map = function(phy) {
    num_nodes = length(phy[['tip.label']]) + phy[['Nnode']]
    depth_map = integer(num_nodes)
    edge_preorder = ape::reorder.phylo(phy, order='cladewise')[['edge']]
    for (i in seq_len(nrow(edge_preorder))) {
        parent_num = edge_preorder[i,1]
        child_num = edge_preorder[i,2]
        depth_map[[child_num]] = depth_map[[parent_num]] + 1L
    }
    return(depth_map)
}

get_lca_num = function(node_a, node_b, parent_map, depth_map) {
    if (is.na(node_a) || is.na(node_b)) {
        return(NA_integer_)
    }
    node_a = as.integer(node_a)
    node_b = as.integer(node_b)
    while (depth_map[[node_a]] > depth_map[[node_b]]) {
        node_a = parent_map[[node_a]]
    }
    while (depth_map[[node_b]] > depth_map[[node_a]]) {
        node_b = parent_map[[node_b]]
    }
    while (node_a != node_b) {
        node_a = parent_map[[node_a]]
        node_b = parent_map[[node_b]]
    }
    return(as.integer(node_a))
}

encode_clade_key = function(tip_values) {
    tip_values = sort(unique(as.character(tip_values)))
    encoded = paste0(nchar(tip_values, type='bytes'), ':', tip_values)
    return(paste(encoded, collapse=''))
}

get_internal_clade_keys = function(phy, tip_values=NULL) {
    num_tips = length(phy[['tip.label']])
    num_nodes = num_tips + phy[['Nnode']]
    if (is.null(tip_values)) {
        tip_values = phy[['tip.label']]
    }
    if (length(tip_values) != num_tips) {
        stop('tip_values length does not match the number of tips in the tree.')
    }
    descendant_values = vector('list', num_nodes)
    descendant_values[seq_len(num_tips)] = as.list(as.character(tip_values))
    edge_postorder = ape::reorder.phylo(phy, order='postorder')[['edge']]
    for (i in seq_len(nrow(edge_postorder))) {
        parent_num = edge_postorder[i,1]
        child_num = edge_postorder[i,2]
        descendant_values[[parent_num]] = c(
            descendant_values[[parent_num]],
            descendant_values[[child_num]]
        )
    }
    internal_nodes = seq.int(num_tips + 1L, num_nodes)
    keys = vapply(
        descendant_values[internal_nodes],
        encode_clade_key,
        character(1)
    )
    names(keys) = as.character(internal_nodes)
    return(keys)
}

get_canonical_subtree_ids = function(phy, leaf_ids, interner) {
    num_tips = length(phy[['tip.label']])
    num_nodes = num_tips + phy[['Nnode']]
    if (length(leaf_ids) != num_tips || any(is.na(leaf_ids))) {
        stop('Every tree tip must have one canonical leaf ID.')
    }
    node_ids = rep(NA_integer_, num_nodes)
    node_ids[seq_len(num_tips)] = as.integer(leaf_ids)
    child_ids = vector('list', num_nodes)
    expected_children = tabulate(phy[['edge']][,1], nbins=num_nodes)
    edge_postorder = ape::reorder.phylo(phy, order='postorder')[['edge']]
    for (i in seq_len(nrow(edge_postorder))) {
        parent_num = edge_postorder[i,1]
        child_num = edge_postorder[i,2]
        child_ids[[parent_num]] = c(child_ids[[parent_num]], node_ids[[child_num]])
        if (length(child_ids[[parent_num]]) == expected_children[[parent_num]]) {
            tuple_key = paste(sort(child_ids[[parent_num]]), collapse=',')
            if (!exists(tuple_key, envir=interner$ids, inherits=FALSE)) {
                interner$next_id = interner$next_id + 1L
                assign(tuple_key, interner$next_id, envir=interner$ids)
            }
            node_ids[[parent_num]] = get(tuple_key, envir=interner$ids, inherits=FALSE)
        }
    }
    internal_nodes = seq.int(num_tips + 1L, num_nodes)
    out = node_ids[internal_nodes]
    names(out) = as.character(internal_nodes)
    out
}

get_ancestor_constraint_minima = function(phy, node_nums, values) {
    total_nodes = length(phy[['tip.label']]) + phy[['Nnode']]
    if (length(node_nums) != length(values)) {
        stop('node_nums and values must have the same length.')
    }
    row_by_node = integer(total_nodes)
    row_by_node[as.integer(node_nums)] = seq_along(node_nums)
    minimum = rep(Inf, total_nodes)
    source = rep(NA_integer_, total_nodes)
    edge_preorder = ape::reorder.phylo(phy, order='cladewise')[['edge']]
    for (i in seq_len(nrow(edge_preorder))) {
        parent_num = edge_preorder[i,1]
        child_num = edge_preorder[i,2]
        minimum[[child_num]] = minimum[[parent_num]]
        source[[child_num]] = source[[parent_num]]
        parent_row = row_by_node[[parent_num]]
        if ((parent_row > 0L) && (values[[parent_row]] < minimum[[child_num]])) {
            minimum[[child_num]] = values[[parent_row]]
            source[[child_num]] = parent_num
        }
    }
    list(minimum=minimum, source=source)
}

match_internal_nodes_by_clade = function(phy_from, phy_to) {
    all_tips = sort(unique(c(phy_from[['tip.label']], phy_to[['tip.label']])))
    leaf_id_map = stats::setNames(seq_along(all_tips), all_tips)
    interner = new.env(parent=emptyenv())
    interner$ids = new.env(hash=TRUE, parent=emptyenv())
    interner$next_id = as.integer(length(all_tips))
    from_ids = get_canonical_subtree_ids(
        phy_from,
        unname(leaf_id_map[phy_from[['tip.label']]]),
        interner
    )
    to_ids = get_canonical_subtree_ids(
        phy_to,
        unname(leaf_id_map[phy_to[['tip.label']]]),
        interner
    )
    match(unname(to_ids), unname(from_ids))
}

get_descendant_tip_labels = function(phy, node_num) {
    num_tips = length(phy[['tip.label']])
    if (node_num <= num_tips) {
        return(phy[['tip.label']][[node_num]])
    }
    clade = ape::extract.clade(
        phy=phy,
        node=node_num,
        root.edge=0,
        interactive=FALSE
    )
    return(clade[['tip.label']])
}

map_gene_nodes_to_species_nodes = function(
    gn_tree,
    sp_tree,
    species_parser,
    species_tree_labels=NULL
) {
    if (is.null(species_tree_labels)) {
        species_tree_labels = sp_tree[['tip.label']]
    }
    gene_tip_species = species_parser_get_gene_species(
        species_parser=species_parser,
        tip_labels=gn_tree[['tip.label']],
        species_tree_labels=species_tree_labels,
        strict=TRUE
    )
    species_tip_labels = resolve_species_tree_tips(
        species_parser,
        gene_tip_species,
        species_tree_labels
    )
    if (any(is.na(species_tip_labels))) {
        stop(
            'Species tree tip mapping failed for gene-tree species: ',
            paste(unique(gene_tip_species[is.na(species_tip_labels)]), collapse=', ')
        )
    }
    mapped_nodes = rep(NA_integer_, length(gn_tree[['tip.label']]) + gn_tree[['Nnode']])
    mapped_nodes[seq_along(gn_tree[['tip.label']])] = get_node_num_by_name(
        sp_tree,
        species_tip_labels
    )
    sp_parent_map = get_parent_map(sp_tree)
    sp_depth_map = get_node_depth_map(sp_tree)
    gene_edge_postorder = ape::reorder.phylo(gn_tree, order='postorder')[['edge']]
    for (i in seq_len(nrow(gene_edge_postorder))) {
        parent_num = gene_edge_postorder[i,1]
        child_num = gene_edge_postorder[i,2]
        child_species_node = mapped_nodes[[child_num]]
        if (is.na(mapped_nodes[[parent_num]])) {
            mapped_nodes[[parent_num]] = child_species_node
        } else {
            mapped_nodes[[parent_num]] = get_lca_num(
                mapped_nodes[[parent_num]],
                child_species_node,
                parent_map=sp_parent_map,
                depth_map=sp_depth_map
            )
        }
    }
    return(mapped_nodes)
}

# copied from rkftools https://github.com/kfuku52/rkftools
