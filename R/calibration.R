# Derived feasible ranges are separate from the immutable input calibrations.
# A postorder lower-bound pass and a preorder upper-bound pass suffice on a tree.
calibration_feasible_ranges = function(phy, calibration, min_edge=NULL) {
    total = length(phy$tip.label) + phy$Nnode
    required = c('node', 'age.min', 'age.max')
    if (!all(required %in% names(calibration))) stop('Calibration table is missing node/age.min/age.max columns.')
    nodes = calibration$node
    if (anyNA(nodes) || any(!is.finite(nodes)) || any(nodes != as.integer(nodes)) ||
        any(nodes < 1L | nodes > total) || anyDuplicated(nodes)) {
        stop('Calibration table contains invalid or duplicated node numbers.')
    }
    if (any(!is.finite(calibration$age.min)) || any(!is.finite(calibration$age.max)) ||
        any(calibration$age.min < 0) || any(calibration$age.max < calibration$age.min)) {
        stop('Calibration table contains invalid age bounds.')
    }
    scale = max(c(1, calibration$age.max))
    if (is.null(min_edge)) min_edge = scale * .Machine$double.eps * 16
    if (length(min_edge) != 1L || !is.finite(min_edge) || min_edge <= 0) stop('min_edge must be positive and finite.')
    lower = numeric(total)
    upper = rep(Inf, total)
    upper[seq_along(phy$tip.label)] = 0
    lower[nodes] = calibration$age.min
    upper[nodes] = pmin(upper[nodes], calibration$age.max)
    postorder = ape::reorder.phylo(phy, 'postorder')$edge
    for (i in seq_len(nrow(postorder))) {
        parent = postorder[i, 1]
        child = postorder[i, 2]
        lower[[parent]] = max(lower[[parent]], lower[[child]] + min_edge)
    }
    preorder = ape::reorder.phylo(phy, 'cladewise')$edge
    for (i in seq_len(nrow(preorder))) {
        parent = preorder[i, 1]
        child = preorder[i, 2]
        upper[[child]] = min(upper[[child]], upper[[parent]] - min_edge)
    }
    bad = which(lower > upper)
    if (length(bad)) {
        stop('Calibration bounds are temporally infeasible at node(s): ',
             format_limited_values(get_node_name_by_num(phy, bad)),
             '. Bounds are preserved; revise the input constraints',
             ' or the requested minimum branch length.')
    }
    list(lower=lower, upper=upper, min_edge=min_edge, preorder=preorder)
}

dated_tree_validation = function(tree, reference=NULL, calibration=NULL, tolerance=1e-7) {
    invalid = function(reason) list(valid=FALSE, reason=reason)
    if (!inherits(tree, 'phylo') || is.null(tree$edge) || ncol(tree$edge) != 2L) {
        return(invalid('No valid phylo edge table was returned.'))
    }
    tips = length(tree$tip.label)
    total = tips + tree$Nnode
    if (tips < 2L || anyNA(tree$tip.label) || anyDuplicated(tree$tip.label) ||
        any(!nzchar(tree$tip.label)) || anyNA(tree$edge) ||
        any(tree$edge < 1L | tree$edge > total) || any(tree$edge != as.integer(tree$edge)) ||
        nrow(tree$edge) != total - 1L || anyDuplicated(tree$edge[, 2]) ||
        any(tree$edge[, 1] <= tips) || length(get_root_num(tree)) != 1L) {
        return(invalid('Tree node numbers, tip labels, or rooted structure are invalid.'))
    }
    # Detect disconnected cycles before entering ape's C-level traversals.
    children = split(tree$edge[, 2], factor(tree$edge[, 1], levels=seq_len(total)))
    queue = integer(total)
    queue[[1]] = get_root_num(tree)
    seen = logical(total)
    end = 1L
    cursor = 1L
    while (cursor <= end) {
        node = queue[[cursor]]
        if (seen[[node]]) return(invalid('Tree contains a cycle.'))
        seen[[node]] = TRUE
        child = children[[node]]
        if (length(child)) {
            if (end + length(child) > total) return(invalid('Tree contains a cycle.'))
            queue[seq.int(end + 1L, end + length(child))] = child
            end = end + length(child)
        }
        cursor = cursor + 1L
    }
    if (!all(seen)) return(invalid('Tree contains disconnected nodes.'))
    if (length(tree$edge.length) != nrow(tree$edge) ||
        any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) {
        return(invalid('Edge lengths must be finite and positive.'))
    }
    depths = ape::node.depth.edgelength(tree)
    height = max(depths[seq_len(tips)])
    eps = tolerance * max(1, height)
    if (diff(range(depths[seq_len(tips)])) > eps) return(invalid('Dated tree is not ultrametric.'))
    ages = height - depths
    mapped_nodes = seq_len(total)
    if (!is.null(reference)) {
        if (!setequal(reference$tip.label, tree$tip.label) || reference$Nnode != tree$Nnode) {
            return(invalid('Output tip set or topology differs from the input gene tree.'))
        }
        internal = match_internal_nodes_by_clade(tree, reference)
        if (anyNA(internal)) return(invalid('Output topology differs from the input gene tree.'))
        mapped_nodes = c(match(reference$tip.label, tree$tip.label), tips + internal)
    }
    if (!is.null(calibration) && nrow(calibration)) {
        observed = ages[mapped_nodes[calibration$node]]
        hard = if ('soft.bounds' %in% names(calibration)) !calibration$soft.bounds %in% TRUE else rep(TRUE, nrow(calibration))
        bad = hard & (is.na(observed) | observed < calibration$age.min - eps | observed > calibration$age.max + eps)
        if (any(bad)) {
            return(invalid(paste('Estimated ages violate original calibration bounds at node(s):',
                                paste(calibration$node[bad], collapse=', '))))
        }
    }
    list(valid=TRUE, reason='', ages=ages, mapped_nodes=mapped_nodes)
}

validate_dated_tree = function(tree, reference=NULL, calibration=NULL, tolerance=1e-7) {
    result = dated_tree_validation(tree, reference, calibration, tolerance)
    if (!result$valid) stop('Invalid dated tree: ', result$reason)
    invisible(result)
}

transfer_dated_ages = function(reference, dated, tolerance=1e-7) {
    validation = validate_dated_tree(dated, reference, tolerance=tolerance)
    ages = validation$ages[validation$mapped_nodes]
    out = reference
    out$edge.length = ages[out$edge[, 1]] - ages[out$edge[, 2]]
    out
}

transfer_species_ages = function(gn_tree, sp_tree, species_parser) {
    species = species_parser_get_gene_species(species_parser, gn_tree$tip.label, sp_tree$tip.label, strict=TRUE)
    if (anyDuplicated(species)) stop('allS requires exactly one matching tip per species; annotate duplications.')
    species_tips = resolve_species_tree_tips(species_parser, species, sp_tree$tip.label)
    if (anyNA(species_tips)) stop('Species tree tip mapping failed for allS.')
    pruned = if (length(species_tips) < length(sp_tree$tip.label)) {
        ape::keep.tip(sp_tree, species_tips)
    } else sp_tree
    pruned$tip.label = gn_tree$tip.label[match(pruned$tip.label, species_tips)]
    # Keep the gene tree's node numbers as well as its topology and labels.
    transfer_dated_ages(gn_tree, pruned)
}

annotate_calibration_fit = function(calibration, tree, backend, selected) {
    depths = ape::node.depth.edgelength(tree)
    ages = max(depths) - depths
    eps = max(1, max(ages)) * if (backend == 'mcmctree') 1e-5 else 1e-7
    calibration$estimated_age = ages[calibration$node]
    calibration$within_original_bounds = calibration$estimated_age >= calibration$age.min - eps &
        calibration$estimated_age <= calibration$age.max + eps
    calibration$constraint_status = ifelse(calibration$node %in% selected, 'retained',
        ifelse(grepl('^S|R', calibration$event), 'dropped', 'not-selected'))
    calibration$bound_policy = if (backend == 'mcmctree') 'PAML-soft-prior' else 'hard'
    calibration$soft.bounds = backend == 'mcmctree'
    calibration
}
