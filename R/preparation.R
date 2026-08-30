read_radte_species_context = function(args, species_parser) {
    sp_file = args$species_tree
    species_node_bounds_tsv = args$species_node_bounds_tsv
    cat('\nStart: species tree processing', '\n')
    sp_tree = ape::read.tree(sp_file)
    sp_tree = validate_species_tree_labels(sp_tree)
    sp_tip_species = validate_species_tip_parser_labels(sp_tree, species_parser)
    if (contains_polytomy(sp_tree)) {
        stop('Input species tree contains polytomy. A completely bifurcated tree is expected as input.')
    }
    validate_tree_edge_lengths(sp_tree, 'species tree')
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

    list(tree=sp_tree, nodes=sp_node_table, has_intervals=has_species_node_interval_input)
}

read_radte_gene_context = function(args, mode, sp_tree, sp_node_table, species_parser) {
    generax_file = args$generax_nhx
    gn_file = args$gene_tree
    parsable_file = args$notung_parsable
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
        validate_tree_edge_lengths(gn_tree, 'gene tree')
        gn_tree = pad_branch_length(gn_tree, pad_size=args[['pad_short_edge']])
        #gn_tree = adjust_branch_length_order(gn_tree, min_bl=args[['pad_short_edge']])
        cat('Minimum branch length in gene tree:', min(gn_tree[['edge.length']]), '\n')
        cols = c('event', 'gn_node', 'gn_node_num', 'lower_sp_node', 'upper_sp_node')
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
        gn_node_table[,'gn_node_num'] = get_node_num_by_name(gn_tree, gn_node_table[['gn_node']])
        gn_node_table = data.frame(gn_node_table[,cols], stringsAsFactors=FALSE)
        validate_gn_species_nodes(gn_node_table, sp_tree)
    } else if (mode=='notung') {
        cat('Reading NOTUNG tree.\n')
        gn_tree = ape::read.tree(gn_file)
        if (contains_polytomy(gn_tree)) {
            stop('Input tree contains polytomy. A completely bifurcated tree is expected as input.')
        }
        gn_tree = validate_gene_tree_labels(gn_tree)
        validate_gene_tip_labels(gn_tree[['tip.label']], species_parser=species_parser, species_tree_labels=sp_tree[['tip.label']])
        validate_gene_species_coverage(gn_tree, sp_tree, species_parser=species_parser, species_tree_labels=sp_tree[['tip.label']])
        validate_tree_edge_lengths(gn_tree, 'gene tree')
        gn_tree = pad_branch_length(gn_tree, pad_size=args[['pad_short_edge']])

        gn_node_table = read_notung_parsable(file=parsable_file, mode='D')
        gn_node_table = merge(gn_node_table, data.frame(lower_age=NA, upper_age=NA, spp=NA), all=TRUE)
        check_gn_node_name_uniqueness(gn_node_table, gn_tree)
        if (nrow(gn_node_table) > 0) {
            gn_node_table$gn_node_num = get_node_num_by_name(
                gn_tree,
                gn_node_table[['gn_node']]
            )
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

    gn_node_table = resolve_gene_node_bounds(gn_node_table, sp_tree, sp_node_table, args$max_age)
    gn_node_table = annotate_species_constraint_groups(gn_node_table, gn_tree)
    validate_gn_node_table(gn_node_table)
    validate_duplication_nodes_internal(gn_node_table, gn_tree)
    root_num = get_root_num(gn_tree)
    root_row = which(gn_node_table$gn_node_num == root_num)
    if (length(root_row) != 1L) stop('Gene node table root event mapping failed: expected exactly one row for root node.')
    gn_node_table$event[root_row] = ensure_root_event_tag(gn_node_table$event[root_row])
    cat('End: gene tree processing\n\n')
    list(tree=gn_tree, nodes=gn_node_table)
}

resolve_gene_node_bounds = function(nodes, sp_tree, species_nodes, max_age) {
    validate_gn_species_nodes(nodes, sp_tree)
    lower = match(nodes$lower_sp_node, species_nodes$node)
    upper = match(nodes$upper_sp_node, species_nodes$node)
    above_root = grepl('^D', nodes$event) & is.na(upper) & lower == get_root_num(sp_tree)
    above_root[is.na(above_root)] = FALSE
    unresolved = is.na(lower) | (is.na(upper) & !above_root)
    if (any(unresolved)) stop('Gene node table contains unresolved age bound(s) for node(s): ',
                              paste(nodes$gn_node[unresolved], collapse=', '),
                              '. A missing upper species node is allowed only for a duplication above the species-tree root.')
    nodes$lower_age = species_nodes$age_min[lower]
    nodes$upper_age = species_nodes$age_max[upper]
    nodes$upper_age[above_root] = max_age
    if (any(nodes$lower_age > nodes$upper_age)) stop('Gene node table contains reversed age bounds.')
    nodes
}

prepare_radte_calibrations = function(gene_context) {
    gn_tree = gene_context$tree
    gn_node_table = gene_context$nodes
    root_num = get_root_num(gn_tree)
    gn_node_table_calibration = gn_node_table[gn_node_table$gn_node_num > ape::Ntip(gn_tree), , drop=FALSE]
    # Calibration table
    calibration_table = data.frame(
        node=as.integer(gn_node_table_calibration$gn_node_num),
        age.min=as.numeric(gn_node_table_calibration$lower_age),
        age.max=as.numeric(gn_node_table_calibration$upper_age),
        soft.bounds=FALSE,
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

    calibration_feasible_ranges(gn_tree, calibration_tables[['RS']])

    list(all=calibration_table, tables=calibration_tables)
}
