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

is_nonempty_scalar_string = function(value) {
    if (is.null(value) || (length(value) != 1) || is.na(value)) {
        return(FALSE)
    }
    return(nchar(trimws(as.character(value))) > 0)
}

get_species_parser_names = function() {
    return(c('legacy', 'taxonomic', 'regex', 'map'))
}

detect_species_map_has_header = function(file) {
    lines = readLines(file, warn=FALSE)
    lines = lines[nchar(trimws(lines)) > 0]
    if (length(lines) == 0) {
        stop('Species map TSV is empty.')
    }
    first_fields = strsplit(lines[[1]], '\t', fixed=TRUE)[[1]]
    first_fields_norm = gsub('[^a-z0-9]+', '', tolower(first_fields))
    label_candidates = c('label', 'tip', 'tiplabel', 'genetip', 'genelabel', 'input', 'query', 'name', 'alias')
    species_candidates = c('species', 'specieslabel', 'canonicalspecies', 'mappedspecies', 'canonical')
    return(
        any(first_fields_norm %in% label_candidates) &&
        any(first_fields_norm %in% species_candidates)
    )
}

find_species_map_column = function(colnames_norm, candidates) {
    matched = which(colnames_norm %in% candidates)
    if (length(matched) == 0) {
        return(NA_integer_)
    }
    return(matched[[1]])
}

read_species_map_tsv = function(file) {
    if (!is_nonempty_scalar_string(file)) {
        stop('--species-map-tsv should be a non-empty path.')
    }
    file = as.character(file)
    if (!file.exists(file)) {
        stop('Species map TSV was not found: ', file)
    }

    has_header = detect_species_map_has_header(file)
    map_df = utils::read.delim(
        file=file,
        sep='\t',
        header=has_header,
        stringsAsFactors=FALSE,
        colClasses='character',
        na.strings='',
        check.names=FALSE,
        quote='',
        comment.char=''
    )
    if (ncol(map_df) < 2) {
        stop('Species map TSV should contain at least two tab-delimited columns.')
    }

    if (!has_header) {
        colnames(map_df)[1:2] = c('label', 'species')
        if (ncol(map_df) >= 3) {
            colnames(map_df)[[3]] = 'taxonomy_query'
        }
    } else {
        colnames_norm = gsub('[^a-z0-9]+', '', tolower(colnames(map_df)))
        label_idx = find_species_map_column(
            colnames_norm,
            c('label', 'tip', 'tiplabel', 'genetip', 'genelabel', 'input', 'query', 'name', 'alias')
        )
        species_idx = find_species_map_column(
            colnames_norm,
            c('species', 'specieslabel', 'canonicalspecies', 'mappedspecies', 'canonical')
        )
        taxonomy_idx = find_species_map_column(
            colnames_norm,
            c('taxonomyquery', 'scientificname', 'scientific', 'taxonomy', 'taxquery')
        )
        if (is.na(label_idx) || is.na(species_idx)) {
            stop(
                'Species map TSV should contain label/species columns. ',
                'Accepted names include: label, gene_tip, species, species_label.'
            )
        }
        selected = data.frame(
            label=map_df[[label_idx]],
            species=map_df[[species_idx]],
            stringsAsFactors=FALSE
        )
        if (!is.na(taxonomy_idx)) {
            selected[['taxonomy_query']] = map_df[[taxonomy_idx]]
        }
        map_df = selected
    }

    map_df[['label']] = trimws(as.character(map_df[['label']]))
    map_df[['species']] = trimws(as.character(map_df[['species']]))
    if ('taxonomy_query' %in% colnames(map_df)) {
        map_df[['taxonomy_query']] = trimws(as.character(map_df[['taxonomy_query']]))
    }
    is_invalid = (
        is.na(map_df[['label']]) | (map_df[['label']] == '') |
        is.na(map_df[['species']]) | (map_df[['species']] == '')
    )
    if (any(is_invalid)) {
        stop('Species map TSV contains empty label/species value(s).')
    }
    if (any(duplicated(map_df[['label']]))) {
        duplicated_labels = unique(map_df[['label']][duplicated(map_df[['label']])])
        stop('Species map TSV contains duplicated label(s): ', paste(duplicated_labels, collapse=', '))
    }
    return(map_df)
}

read_species_node_bounds_tsv = function(file) {
    if (!is_nonempty_scalar_string(file)) {
        stop('--species_node_bounds_tsv should be a non-empty path.')
    }
    file = as.character(file)
    if (!file.exists(file)) {
        stop('Species node bounds TSV was not found: ', file)
    }

    bounds_df = utils::read.delim(
        file=file,
        sep='\t',
        header=TRUE,
        stringsAsFactors=FALSE,
        colClasses='character',
        na.strings='',
        check.names=FALSE,
        quote='',
        comment.char=''
    )
    if (ncol(bounds_df) < 2) {
        stop('Species node bounds TSV should contain at least two tab-delimited columns with a header.')
    }

    colnames_norm = gsub('[^a-z0-9]+', '', tolower(colnames(bounds_df)))
    node_idx = find_species_map_column(
        colnames_norm,
        c('node', 'nodelabel', 'speciesnode', 'speciesnodelabel', 'spnode')
    )
    age_idx = find_species_map_column(
        colnames_norm,
        c('age', 'exactage', 'pointage')
    )
    min_idx = find_species_map_column(
        colnames_norm,
        c('agemin', 'minage', 'lowerage', 'lowerbound', 'minimumage')
    )
    max_idx = find_species_map_column(
        colnames_norm,
        c('agemax', 'maxage', 'upperage', 'upperbound', 'maximumage')
    )
    if (is.na(node_idx)) {
        stop(
            'Species node bounds TSV should contain a node column. ',
            'Accepted names include: node, node_label, species_node.'
        )
    }
    has_exact_age = !is.na(age_idx)
    has_interval_age = (!is.na(min_idx) || !is.na(max_idx))
    if (has_exact_age && has_interval_age) {
        stop('Species node bounds TSV should use either age or age_min/age_max columns, not both.')
    }
    if ((!has_exact_age) && (is.na(min_idx) || is.na(max_idx))) {
        stop(
            'Species node bounds TSV should contain either an age column or both age_min and age_max columns.'
        )
    }

    selected = data.frame(
        node=trimws(as.character(bounds_df[[node_idx]])),
        stringsAsFactors=FALSE
    )
    if (has_exact_age) {
        exact_age = suppressWarnings(as.numeric(bounds_df[[age_idx]]))
        selected[['age_min']] = exact_age
        selected[['age_max']] = exact_age
    } else {
        selected[['age_min']] = suppressWarnings(as.numeric(bounds_df[[min_idx]]))
        selected[['age_max']] = suppressWarnings(as.numeric(bounds_df[[max_idx]]))
    }

    has_invalid_node = is.na(selected[['node']]) | (selected[['node']] == '')
    has_invalid_age = (
        is.na(selected[['age_min']]) | !is.finite(selected[['age_min']]) |
        is.na(selected[['age_max']]) | !is.finite(selected[['age_max']])
    )
    if (any(has_invalid_node)) {
        stop('Species node bounds TSV contains empty node value(s).')
    }
    if (any(has_invalid_age)) {
        stop('Species node bounds TSV contains missing or non-finite age bound(s).')
    }
    if (any(selected[['age_min']] < 0) || any(selected[['age_max']] < 0)) {
        stop('Species node bounds TSV contains negative age bound(s).')
    }
    if (any(selected[['age_max']] < selected[['age_min']])) {
        stop('Species node bounds TSV contains age_max younger than age_min.')
    }
    if (any(duplicated(selected[['node']]))) {
        duplicated_nodes = unique(selected[['node']][duplicated(selected[['node']])])
        stop('Species node bounds TSV contains duplicated node(s): ', paste(duplicated_nodes, collapse=', '))
    }
    return(selected)
}

merge_species_node_bounds = function(sp_node_table, bounds_df) {
    if (nrow(bounds_df) == 0) {
        return(sp_node_table)
    }
    missing_nodes = setdiff(bounds_df[['node']], sp_node_table[['node']])
    if (length(missing_nodes) > 0) {
        stop(
            'Species node bounds TSV contains node(s) not found in the species tree: ',
            paste(missing_nodes, collapse=', ')
        )
    }
    row_idx = match(bounds_df[['node']], sp_node_table[['node']])
    sp_node_table[row_idx, 'age_min'] = bounds_df[['age_min']]
    sp_node_table[row_idx, 'age_max'] = bounds_df[['age_max']]

    point_age = suppressWarnings(as.numeric(sp_node_table[['age']]))
    is_outside = (point_age < sp_node_table[['age_min']]) | (point_age > sp_node_table[['age_max']])
    is_outside[is.na(is_outside)] = FALSE
    if (any(is_outside)) {
        invalid_nodes = sp_node_table[is_outside, 'node']
        stop(
            'Species node bounds TSV is inconsistent with the species-tree branch-length ages for node(s): ',
            paste(invalid_nodes, collapse=', '),
            '. The point age from the tree must fall within [age_min, age_max].'
        )
    }
    return(sp_node_table)
}

build_species_parser = function(parser_name='legacy', species_regex=NULL, species_map_tsv=NULL, sep='_') {
    if (!is_nonempty_scalar_string(parser_name)) {
        stop('--species-parser should be a non-empty string.')
    }
    parser_name = tolower(trimws(as.character(parser_name)))
    supported = get_species_parser_names()
    if (!parser_name %in% supported) {
        stop('--species-parser should be one of: ', paste(supported, collapse=', '))
    }

    has_species_regex = is_nonempty_scalar_string(species_regex)
    has_species_map = is_nonempty_scalar_string(species_map_tsv)
    if ((parser_name == 'regex') && (!has_species_regex)) {
        stop('--species-regex is required when --species-parser=regex.')
    }
    if ((parser_name != 'regex') && has_species_regex) {
        stop('--species-regex can only be used with --species-parser=regex.')
    }
    if ((parser_name == 'map') && (!has_species_map)) {
        stop('--species-map-tsv is required when --species-parser=map.')
    }
    if ((parser_name != 'map') && has_species_map) {
        stop('--species-map-tsv can only be used with --species-parser=map.')
    }

    species_map = NULL
    if (parser_name == 'map') {
        species_map = read_species_map_tsv(species_map_tsv)
    }
    if (parser_name != 'regex') {
        species_regex = NULL
    }
    return(
        list(
            name=parser_name,
            sep=sep,
            species_regex=species_regex,
            species_map=species_map
        )
    )
}

get_invalid_tip_label_message = function(species_parser, invalid_labels) {
    max_show = min(5, length(invalid_labels))
    shown_labels = paste(invalid_labels[seq_len(max_show)], collapse=', ')
    extra_suffix = ''
    if (length(invalid_labels) > max_show) {
        extra_suffix = paste0(' ... (', length(invalid_labels) - max_show, ' more)')
    }
    parser_name = species_parser[['name']]
    if (parser_name == 'legacy') {
        return(
            paste0(
                'Input gene tree tip label(s) must follow GENUS_SPECIES_GENEID format. Invalid label(s): ',
                shown_labels,
                extra_suffix
            )
        )
    }
    if (parser_name == 'regex') {
        return(
            paste0(
                'Input gene tree tip label(s) are invalid for --species-parser=regex. ',
                'They must match --species-regex and produce a non-empty species label. Invalid label(s): ',
                shown_labels,
                extra_suffix
            )
        )
    }
    if (parser_name == 'map') {
        return(
            paste0(
                'Input gene tree tip label(s) are invalid for --species-parser=map. ',
                'They must be present in --species-map-tsv. Invalid label(s): ',
                shown_labels,
                extra_suffix
            )
        )
    }
    return(
        paste0(
            'Input gene tree tip label(s) are invalid for --species-parser=',
            parser_name,
            '. Invalid label(s): ',
            shown_labels,
            extra_suffix
        )
    )
}

infer_taxonomic_species_token_count = function(tokens) {
    if (length(tokens) < 2) {
        return(0)
    }
    token2 = tolower(tokens[[2]])
    token3 = ''
    if (length(tokens) >= 3) {
        token3 = tolower(tokens[[3]])
    }
    genus_only_placeholders = c('sp', 'spp')
    single_token_qualifiers = c('cf', 'aff', 'nr', 'x')
    paired_token_qualifiers = c(
        'subsp', 'ssp', 'subspecies',
        'var', 'variety',
        'forma', 'form', 'f',
        'strain', 'substrain',
        'serovar', 'serotype', 'serogroup',
        'pathovar', 'pv',
        'biovar', 'biotype', 'chemovar', 'morphovar',
        'cultivar', 'cv',
        'isolate',
        'group', 'subgroup', 'complex', 'clade', 'lineage',
        'section', 'series', 'ecotype', 'breed'
    )
    if ((length(tokens) >= 5) && (token3 %in% paired_token_qualifiers)) {
        return(4)
    }
    if ((length(tokens) >= 4) && ((token2 %in% c(genus_only_placeholders, single_token_qualifiers)) || (token3 %in% c(single_token_qualifiers, genus_only_placeholders)))) {
        return(3)
    }
    return(2)
}

match_species_prefix_from_tree = function(tip_label, species_tree_labels) {
    if (length(species_tree_labels) == 0) {
        return(character(0))
    }
    matches = species_tree_labels[
        (tip_label == species_tree_labels) |
        startsWith(tip_label, paste0(species_tree_labels, '_'))
    ]
    if (length(matches) == 0) {
        return(character(0))
    }
    match_lengths = nchar(matches)
    longest_matches = unique(matches[match_lengths == max(match_lengths)])
    return(longest_matches)
}

species_parser_get_gene_species = function(species_parser, tip_labels, species_tree_labels=NULL, strict=TRUE) {
    if (is.null(species_parser)) {
        species_parser = build_species_parser('legacy')
    }
    tip_labels = as.character(tip_labels)
    if (length(tip_labels) == 0) {
        return(character(0))
    }

    parser_name = species_parser[['name']]
    if (parser_name == 'legacy') {
        split_labels = strsplit(tip_labels, species_parser[['sep']], fixed=TRUE)
        species_names = vapply(split_labels, function(items) {
            if (length(items) < 3) {
                return(NA_character_)
            }
            genus = items[[1]]
            species = items[[2]]
            gene_id = paste(items[3:length(items)], collapse=species_parser[['sep']])
            if ((nchar(genus) == 0) || (nchar(species) == 0) || (nchar(gene_id) == 0)) {
                return(NA_character_)
            }
            return(paste0(genus, species_parser[['sep']], species))
        }, character(1))
    } else if (parser_name == 'taxonomic') {
        split_labels = strsplit(tip_labels, '_', fixed=TRUE)
        species_names = vapply(seq_along(split_labels), function(i) {
            items = split_labels[[i]]
            if (length(items) < 3) {
                return(NA_character_)
            }
            if (any(items == '')) {
                return(NA_character_)
            }
            if (length(species_tree_labels) > 0) {
                matched_species = match_species_prefix_from_tree(tip_labels[[i]], species_tree_labels)
                matched_species = matched_species[nchar(matched_species) < nchar(tip_labels[[i]])]
                if (length(matched_species) == 1) {
                    return(matched_species[[1]])
                }
                if (length(matched_species) > 1) {
                    return(NA_character_)
                }
            }
            species_token_count = infer_taxonomic_species_token_count(items)
            if ((species_token_count < 2) || (length(items) <= species_token_count)) {
                return(NA_character_)
            }
            gene_id = paste(items[(species_token_count + 1):length(items)], collapse='_')
            if (nchar(gene_id) == 0) {
                return(NA_character_)
            }
            return(paste(items[seq_len(species_token_count)], collapse='_'))
        }, character(1))
    } else if (parser_name == 'regex') {
        regex_match = regexec(species_parser[['species_regex']], tip_labels, perl=TRUE)
        regex_items = regmatches(tip_labels, regex_match)
        species_names = vapply(regex_items, function(items) {
            if (length(items) == 0) {
                return(NA_character_)
            }
            species_label = items[[1]]
            if (length(items) >= 2) {
                species_label = items[[2]]
            }
            if (is.na(species_label) || (nchar(species_label) == 0)) {
                return(NA_character_)
            }
            return(species_label)
        }, character(1))
    } else if (parser_name == 'map') {
        species_map = species_parser[['species_map']]
        map_idx = match(tip_labels, species_map[['label']])
        species_names = species_map[['species']][map_idx]
        species_names[is.na(map_idx)] = NA_character_
    } else {
        stop('Unsupported species parser: ', parser_name)
    }

    if (strict && any(is.na(species_names))) {
        invalid_labels = tip_labels[is.na(species_names)]
        stop(get_invalid_tip_label_message(species_parser, invalid_labels))
    }
    return(species_names)
}

species_parser_get_species_tip_labels = function(species_parser, tip_labels, strict=TRUE) {
    if (is.null(species_parser)) {
        species_parser = build_species_parser('legacy')
    }
    tip_labels = as.character(tip_labels)
    if (length(tip_labels) == 0) {
        return(character(0))
    }

    parser_name = species_parser[['name']]
    if (parser_name == 'map') {
        species_map = species_parser[['species_map']]
        map_idx = match(tip_labels, species_map[['label']])
        species_names = tip_labels
        is_mapped = !is.na(map_idx)
        species_names[is_mapped] = species_map[['species']][map_idx[is_mapped]]
        known_species = unique(species_map[['species']])
        is_valid = is_mapped | (tip_labels %in% known_species)
        if (strict && any(!is_valid)) {
            invalid_labels = tip_labels[!is_valid]
            stop(
                'Input species tree tip label(s) are invalid for --species-parser=map. ',
                'They must be present in --species-map-tsv or match mapped species labels. Invalid label(s): ',
                format_limited_values(invalid_labels, max_items=5)
            )
        }
        species_names[!is_valid] = NA_character_
        return(species_names)
    }
    return(tip_labels)
}

validate_species_tip_parser_labels = function(sp_tree, species_parser) {
    parsed_tip_labels = species_parser_get_species_tip_labels(species_parser, sp_tree[['tip.label']], strict=TRUE)
    if (any(duplicated(parsed_tip_labels))) {
        duplicated_labels = unique(parsed_tip_labels[duplicated(parsed_tip_labels)])
        stop(
            'Input species tree contains duplicated species label(s) after applying the species parser: ',
            paste(duplicated_labels, collapse=', ')
        )
    }
    return(parsed_tip_labels)
}

species_parser_taxonomy_query = function(species_parser, species_labels) {
    if (is.null(species_parser)) {
        species_parser = build_species_parser('legacy')
    }
    species_labels = as.character(species_labels)
    if (length(species_labels) == 0) {
        return(character(0))
    }

    taxonomy_names = gsub('_', ' ', species_labels, fixed=TRUE)
    if (species_parser[['name']] == 'map') {
        species_map = species_parser[['species_map']]
        if ('taxonomy_query' %in% colnames(species_map)) {
            lookup_idx = match(species_labels, species_map[['species']])
            has_lookup = !is.na(lookup_idx)
            mapped_taxonomy = species_map[['taxonomy_query']][lookup_idx[has_lookup]]
            valid_taxonomy = !is.na(mapped_taxonomy) & (mapped_taxonomy != '')
            if (any(valid_taxonomy)) {
                matched_positions = which(has_lookup)
                taxonomy_names[matched_positions[valid_taxonomy]] = mapped_taxonomy[valid_taxonomy]
            }
        }
    }
    return(taxonomy_names)
}

resolve_species_tree_tips = function(species_parser, species_labels, species_tree_labels) {
    species_tree_labels = as.character(species_tree_labels)
    parsed_species_tree = species_parser_get_species_tip_labels(species_parser, species_tree_labels, strict=TRUE)
    if (any(duplicated(parsed_species_tree))) {
        duplicated_labels = unique(parsed_species_tree[duplicated(parsed_species_tree)])
        stop(
            'Species tree tip mapping is ambiguous after applying the species parser: ',
            paste(duplicated_labels, collapse=', ')
        )
    }
    resolved = vapply(species_labels, function(species_label) {
        matched_tip = species_tree_labels[parsed_species_tree == species_label]
        if (length(matched_tip) != 1) {
            return(NA_character_)
        }
        return(matched_tip[[1]])
    }, character(1))
    return(resolved)
}

# copied from rkftools https://github.com/kfuku52/rkftools
get_species_names = function(phy, sep='_', species_parser=NULL, species_tree_labels=NULL) {
    if (is.null(species_parser)) {
        species_parser = build_species_parser('legacy', sep=sep)
    }
    return(
        species_parser_get_gene_species(
            species_parser=species_parser,
            tip_labels=phy[['tip.label']],
            species_tree_labels=species_tree_labels,
            strict=TRUE
        )
    )
}

# copied from rkftools https://github.com/kfuku52/rkftools
leaf2species = function(leaf_names, species_parser=NULL, species_tree_labels=NULL) {
    if (is.null(species_parser)) {
        species_parser = build_species_parser('legacy')
    }
    if (length(leaf_names) == 0) {
        return(character(0))
    }
    species_names = species_parser_get_gene_species(
        species_parser=species_parser,
        tip_labels=leaf_names,
        species_tree_labels=species_tree_labels,
        strict=FALSE
    )
    invalid_labels = leaf_names[is.na(species_names)]
    if (length(invalid_labels) > 0) {
        for (invalid_label in invalid_labels) {
            warning('leaf name could not be interpreted by species parser: ', invalid_label, '\n')
        }
    }
    species_names = species_names[!is.na(species_names)]
    return(species_parser_taxonomy_query(species_parser, species_names))
}

unquote_newick_labels = function(labels) {
    if (is.null(labels)) return(NULL)
    quoted = !is.na(labels) & startsWith(labels, "'") & endsWith(labels, "'") & nchar(labels) >= 2L
    labels[quoted] = gsub("''", "'", substring(labels[quoted], 2L, nchar(labels[quoted]) - 1L), fixed=TRUE)
    labels
}

validate_species_tree_labels = function(sp_tree) {
    sp_tree$tip.label = unquote_newick_labels(sp_tree$tip.label)
    sp_tree$node.label = unquote_newick_labels(sp_tree$node.label)
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

validate_gene_tip_labels = function(tip_labels, species_parser=NULL, species_tree_labels=NULL) {
    if (is.null(species_parser)) {
        species_parser = build_species_parser('legacy')
    }
    invisible(
        species_parser_get_gene_species(
            species_parser=species_parser,
            tip_labels=tip_labels,
            species_tree_labels=species_tree_labels,
            strict=TRUE
        )
    )
}

validate_gene_species_coverage = function(gn_tree, sp_tree, species_parser=NULL, species_tree_labels=NULL) {
    if (is.null(species_parser)) {
        species_parser = build_species_parser('legacy')
    }
    if (is.null(species_tree_labels)) {
        species_tree_labels = sp_tree[['tip.label']]
    }
    gn_species = unique(
        get_species_names(
            phy=gn_tree,
            species_parser=species_parser,
            species_tree_labels=species_tree_labels
        )
    )
    sp_species = unique(
        species_parser_get_species_tip_labels(
            species_parser=species_parser,
            tip_labels=species_tree_labels,
            strict=TRUE
        )
    )
    missing_species = setdiff(gn_species, sp_species)
    if (length(missing_species) > 0) {
        stop(
            'Species in the gene tree were not found in the species tree: ',
            paste(missing_species, collapse=', ')
        )
    }
}

validate_gene_tree_labels = function(gn_tree) {
    gn_tree$tip.label = unquote_newick_labels(gn_tree$tip.label)
    gn_tree$node.label = unquote_newick_labels(gn_tree$node.label)
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

annotate_species_constraint_groups = function(gn_node_table, gn_tree) {
    if (nrow(gn_node_table) == 0) {
        gn_node_table[['constraint_sp_node']] = character(0)
        gn_node_table[['shared_speciation_group']] = character(0)
        return(gn_node_table)
    }

    gn_node_table[['constraint_sp_node']] = NA_character_
    is_exact_speciation = (
        grepl('^S', gn_node_table[['event']]) &
        !is.na(gn_node_table[['lower_sp_node']]) &
        !is.na(gn_node_table[['upper_sp_node']]) &
        (gn_node_table[['lower_sp_node']] == gn_node_table[['upper_sp_node']])
    )
    is_exact_speciation[is.na(is_exact_speciation)] = FALSE
    gn_node_table[is_exact_speciation, 'constraint_sp_node'] = gn_node_table[is_exact_speciation, 'lower_sp_node']

    gn_node_table[['shared_speciation_group']] = NA_character_
    is_internal_exact_speciation = is_exact_speciation & (gn_node_table[['gn_node_num']] > ape::Ntip(gn_tree))
    constraint_nodes = gn_node_table[is_internal_exact_speciation, 'constraint_sp_node']
    if (length(constraint_nodes) > 0) {
        group_sizes = table(constraint_nodes)
        shared_nodes = names(group_sizes[group_sizes >= 2])
        is_shared = is_internal_exact_speciation & (gn_node_table[['constraint_sp_node']] %in% shared_nodes)
        is_shared[is.na(is_shared)] = FALSE
        gn_node_table[is_shared, 'shared_speciation_group'] = gn_node_table[is_shared, 'constraint_sp_node']
    }
    return(gn_node_table)
}

calibration_table_has_interval_bounds = function(calibration_table) {
    if (nrow(calibration_table) == 0) {
        return(FALSE)
    }
    age_min = suppressWarnings(as.numeric(calibration_table[['age.min']]))
    age_max = suppressWarnings(as.numeric(calibration_table[['age.max']]))
    is_interval = (!is.na(age_min)) & (!is.na(age_max)) & (abs(age_max - age_min) > 1e-12)
    return(any(is_interval))
}

gn_node_table_has_shared_interval_speciation = function(gn_node_table) {
    required_cols = c('shared_speciation_group', 'lower_age', 'upper_age')
    if ((nrow(gn_node_table) == 0) || any(!required_cols %in% colnames(gn_node_table))) {
        return(FALSE)
    }
    is_shared = !is.na(gn_node_table[['shared_speciation_group']]) & (gn_node_table[['shared_speciation_group']] != '')
    if (!any(is_shared)) {
        return(FALSE)
    }
    age_min = suppressWarnings(as.numeric(gn_node_table[['lower_age']]))
    age_max = suppressWarnings(as.numeric(gn_node_table[['upper_age']]))
    is_interval = (!is.na(age_min)) & (!is.na(age_max)) & (abs(age_max - age_min) > 1e-12)
    return(any(is_shared & is_interval))
}
