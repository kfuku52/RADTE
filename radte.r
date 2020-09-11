#!/usr/bin/env Rscript
library(ape)
library(rkftools)

run_mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')
cat('RADTE run_mode:', run_mode, '\n')
if (run_mode=='batch') {
    cat('arguments:\n')
    args = commandArgs(trailingOnly=TRUE)
    args = rkftools::get_parsed_args(args, print=TRUE)
} else if (run_mode=='debug') {
    test_type = 'generax'
    args = list()
    args[['max_age']] = 1000
    args[['chronos_lambda']] = 1
    args[['chronos_model']] = 'discrete'
    args[['pad_short_edge']] = 0.001
    if (test_type=='notung') {
        args[['species_tree']] = "/Users/kef74yk/Dropbox_p/repos/RADTE/data/phyloxml1/input/dated_species_tree.generax.nwk"
        args[['notung_parsable']] = ""
    }
    if (test_type=='generax') {
        dir_work = '/Users/kef74yk/Downloads/'
        setwd(dir_work)
        args[['species_tree']] = '/Users/kef74yk/Dropbox (Personal)/repos/gfe_pipeline/gfe_data/species_tree/dated_species_tree.pruned.nwk'
        args[['generax_nhx']] = '/Users/kef74yk/Downloads/orthogroup/generax.tree/OG0000002.generax.nhx'
        #args[['species_tree']] = paste0(dir_work, "species_tree.nwk")
        #args[['generax_nhx']] = paste0(dir_work, "gene_tree.nhx")
    }
}

if ('generax_nhx' %in% names(args)) {
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
sp_file = args[['species_tree']]
max_age = as.numeric(args[['max_age']])
chronos_lambda = as.numeric(args[['chronos_lambda']])
chronos_model = args[['chronos_model']]

cat('\nStart: species tree processing', '\n')
tree_text0 = scan(sp_file, what=character(), sep="\n", blank.lines.skip=FALSE)
tree_text1 = gsub("'([0-9]+)'", "PLACEHOLDER\\1", tree_text0)
sp_tree = read.tree(text=tree_text1)
if (all(is.na(sp_tree$node.label))) {
    sp_tree
} else {
    sp_tree[['node.label']] = sub('PLACEHOLDER', '', sp_tree[['node.label']])
}
if (length(args[['pad_short_edge']])) {
    sp_tree = rkftools::pad_short_edges(sp_tree, threshold=args[['pad_short_edge']], external_only=FALSE)
}
sp_tree = rkftools::force_ultrametric(sp_tree, stop_if_larger_change=0.01)
root_depth = max(node.depth.edgelength(sp_tree))
sp_node_ages = abs(node.depth.edgelength(sp_tree) - root_depth)
sp_node_names = c(sp_tree$tip.label, sp_tree$node.label)
sp_node_table = data.frame(node=sp_node_names, age=sp_node_ages, spp=NA)
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
    if (length(gregexpr('\\(', treetext)[[1]])-length(gregexpr('\\)', treetext)[[1]])==-1) {
        cat('Number of parentheses in the .nhx is not consistent. Trying to fix.')
        treetext <- gsub("\\);", ";", treetext)
    }
    write(treetext, 'tmp.treetext.txt')
    nhxtree = treeio::read.nhx('tmp.treetext.txt')
    file.remove('tmp.treetext.txt')
    return(nhxtree)
}

if (mode=='generax') {
    cat('Reading GeneRax tree.\n')
    nhxtree = read_generax_nhx(generax_file)

    gn_tree = nhxtree@phylo
    
    cols = c('event', 'gn_node', 'gn_node_num', 'lower_sp_node', 'upper_sp_node', 'lower_age', 'upper_age')
    gn_node_table = data.frame(nhxtree@data, stringsAsFactors=FALSE)
    gn_node_table[,'event'] = 'S'
    gn_node_table[is.na(gn_node_table[['D']]),'D'] = 'N'
    gn_node_table[(gn_node_table[['D']]=='Y'),'event'] = 'D'
    colnames(gn_node_table) = sub('^S$', 'lower_sp_node', colnames(gn_node_table))
    gn_node_table[,'upper_sp_node'] = gn_node_table[['lower_sp_node']]
    gn_node_table[,'gn_node'] = c(gn_tree[['tip.label']], gn_tree[['node.label']])
    #gn_node_table = gn_node_table[(gn_node_table[['gn_node']]!='root'),] ### Remove root node
    gn_node_table[(gn_node_table[['event']]=='D'),'upper_sp_node'] = NA
    for (sp_node in unique(gn_node_table[['lower_sp_node']])) {
        node_num = rkftools::get_node_num_by_name(sp_tree, sp_node)
        parent_num = rkftools::get_parent_num(sp_tree, node_num)
        parent_name = rkftools::get_node_name_by_num(sp_tree, parent_num)
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
    gn_node_table[,'gn_node_num'] = rkftools::get_node_num_by_name(gn_tree, gn_node_table[['gn_node']])
    gn_node_table = data.frame(gn_node_table[,cols], stringsAsFactors=FALSE)
}

if (mode=='notung') {
    cat('Reading NOTUNG tree.\n')
    gn_tree = read.tree(gn_file)
    gn_tree$node.label = gsub("\\'", "",gn_tree$node.label)

    gn_node_table = read_notung_parsable(file=parsable_file, mode='D')
    gn_node_table = merge(gn_node_table, data.frame(lower_age=NA, upper_age=NA, spp=NA), all=TRUE)
    if (nrow(gn_node_table) > 0) {
        print(gn_node_table[,'gn_node'])
        gn_node_nums = sapply(gn_node_table[,'gn_node'], function(x){rkftools::get_node_num_by_name(gn_tree, x)})
        print(str(gn_node_nums))
        gn_node_table$gn_node_num = gn_node_nums
        for (i in 1:nrow(gn_node_table)) {
            if (any(sp_node_table$node==gn_node_table$lower_sp_node[i])) {
                gn_node_table$lower_age[i] = sp_node_table$age[sp_node_table$node==gn_node_table$lower_sp_node[i]]
            }
            if (any(sp_node_table$node==gn_node_table$upper_sp_node[i])) {
                gn_node_table$upper_age[i] = sp_node_table$age[sp_node_table$node==gn_node_table$upper_sp_node[i]]
            }
        }
    } else {
        gn_node_table = gn_node_table[0,]
    }

    for (gn_sub in ape::subtrees(gn_tree)) {
        root_node = gn_sub$node.label[1]
        if (! root_node %in% gn_node_table$gn_node) {
            root_num = rkftools::get_node_num_by_name(gn_tree, root_node)
            node_spp = c()
            for (gn_gene in gn_sub$tip.label) {
                pos_underbar = gregexpr("_", gn_gene)[[1]]
                node_sp = substring(gn_gene, 1, pos_underbar[length(pos_underbar)]-1)
                node_spp = c(node_spp, node_sp)
            }
            node_age = min(sp_node_table$age[sapply(sp_node_table$spp, function(x){all(node_spp %in% x)})])
            sp_node = sp_node_table$node[sp_node_table$age==node_age]
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

# Calibration node check
if ((sum(gn_node_table[['event']]=="D") > 0)&(any(is.na(gn_node_table[['upper_age']])))) {
    gn_spp = unique(get_species_names(gn_tree))
    num_sp = length(gn_spp)
    cat('# species in the gene tree:', num_sp, '\n')
    cat('Species in the gene tree:', paste(gn_spp, collapse=', '), '\n')
    num_sp_gntree = max(1, ape::getMRCA(sp_tree, gn_spp))
    if (num_sp_gntree==rkftools::get_root_num(sp_tree)) {
        divtime_max = max_age
        divtime_min = max(ape::node.depth.edgelength(sp_tree))
    } else {
        if (length(gn_spp)==1) {
            num_mrca = rkftools::get_node_num_by_name(sp_tree, gn_spp)
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
root_num = rkftools::get_root_num(gn_tree)
gn_node_table[(gn_node_table$gn_node_num==root_num),'event'] = paste0(gn_node_table[(gn_node_table$gn_node_num==root_num),'event'], '(R)')

droppable_nodes = c()
flag_first = TRUE
for (gn_node_num in gn_node_table[['gn_node_num']]) {
    if (gn_node_num==root_num) {
        next
    }
    child_lower = gn_node_table[(gn_node_table[['gn_node_num']]==gn_node_num),'lower_age']
    child_upper = gn_node_table[(gn_node_table[['gn_node_num']]==gn_node_num),'upper_age']
    ancestor_nums = get_ancestor_num(gn_tree, gn_node_num)
    ancestor_lower = min(gn_node_table[(gn_node_table[['gn_node_num']]%in%ancestor_nums),'lower_age'])
    ancestor_upper = min(gn_node_table[(gn_node_table[['gn_node_num']]%in%ancestor_nums),'upper_age'])
    is_same_constraint = (child_lower>=ancestor_lower) & (child_upper>=ancestor_upper)
    is_same_constraint = ifelse(length(is_same_constraint)==0, FALSE, is_same_constraint)
    if (is_same_constraint) {
        if (flag_first) {
            cat('calibration node removed because of the constraint identical to or greater than one of the upper nodes (name/id/lower/upper):\n')
            flag_first = FALSE
        }
        droppable_name = rkftools::get_node_name_by_num(gn_tree, gn_node_num)
        cat(paste(c(droppable_name, gn_node_num, child_upper, child_lower), collapse='/'), '\n')
        droppable_nodes = c(droppable_nodes, gn_node_num)
    }
}
cat('\n')
gn_node_table_dropped = gn_node_table[(!gn_node_table[['gn_node_num']] %in% droppable_nodes), ]
gn_node_table_dropped = gn_node_table_dropped[(gn_node_table_dropped[,'gn_node_num']>ape::Ntip(gn_tree)),] # Drop leaves
num_constrained_speciation = sum(grepl('^S', gn_node_table_dropped[,'event']))
num_constrained_duplication = sum(grepl('^D', gn_node_table_dropped[,'event']))
cat('Number of constrained speciation nodes:', num_constrained_speciation, '\n')
cat('Number of constrained duplication nodes:', num_constrained_duplication, '\n')

# Calibration table
calibration_table = data.frame(
    node=as.integer(gn_node_table_dropped$gn_node_num),
    age.min=as.numeric(gn_node_table_dropped$lower_age),
    age.max=as.numeric(gn_node_table_dropped$upper_age),
    soft.bounds=NA,
    stringsAsFactors=FALSE
)

calibration_table_R = calibration_table[(calibration_table$node==root_num),]

if ("D" %in% gn_node_table_dropped$event) {
    D_nodes = gn_node_table_dropped[(gn_node_table_dropped$event=='D'),'gn_node_num']
    calibration_table_D = calibration_table[calibration_table$node %in% D_nodes,]
} else {
    calibration_table_D = NA
}

if ("S" %in% gn_node_table_dropped$event) {
    S_nodes = gn_node_table_dropped[(gn_node_table_dropped$event=='S'),'gn_node_num']
    calibration_table_S = calibration_table[calibration_table$node %in% S_nodes,]
} else {
    calibration_table_S = NA
}

calibration_tables = list(
    'RDS' = calibration_table,
    'RS' = rbind(calibration_table_R, calibration_table_S),
    'RD' = rbind(calibration_table_R, calibration_table_D),
    'S' = calibration_table_S,
    'D' = calibration_table_D,
    'R' = calibration_table_R
)

# chronos
chronos_out = NULL
chronos_control = chronos.control()
chronos_control$iter.max = 100000
chronos_control$eval.max = 100000
chronos_control$dual.iter.max = 200

if (all(gn_node_table$lower_age==gn_node_table$upper_age)) {
    # Gene tree without duplication nodes
    calibrated_node = "allS"
    cat("Constrained nodes:", calibrated_node, '\n')
    cat("All nodes are speciation nodes. Transferring node ages from species tree without age inference by chronos.", '\n')
    dup_constraint = NA
    gn_spp = c()
    for (gn_gene in gn_tree$tip.label) {
        pos_underbar = gregexpr("_", gn_gene)[[1]]
        gn_sp = substring(gn_gene, 1, pos_underbar[length(pos_underbar)]-1)
        gn_spp = c(gn_spp, gn_sp)
    }
    drop_spp = sp_tree$tip.label[! sp_tree$tip.label %in% gn_spp]
    if (length(drop_spp) > 0) {
        chronos_out = drop.tip(phy=sp_tree, tip=drop_spp, trim.internal = TRUE)
    } else {
        chronos_out = sp_tree
    }
    gn_tip_index = c()
    for (sp in chronos_out$tip.label) {
        gn_tip_index = c(gn_tip_index, grep(sp, gn_tree$tip.label))
    }
    chronos_out$tip.label = gn_tree$tip.label[gn_tip_index]
    chronos_out = transfer_node_labels(phy_from=gn_tree, phy_to=chronos_out)
    current_calibration_table = calibration_table_S
} else {
    # Gene tree with duplication nodes
    chronos_out = 'PLACEHOLDER'
    class(chronos_out) = 'try-error'
    for (cn in c('RDS','RS','RD','S','D','R')) {
        if ("try-error" %in% class(chronos_out)) {
            calibrated_node = cn
            cat("chronos, calibrated nodes:", calibrated_node, '\n')
            current_calibration_table = calibration_tables[[calibrated_node]]
            chronos_out = try(
                chronos(
                    gn_tree, 
                    lambda=chronos_lambda, 
                    model=chronos_model, 
                    calibration=current_calibration_table, 
                    control=chronos_control
                )
            )
        }
    }
}

save_tree_pdf = function(phy, file, show.age=FALSE, edge_colors=list()) {
    phy = ape::ladderize(phy)
    if (show.age) {
        root_depth = max(node.depth.edgelength(phy))
        node_ages = abs(node.depth.edgelength(phy) - root_depth)
        int_node_ages = node_ages[(length(phy$tip.label)+1):length(node_ages)]
        phy$node.label = paste(phy$node.label, as.character(round(int_node_ages, digits=1)))
    }
    ec = rep('black', nrow(phy[['edge']]))
    if (length(edge_colors)!=0) {
        for (col in names(edge_colors)) {
            ec[(phy[['edge']][,2]%in%edge_colors[[col]])] = col
        }
    }
    pdf(file, height=max(3, length(phy$tip.label)/5+1), width=7.2)
    plot(phy, show.node.label=TRUE, show.tip.label=TRUE, cex=0.5, label.offset=0, edge.color=ec)
    invisible(dev.off())
}

if ("try-error" %in% class(chronos_out)) {
    cat('All attempts for divergence time estimation were failed. Exiting.\n')
    q('no')
} else {
    cat('Writing output files.\n')
    chronos_out2 = chronos_out
    num_neg = 1
    counter = 1
    if (length(args[['pad_short_edge']])) {
        while ((num_neg>0)&(counter<100)) {
            cat(paste0(counter, 'th round of padding started.\n'))
            chronos_out2 = rkftools::pad_short_edges(chronos_out2, threshold=args[['pad_short_edge']], external_only=FALSE)
            chronos_out2 = rkftools::force_ultrametric(chronos_out2, stop_if_larger_change=0.01)
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
    
    ec = list('red'=droppable_nodes, 'blue'=current_calibration_table[['node']])
    save_tree_pdf(phy=gn_tree, file="radte_gene_tree_input.pdf", show.age=FALSE, edge_colors=ec)
    save_tree_pdf(phy=chronos_out2, file="radte_gene_tree_output.pdf", show.age=TRUE, edge_colors=ec)
    save_tree_pdf(phy=sp_tree, file="radte_species_tree.pdf", show.age=TRUE)

    cat('Calibrated nodes:', calibrated_node, '\n')
    cat('Tree height:', max(ape::node.depth.edgelength(sp_tree)), 'million years', '\n')
    cat('Completed: RADTE divergence time estimation', '\n')
}


