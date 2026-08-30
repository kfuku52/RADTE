save_tree_pdf = function(phy, file, show.age=FALSE, edge_colors=list()) {
    phy = ape::ladderize(phy)
    if (show.age) {
        root_depth = max(ape::node.depth.edgelength(phy))
        node_ages = abs(ape::node.depth.edgelength(phy) - root_depth)
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
    grDevices::pdf(file, height=max(3, length(phy$tip.label)/5+1), width=7.2)
    pdf_device = grDevices::dev.cur()
    on.exit({
        open_devices = grDevices::dev.list()
        if (!is.null(open_devices) && pdf_device %in% open_devices) {
            grDevices::dev.off(pdf_device)
        }
    }, add=TRUE)
    graphics::plot(phy, show.node.label=FALSE, show.tip.label=TRUE, cex=0.5, label.offset=0,
         edge.color=ec2, root.edge=TRUE)
    ape::nodelabels(text=phy[['node.label']], col=node_colors, bg='white', cex=0.5)
    invisible(grDevices::dev.off(pdf_device))
}

atomic_write_file = function(file, writer) {
    parent_dir = dirname(file)
    dir.create(parent_dir, recursive=TRUE, showWarnings=FALSE)
    tmp_file = tempfile(pattern=paste0('.', basename(file), '-'), tmpdir=parent_dir)
    on.exit(unlink(tmp_file), add=TRUE)
    writer(tmp_file)
    if (!file.exists(tmp_file)) {
        stop('Atomic writer did not create its temporary output: ', tmp_file)
    }
    backup_file = NULL
    if (dir.exists(file)) stop('Output path is an existing directory: ', file)
    if (file.exists(file) || (!is.na(Sys.readlink(file)) && nzchar(Sys.readlink(file)))) {
        backup_file = tempfile(pattern=paste0('.', basename(file), '-backup-'), tmpdir=parent_dir)
        if (!file.rename(file, backup_file)) {
            stop('Failed to stage the existing output for atomic replacement: ', file)
        }
    }
    if (!file.rename(tmp_file, file)) {
        if (!is.null(backup_file) && file.exists(backup_file)) {
            file.rename(backup_file, file)
        }
        stop('Failed to atomically replace output file: ', file)
    }
    if (!is.null(backup_file)) {
        unlink(backup_file)
    }
    invisible(file)
}

atomic_write_lines = function(lines, file) {
    atomic_write_file(file, function(tmp) writeLines(lines, con=tmp, useBytes=TRUE))
}

atomic_write_table = function(value, file, ...) {
    atomic_write_file(file, function(tmp) utils::write.table(value, file=tmp, ...))
}

atomic_write_tree = function(tree, file) {
    atomic_write_file(file, function(tmp) ape::write.tree(tree, file=tmp))
}

atomic_save_tree_pdf = function(phy, file, show.age=FALSE, edge_colors=list()) {
    atomic_write_file(
        file,
        function(tmp) save_tree_pdf(phy=phy, file=tmp, show.age=show.age, edge_colors=edge_colors)
    )
}

atomic_copy_file = function(source, target) {
    atomic_write_file(target, function(tmp) {
        if (!file.copy(source, tmp, overwrite=TRUE)) {
            stop('Failed to copy output artifact: ', source)
        }
    })
}

compute_file_sha256 = function(file) {
    if (!is_nonempty_scalar_string(file) || !file.exists(file)) {
        return(NA_character_)
    }
    file_abs = normalizePath(file, mustWork=TRUE)
    sha256sum_bin = Sys.which('sha256sum')
    if (nzchar(sha256sum_bin)) {
        out = suppressWarnings(system2(sha256sum_bin, shQuote(file_abs), stdout=TRUE, stderr=FALSE))
        if (length(out) > 0 && grepl('^[0-9a-fA-F]{64}', out[[1]])) {
            return(tolower(sub('^([0-9a-fA-F]{64}).*$', '\\1', out[[1]])))
        }
    }
    shasum_bin = Sys.which('shasum')
    if (nzchar(shasum_bin)) {
        out = suppressWarnings(system2(shasum_bin, c('-a', '256', shQuote(file_abs)), stdout=TRUE, stderr=FALSE))
        if (length(out) > 0 && grepl('^[0-9a-fA-F]{64}', out[[1]])) {
            return(tolower(sub('^([0-9a-fA-F]{64}).*$', '\\1', out[[1]])))
        }
    }
    certutil_bin = Sys.which('certutil')
    if (nzchar(certutil_bin)) {
        out = suppressWarnings(system2(certutil_bin, c('-hashfile', shQuote(file_abs), 'SHA256'), stdout=TRUE, stderr=FALSE))
        candidates = gsub('[^0-9a-fA-F]', '', out)
        candidates = candidates[nchar(candidates) == 64]
        if (length(candidates) > 0) {
            return(tolower(candidates[[1]]))
        }
    }
    NA_character_
}

write_run_manifest = function(entries, file) {
    manifest = data.frame(
        key=names(entries),
        value=vapply(entries, function(value) {
            if (length(value) == 0L || all(is.na(value))) return('')
            paste(as.character(value), collapse=',')
        }, character(1)),
        stringsAsFactors=FALSE
    )
    atomic_write_table(manifest, file=file, sep='\t', quote=FALSE, row.names=FALSE, na='')
}

radte_output_suffixes = function() {
    c('_calibrated_nodes.txt', '_gene_tree_output.nwk', '_calibration_used.tsv',
      '_calibration_all.tsv', '_gene_tree.tsv', '_species_tree.tsv',
      '_shared_speciation_ages.tsv', '_gene_tree_input.pdf', '_gene_tree_output.pdf',
      '_species_tree.pdf', '_run_manifest.tsv',
      paste0('_mcmctree_', c('mcmctree.ctl', 'input.trees', 'out.txt', 'mcmc.txt',
                            'FigTree.tre', 'mcmctree.stdout.log', 'mcmctree.stderr.log')))
}

begin_output_transaction = function(outdir, prefix, inputs=character()) {
    transaction = new.env(parent=emptyenv())
    transaction$outdir = outdir
    transaction$prefix = prefix
    transaction$inputs = inputs
    transaction$input_hashes = vapply(inputs, compute_file_sha256, character(1))
    transaction$targets = file.path(outdir, paste0(prefix, radte_output_suffixes()))
    assert_no_input_output_collision(inputs, transaction$targets)
    if (any(dir.exists(transaction$targets))) stop('Output path is an existing directory.')
    transaction$lock = file.path(outdir, paste0('.', prefix, '.radte-lock'))
    if (!dir.create(transaction$lock, showWarnings=FALSE)) {
        stop('Another run owns this output prefix, or a stale lock exists: ', transaction$lock)
    }
    transaction$stage = tempfile(pattern=paste0('.', prefix, '-run-'), tmpdir=outdir)
    if (!dir.create(transaction$stage)) {
        unlink(transaction$lock, recursive=TRUE)
        stop('Could not create output staging directory.')
    }
    transaction$run_id = basename(transaction$stage)
    transaction$preserve = FALSE
    writeLines(c(paste0('pid=', Sys.getpid()), paste0('run_id=', transaction$run_id)),
               file.path(transaction$lock, 'owner.txt'))
    transaction
}

cleanup_output_transaction = function(transaction) {
    if (!isTRUE(transaction$preserve)) unlink(transaction$stage, recursive=TRUE)
    unlink(transaction$lock, recursive=TRUE)
    invisible(NULL)
}

commit_output_transaction = function(transaction, manifest_name) {
    sources = list.files(transaction$stage, full.names=TRUE)
    sources = c(sources[basename(sources) != manifest_name], file.path(transaction$stage, manifest_name))
    if (!all(file.exists(sources))) stop('Output manifest or staged artifact is missing.')
    targets = file.path(transaction$outdir, basename(sources))
    if (!all(targets %in% transaction$targets)) stop('Unexpected output filename in staging directory.')
    assert_no_input_output_collision(transaction$inputs, transaction$targets)
    current_hashes = vapply(transaction$inputs, compute_file_sha256, character(1))
    if (!identical(current_hashes, transaction$input_hashes)) stop('Input files changed during the run; outputs were not published.')
    backup_dir = file.path(transaction$stage, '.backup')
    dir.create(backup_dir)
    backed_up = character()
    published = character()
    complete = FALSE
    on.exit({
        if (!complete) {
            unlink(published)
            restored = vapply(backed_up, function(target) {
                file.rename(file.path(backup_dir, basename(target)), target)
            }, logical(1))
            if (!all(restored)) {
                transaction$preserve = TRUE
                warning('Output rollback was incomplete; recover previous artifacts from: ', backup_dir)
            }
        }
    }, add=TRUE)
    # Preserve the previous complete generation until every writer (including
    # PDF generation) succeeds. Manifest is always the final published file.
    suspendInterrupts({
        for (target in transaction$targets) {
            if (file.exists(target) || (!is.na(Sys.readlink(target)) && nzchar(Sys.readlink(target)))) {
                if (dir.exists(target) || !file.rename(target, file.path(backup_dir, basename(target)))) {
                    stop('Could not back up output: ', target)
                }
                backed_up = c(backed_up, target)
            }
        }
        for (i in seq_along(sources)) {
            if (!file.rename(sources[[i]], targets[[i]])) stop('Could not publish output: ', targets[[i]])
            published = c(published, targets[[i]])
        }
        complete = TRUE
    })
    invisible(targets)
}
