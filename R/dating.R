run_radte_dating = function(args, gene_context, species_context, species_parser, calibration_tables) {
    gn_tree = gene_context$tree
    gn_node_table = gene_context$nodes
    sp_tree = species_context$tree
    root_num = get_root_num(gn_tree)
    calibration_table_R = calibration_tables$R
    calibration_table_S = calibration_tables$S
    has_species_node_interval_input = species_context$has_intervals
    has_shared_interval_speciation_input = gn_node_table_has_shared_interval_speciation(gn_node_table)
    dating_backend = args$dating_backend
    random_seed = args$seed
    chronos_iter_max = args$chronos_iter_max
    chronos_eval_max = args$chronos_eval_max
    chronos_dual_iter_max = args$chronos_dual_iter_max
    chronos_high_control_fallback = args$chronos_high_control_fallback
    chronos_lambda = args$chronos_lambda
    chronos_model = args$chronos_model
    allow_constraint_drop = args$allow_constraint_drop
    chronos_attempt_timeout_sec = args$chronos_attempt_timeout_sec
    chronos_total_timeout_sec = args$chronos_total_timeout_sec
    mcmctree_seqfile = args$mcmctree_seqfile
    mcmctree_bin = args$mcmctree_bin
    mcmctree_workdir = args$mcmctree_workdir
    mcmctree_timeout_sec = args$mcmctree_timeout_sec
    mcmctree_usedata = args$mcmctree_usedata
    mcmctree_seqtype = args$mcmctree_seqtype
    mcmctree_clock = args$mcmctree_clock
    mcmctree_model = args$mcmctree_model
    mcmctree_burnin = args$mcmctree_burnin
    mcmctree_sampfreq = args$mcmctree_sampfreq
    mcmctree_nsample = args$mcmctree_nsample
    mcmctree_ncatG = args$mcmctree_ncatG

    # dating
    chronos_out = NULL
    current_calibration_table = calibration_tables[['RS']]
    dating_backend_used = dating_backend
    mcmctree_workdir_used = NA_character_
    mcmctree_executable_used = NA_character_
    mcmctree_banner = NA_character_
    mcmctree_mirror_table = NULL
    shared_speciation_age_summary = empty_shared_speciation_age_summary()
    chronos_model_used = NA_character_
    chronos_lambda_used = NA_real_
    chronos_seed_used = NA_integer_
    has_duplication_event = any(grepl('^D', gn_node_table$event))
    has_interval_calibration = has_species_node_interval_input
    has_shared_interval_speciation = has_shared_interval_speciation_input

    if (dating_backend == 'chronos') {
        if (has_shared_interval_speciation) {
            cat(
                'chronos backend note: mirrored speciation nodes share the same interval bounds, ',
                'but chronos cannot force separate internal nodes to have identical estimated ages unless the bounds are exact.\n',
                sep=''
            )
        }
        chronos_control_profiles = make_chronos_control_profiles(
            iter_max=chronos_iter_max,
            eval_max=chronos_eval_max,
            dual_iter_max=chronos_dual_iter_max,
            enable_high_fallback=chronos_high_control_fallback
        )

        if ((!has_duplication_event) && (!has_interval_calibration)) {
            # Gene tree without duplication nodes
            calibrated_node = "allS"
            cat("Constrained nodes:", calibrated_node, '\n')
            cat("All nodes are speciation nodes. Transferring node ages from species tree without age inference by chronos.", '\n')
            chronos_out = transfer_species_ages(gn_tree, sp_tree, species_parser)
            current_calibration_table = rbind(calibration_table_R, calibration_table_S)
        } else {
            if (!has_duplication_event) {
                cat(
                    'All nodes are speciation nodes, but interval age bounds are present. ',
                    'Running chronos to honor age.min/age.max constraints.\n',
                    sep=''
                )
            }
            validate_tree_edge_lengths(gn_tree, 'gene tree for chronos')
            # Normalize edge lengths to prevent numerical overflow in chronos
            gn_tree = normalize_edge_length_range(gn_tree, max_edge = 1000, min_edge = 1e-8)
            chronos_out = structure('PLACEHOLDER', class='try-error')
            calibrated_node = 'RS'
            current_calibration_table = calibration_tables[['RS']]

            chronos_model_used = chronos_model
            chronos_lambda_used = chronos_lambda
            chronos_seed_used = NA_integer_
            max_restarts_main = 3
            max_restarts_fallback = 2
            alternate_strategies = list(list(model=chronos_model, lambda=chronos_lambda, label='requested'))
            if (chronos_model != 'discrete') {
                alternate_strategies[[length(alternate_strategies)+1]] = list(model='discrete', lambda=chronos_lambda, label='model-discrete')
            }
            if (chronos_model == 'discrete') {
                if (!isTRUE(all.equal(chronos_lambda, 0.1))) {
                    alternate_strategies[[length(alternate_strategies)+1]] = list(model='discrete', lambda=0.1, label='lambda0.1')
                }
                if (!isTRUE(all.equal(chronos_lambda, 0))) {
                    alternate_strategies[[length(alternate_strategies)+1]] = list(model='discrete', lambda=0, label='lambda0')
                }
            }
            if (chronos_model != 'relaxed') {
                alternate_strategies[[length(alternate_strategies)+1]] = list(model='relaxed', lambda=chronos_lambda, label='model-relaxed')
            }
            if (chronos_model != 'correlated') {
                alternate_strategies[[length(alternate_strategies)+1]] = list(model='correlated', lambda=chronos_lambda, label='model-correlated')
            }
            chronos_time_budget = create_chronos_time_budget(chronos_total_timeout_sec)

            calibration_sequence = 'RS'
            if (allow_constraint_drop) {
                calibration_sequence = c('RS', 'S', 'R')
            }
            for (control_profile_name in names(chronos_control_profiles)) {
                if (!("try-error" %in% class(chronos_out))) {
                    break
                }
                chronos_control = chronos_control_profiles[[control_profile_name]]
                cat(
                    '\nchronos control profile: ',
                    control_profile_name,
                    ' (iter.max=',
                    chronos_control$iter.max,
                    ', eval.max=',
                    chronos_control$eval.max,
                    ', dual.iter.max=',
                    chronos_control$dual.iter.max,
                    ')\n',
                    sep=''
                )
                if (control_profile_name == 'high-fallback') {
                    cat('Fast chronos control profile was exhausted; enabling high-cost fallback.\n')
                }
                seed_cursor = random_seed
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
                        alternate_strategies=alternate_strategies,
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
    } else if (dating_backend == 'mcmctree') {
        calibrated_node = 'RS'
        current_calibration_table = calibration_tables[['RS']]
        if (nrow(current_calibration_table) == 0) {
            stop('No root/speciation calibration nodes are available for --dating_backend=mcmctree.')
        }
        if ('allow_constraint_drop' %in% attr(args, 'supplied')) {
            cat('--allow_constraint_drop is ignored when --dating_backend=mcmctree.\n')
        }
        chronos_only_options = c(
            'chronos_attempt_timeout_sec',
            'chronos_total_timeout_sec',
            'chronos_iter_max',
            'chronos_eval_max',
            'chronos_dual_iter_max',
            'chronos_high_control_fallback'
        )
        if (any(chronos_only_options %in% attr(args, 'supplied'))) {
            cat('chronos options are ignored when --dating_backend=mcmctree.\n')
        }
        cat('Running MCMCTree backend.\n')
        mcmctree_out = run_mcmctree_backend(
            phy=gn_tree,
            gn_node_table=gn_node_table,
            calibration_table=current_calibration_table,
            root_num=root_num,
            seqfile=mcmctree_seqfile,
            bin=mcmctree_bin,
            workdir=mcmctree_workdir,
            usedata=mcmctree_usedata,
            seqtype=mcmctree_seqtype,
            clock=mcmctree_clock,
            model=mcmctree_model,
            burnin=mcmctree_burnin,
            sampfreq=mcmctree_sampfreq,
            nsample=mcmctree_nsample,
            ncatG=mcmctree_ncatG,
            seed=random_seed,
            timeout_sec=mcmctree_timeout_sec,
            protected_inputs=unlist(attr(args, "input_files"), use.names=FALSE)
        )
        chronos_out = mcmctree_out$tree
        mcmctree_workdir_used = mcmctree_out$workdir
        mcmctree_mirror_table = mcmctree_out$mirror_table
        shared_speciation_age_summary = mcmctree_out$shared_speciation_ages
        mcmctree_executable_used = mcmctree_out$executable
        mcmctree_banner = mcmctree_out$banner
    }

    if (inherits(chronos_out, 'try-error')) {
        stop('All attempts for divergence time estimation failed: ', as.character(chronos_out))
    }
    list(tree=chronos_out, calibration=current_calibration_table, calibrated_nodes=calibrated_node,
         mirror_table=mcmctree_mirror_table, shared_ages=shared_speciation_age_summary,
         metadata=list(dating_backend=dating_backend_used, calibrated_nodes=calibrated_node,
             seed_requested=random_seed,
             effective_seed=if (dating_backend == 'mcmctree') random_seed else chronos_seed_used,
             chronos_model_used=chronos_model_used, chronos_lambda_used=chronos_lambda_used,
             chronos_seed_used=chronos_seed_used, mcmctree_executable=mcmctree_executable_used,
             mcmctree_banner=mcmctree_banner, mcmctree_workdir=mcmctree_workdir_used))
}
