#!/usr/bin/env Rscript

parse_named_args <- function(values) {
    if (length(values) %% 2L != 0L) {
        stop("Arguments must be supplied as --name value pairs.")
    }
    keys <- values[seq.int(1L, length(values), by = 2L)]
    if (any(!startsWith(keys, "--"))) {
        stop("Every argument name must start with --.")
    }
    keys <- gsub("-", "_", substring(keys, 3L), fixed = TRUE)
    stats::setNames(as.list(values[seq.int(2L, length(values), by = 2L)]), keys)
}

require_arg <- function(args, name) {
    value <- args[[name]]
    if (is.null(value) || !nzchar(value)) {
        stop(sprintf("Missing required argument: --%s", gsub("_", "-", name)))
    }
    value
}

parse_bool <- function(value) {
    normalized <- tolower(trimws(value))
    if (normalized %in% c("true", "1", "yes")) return(TRUE)
    if (normalized %in% c("false", "0", "no")) return(FALSE)
    stop(sprintf("Expected a boolean value, received %s.", value))
}

write_csv <- function(table, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (endsWith(path, ".gz")) {
        connection <- gzfile(path, open = "wt")
        on.exit(close(connection), add = TRUE)
        utils::write.csv(table, connection, row.names = FALSE, na = "")
    } else {
        utils::write.csv(table, path, row.names = FALSE, na = "")
    }
}

safe_label <- function(value) {
    cleaned <- gsub("[^A-Za-z0-9]+", "_", value)
    cleaned <- gsub("^_+|_+$", "", cleaned)
    if (!nzchar(cleaned)) "level" else cleaned
}

design_diagnostics <- function(metadata, formula_text) {
    formula <- stats::as.formula(formula_text)
    design <- tryCatch(
        stats::model.matrix(formula, data = metadata),
        error = function(error) error
    )
    if (inherits(design, "error")) {
        return(list(
            estimable = FALSE,
            reason = paste0("model_matrix_error: ", conditionMessage(design)),
            design = NULL,
            rank = NA_integer_,
            columns = NA_integer_,
            residual_df = NA_integer_
        ))
    }
    rank <- qr(design)$rank
    columns <- ncol(design)
    residual_df <- nrow(design) - rank
    age_column <- which(colnames(design) == "age_decade")
    estimable <- length(age_column) == 1L && rank == columns
    reason <- if (estimable) "" else "rank_deficient_or_age_coefficient_absent"
    list(
        estimable = estimable,
        reason = reason,
        design = design,
        rank = rank,
        columns = columns,
        residual_df = residual_df
    )
}

run_model <- function(counts, metadata, design, age_column) {
    y <- edgeR::DGEList(counts = counts)
    keep <- edgeR::filterByExpr(y, design = design)
    if (!any(keep)) stop("No genes passed design-aware filterByExpr.")
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::normLibSizes(y, method = "TMM")
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    test <- edgeR::glmQLFTest(fit, coef = age_column)
    table <- edgeR::topTags(
        test,
        n = Inf,
        adjust.method = "none",
        sort.by = "none"
    )$table
    list(table = table, n_tested = sum(keep))
}

scenario_table <- function(metadata, primary_formula, high_cell_min, run_lobo) {
    scenarios <- data.frame(
        scenario_id = c("high_cell_support", "drop_sex", "drop_batch"),
        scenario_type = c("support", "covariate", "covariate"),
        model_formula = c(
            primary_formula,
            "~ batch + age_decade",
            "~ sex + age_decade"
        ),
        min_cells = c(high_cell_min, 0L, 0L),
        excluded_batch = c("", "", ""),
        stringsAsFactors = FALSE
    )
    if (run_lobo) {
        batches <- sort(unique(as.character(metadata$batch)))
        labels <- vapply(batches, safe_label, character(1L))
        if (anyDuplicated(labels)) {
            stop("Batch labels are not unique after safe-label normalization.")
        }
        leave_out <- data.frame(
            scenario_id = paste0("leave_one_batch_out__", labels),
            scenario_type = "leave_one_batch_out",
            model_formula = primary_formula,
            min_cells = 0L,
            excluded_batch = batches,
            stringsAsFactors = FALSE
        )
        scenarios <- rbind(scenarios, leave_out)
    }
    scenarios
}

empty_sensitivity_table <- function() {
    data.frame(
        gene_id = character(),
        cell_type = character(),
        analysis_tier = character(),
        scenario_id = character(),
        scenario_type = character(),
        excluded_batch = character(),
        model_formula = character(),
        log2_fc_per_10_years = numeric(),
        average_log_cpm = numeric(),
        ql_f_statistic = numeric(),
        p_value = numeric(),
        n_replicates = integer(),
        age_span_years = numeric(),
        stringsAsFactors = FALSE
    )
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
matrix_path <- require_arg(args, "matrix")
profiles_path <- require_arg(args, "profiles")
genes_path <- require_arg(args, "genes")
primary_results_path <- require_arg(args, "primary_results")
sensitivity_out <- require_arg(args, "sensitivity_out")
diagnostics_out <- require_arg(args, "diagnostics_out")
manifest_out <- require_arg(args, "manifest_out")
primary_formula <- require_arg(args, "primary_formula")
min_replicates <- as.integer(require_arg(args, "min_replicates"))
min_age_span <- as.numeric(require_arg(args, "min_age_span"))
min_residual_df <- as.integer(require_arg(args, "min_residual_df"))
candidate_fdr <- as.numeric(require_arg(args, "candidate_fdr"))
high_cell_min <- as.integer(require_arg(args, "high_cell_min"))
run_lobo <- parse_bool(require_arg(args, "run_leave_one_batch_out"))
expected_r_version <- require_arg(args, "expected_r_version")
expected_bioconductor_version <- require_arg(
    args,
    "expected_bioconductor_version"
)
expected_edger_version <- require_arg(args, "expected_edger_version")

if (!is.finite(candidate_fdr) || candidate_fdr <= 0 || candidate_fdr >= 1) {
    stop("candidate-fdr must be strictly between zero and one.")
}
if (high_cell_min < 1L) stop("high-cell-min must be positive.")

suppressPackageStartupMessages({
    library(Matrix)
    library(edgeR)
    library(BiocManager)
})

observed_versions <- c(
    R = as.character(getRversion()),
    Bioconductor = as.character(BiocManager::version()),
    edgeR = as.character(utils::packageVersion("edgeR"))
)
if (
    !startsWith(observed_versions[["R"]], expected_r_version) ||
    observed_versions[["Bioconductor"]] != expected_bioconductor_version ||
    observed_versions[["edgeR"]] != expected_edger_version
) {
    stop("The R/Bioconductor/edgeR runtime does not match the approved contract.")
}

profiles <- utils::read.csv(profiles_path, stringsAsFactors = FALSE)
genes <- utils::read.csv(genes_path, stringsAsFactors = FALSE)
primary_results <- utils::read.csv(
    primary_results_path,
    stringsAsFactors = FALSE
)
required_profile_columns <- c(
    "profile_id", "biological_replicate_id", "cell_type", "age",
    "age_decade", "sex", "batch", "analysis_tier", "n_cells"
)
required_result_columns <- c(
    "gene_id", "cell_type", "analysis_tier", "log2_fc_per_10_years",
    "fdr_global"
)
missing_profiles <- setdiff(required_profile_columns, colnames(profiles))
missing_results <- setdiff(required_result_columns, colnames(primary_results))
if (length(missing_profiles)) {
    stop("Profile metadata lacks required columns: ", paste(missing_profiles, collapse = ", "))
}
if (length(missing_results)) {
    stop("Primary results lack required columns: ", paste(missing_results, collapse = ", "))
}
if (!identical(colnames(genes), "gene_id")) {
    stop("Gene metadata must contain exactly one gene_id column.")
}
if (anyDuplicated(profiles$profile_id) || anyDuplicated(genes$gene_id)) {
    stop("Profile and gene identifiers must be unique.")
}

matrix_connection <- gzfile(matrix_path, open = "rb")
counts_profiles_by_genes <- Matrix::readMM(matrix_connection)
close(matrix_connection)
if (
    nrow(counts_profiles_by_genes) != nrow(profiles) ||
    ncol(counts_profiles_by_genes) != nrow(genes)
) {
    stop("Matrix dimensions do not match profile and gene metadata.")
}
counts <- methods::as(Matrix::t(counts_profiles_by_genes), "CsparseMatrix")
rownames(counts) <- genes$gene_id
colnames(counts) <- profiles$profile_id

candidates <- primary_results[
    is.finite(primary_results$fdr_global) &
        primary_results$fdr_global < candidate_fdr,
    required_result_columns,
    drop = FALSE
]
if (!nrow(candidates)) stop("No primary-analysis candidates passed global FDR.")
if (anyDuplicated(candidates[c("gene_id", "cell_type")])) {
    stop("Primary candidate gene-by-cell-type keys must be unique.")
}

all_sensitivity <- list()
all_diagnostics <- list()
all_manifest <- list()
cell_types <- unique(profiles$cell_type)

for (cell_type_index in seq_along(cell_types)) {
    cell_type <- cell_types[[cell_type_index]]
    candidate_genes <- candidates$gene_id[candidates$cell_type == cell_type]
    if (!length(candidate_genes)) next

    profile_index <- which(profiles$cell_type == cell_type)
    metadata_base <- profiles[profile_index, , drop = FALSE]
    counts_base <- counts[, profile_index, drop = FALSE]
    scenarios <- scenario_table(
        metadata_base,
        primary_formula,
        high_cell_min,
        run_lobo
    )
    all_manifest[[length(all_manifest) + 1L]] <- scenarios

    for (scenario_index in seq_len(nrow(scenarios))) {
        scenario <- scenarios[scenario_index, , drop = FALSE]
        keep_profiles <- rep(TRUE, nrow(metadata_base))
        if (scenario$scenario_id == "high_cell_support") {
            keep_profiles <- metadata_base$n_cells >= high_cell_min
        } else if (scenario$scenario_type == "leave_one_batch_out") {
            keep_profiles <- metadata_base$batch != scenario$excluded_batch
        }

        metadata <- metadata_base[keep_profiles, , drop = FALSE]
        metadata$sex <- droplevels(factor(metadata$sex))
        metadata$batch <- droplevels(factor(metadata$batch))
        n_profiles <- nrow(metadata)
        age_span <- if (n_profiles) diff(range(metadata$age)) else NA_real_
        reason <- ""
        if (n_profiles < min_replicates) {
            reason <- "insufficient_replicates"
        } else if (!is.finite(age_span) || age_span < min_age_span) {
            reason <- "insufficient_age_span"
        }

        design_info <- design_diagnostics(metadata, scenario$model_formula)
        if (
            !nzchar(reason) &&
            (!design_info$estimable || is.na(design_info$residual_df) ||
                design_info$residual_df < min_residual_df)
        ) {
            reason <- if (nzchar(design_info$reason)) {
                design_info$reason
            } else {
                "insufficient_residual_df"
            }
        }

        if (nzchar(reason)) {
            all_diagnostics[[length(all_diagnostics) + 1L]] <- data.frame(
                cell_type = cell_type,
                analysis_tier = metadata_base$analysis_tier[[1L]],
                scenario_id = scenario$scenario_id,
                scenario_type = scenario$scenario_type,
                excluded_batch = scenario$excluded_batch,
                model_formula = scenario$model_formula,
                n_replicates = n_profiles,
                age_span_years = age_span,
                design_columns = design_info$columns,
                design_rank = design_info$rank,
                residual_df = design_info$residual_df,
                n_genes_tested = 0L,
                n_candidates_expected = length(candidate_genes),
                n_candidates_observed = 0L,
                status = "not_estimable",
                reason = reason,
                stringsAsFactors = FALSE
            )
            next
        }

        counts_scenario <- as.matrix(counts_base[, keep_profiles, drop = FALSE])
        age_column <- which(colnames(design_info$design) == "age_decade")
        model <- tryCatch(
            run_model(
                counts_scenario,
                metadata,
                design_info$design,
                age_column
            ),
            error = function(error) error
        )
        rm(counts_scenario)
        gc(verbose = FALSE)

        if (inherits(model, "error")) {
            all_diagnostics[[length(all_diagnostics) + 1L]] <- data.frame(
                cell_type = cell_type,
                analysis_tier = metadata_base$analysis_tier[[1L]],
                scenario_id = scenario$scenario_id,
                scenario_type = scenario$scenario_type,
                excluded_batch = scenario$excluded_batch,
                model_formula = scenario$model_formula,
                n_replicates = n_profiles,
                age_span_years = age_span,
                design_columns = design_info$columns,
                design_rank = design_info$rank,
                residual_df = design_info$residual_df,
                n_genes_tested = 0L,
                n_candidates_expected = length(candidate_genes),
                n_candidates_observed = 0L,
                status = "model_failed",
                reason = conditionMessage(model),
                stringsAsFactors = FALSE
            )
            next
        }

        model_table <- model$table
        model_table$gene_id <- rownames(model_table)
        selected <- model_table[model_table$gene_id %in% candidate_genes, , drop = FALSE]
        result <- data.frame(
            gene_id = selected$gene_id,
            cell_type = cell_type,
            analysis_tier = metadata_base$analysis_tier[[1L]],
            scenario_id = scenario$scenario_id,
            scenario_type = scenario$scenario_type,
            excluded_batch = scenario$excluded_batch,
            model_formula = scenario$model_formula,
            log2_fc_per_10_years = selected$logFC,
            average_log_cpm = selected$logCPM,
            ql_f_statistic = selected$F,
            p_value = selected$PValue,
            n_replicates = n_profiles,
            age_span_years = age_span,
            stringsAsFactors = FALSE
        )
        all_sensitivity[[length(all_sensitivity) + 1L]] <- result
        all_diagnostics[[length(all_diagnostics) + 1L]] <- data.frame(
            cell_type = cell_type,
            analysis_tier = metadata_base$analysis_tier[[1L]],
            scenario_id = scenario$scenario_id,
            scenario_type = scenario$scenario_type,
            excluded_batch = scenario$excluded_batch,
            model_formula = scenario$model_formula,
            n_replicates = n_profiles,
            age_span_years = age_span,
            design_columns = design_info$columns,
            design_rank = design_info$rank,
            residual_df = design_info$residual_df,
            n_genes_tested = model$n_tested,
            n_candidates_expected = length(candidate_genes),
            n_candidates_observed = nrow(result),
            status = "completed",
            reason = "",
            stringsAsFactors = FALSE
        )
        message(sprintf(
            "[pseudobulk_robustness] %s | %s | n=%d | candidates=%d/%d",
            cell_type,
            scenario$scenario_id,
            n_profiles,
            nrow(result),
            length(candidate_genes)
        ))
    }
}

sensitivity_table <- if (length(all_sensitivity)) {
    do.call(rbind, all_sensitivity)
} else {
    empty_sensitivity_table()
}
diagnostics_table <- do.call(rbind, all_diagnostics)
manifest_table <- unique(do.call(rbind, all_manifest))

sensitivity_table <- sensitivity_table[order(
    sensitivity_table$cell_type,
    sensitivity_table$scenario_id,
    sensitivity_table$p_value,
    sensitivity_table$gene_id
), ]
diagnostics_table <- diagnostics_table[order(
    diagnostics_table$cell_type,
    diagnostics_table$scenario_id
), ]
manifest_table <- manifest_table[order(manifest_table$scenario_id), ]

write_csv(sensitivity_table, sensitivity_out)
write_csv(diagnostics_table, diagnostics_out)
write_csv(manifest_table, manifest_out)
message(sprintf(
    "[pseudobulk_robustness] completed %d sensitivity rows across %d model scenarios",
    nrow(sensitivity_table),
    sum(diagnostics_table$status == "completed")
))
