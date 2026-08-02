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
    if (!nzchar(cleaned)) "cell_type" else cleaned
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
    estimable <- (
        length(age_column) == 1L &&
        rank == columns &&
        residual_df >= 0L
    )
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
    if (!any(keep)) {
        stop("No genes passed design-aware filterByExpr.")
    }
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::normLibSizes(y, method = "TMM")
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    test <- edgeR::glmQLFTest(fit, coef = age_column)
    table <- edgeR::topTags(
        test,
        n = Inf,
        adjust.method = "BH",
        sort.by = "none"
    )$table
    list(y = y, test = test, table = table, n_tested = sum(keep))
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
matrix_path <- require_arg(args, "matrix")
profiles_path <- require_arg(args, "profiles")
genes_path <- require_arg(args, "genes")
combined_out <- require_arg(args, "combined_out")
celltype_dir <- require_arg(args, "celltype_dir")
manifest_out <- require_arg(args, "manifest_out")
diagnostics_out <- require_arg(args, "diagnostics_out")
plot_out <- require_arg(args, "plot_out")
runtime_versions_out <- require_arg(args, "runtime_versions_out")
session_info_out <- require_arg(args, "session_info_out")
primary_formula <- require_arg(args, "primary_formula")
min_replicates <- as.integer(require_arg(args, "min_replicates"))
min_age_span <- as.numeric(require_arg(args, "min_age_span"))
min_residual_df <- as.integer(require_arg(args, "min_residual_df"))
expected_r_version <- require_arg(args, "expected_r_version")
expected_bioconductor_version <- require_arg(
    args,
    "expected_bioconductor_version"
)
expected_edger_version <- require_arg(args, "expected_edger_version")

suppressPackageStartupMessages({
    library(Matrix)
    library(edgeR)
    library(BiocManager)
})

runtime_versions <- data.frame(
    component = c("R", "Bioconductor", "edgeR", "Matrix"),
    expected = c(
        expected_r_version,
        expected_bioconductor_version,
        expected_edger_version,
        ""
    ),
    observed = c(
        as.character(getRversion()),
        as.character(BiocManager::version()),
        as.character(utils::packageVersion("edgeR")),
        as.character(utils::packageVersion("Matrix"))
    ),
    stringsAsFactors = FALSE
)
write_csv(runtime_versions, runtime_versions_out)

version_matches <- c(
    startsWith(as.character(getRversion()), expected_r_version),
    as.character(BiocManager::version()) == expected_bioconductor_version,
    as.character(utils::packageVersion("edgeR")) == expected_edger_version
)
if (!all(version_matches)) {
    stop(
        "The R/Bioconductor/edgeR runtime does not match the approved ",
        "version contract; inspect ", runtime_versions_out, "."
    )
}

profiles <- utils::read.csv(
    profiles_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
genes <- utils::read.csv(
    genes_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
required_profile_columns <- c(
    "profile_id",
    "biological_replicate_id",
    "cell_type",
    "age",
    "age_decade",
    "sex",
    "batch",
    "analysis_tier",
    "n_cells",
    "library_size"
)
missing_profile_columns <- setdiff(required_profile_columns, colnames(profiles))
if (length(missing_profile_columns)) {
    stop(
        "Profile metadata lacks required columns: ",
        paste(missing_profile_columns, collapse = ", ")
    )
}
if (!identical(colnames(genes), "gene_id")) {
    stop("Gene metadata must contain exactly one gene_id column.")
}
if (anyDuplicated(profiles$profile_id)) {
    stop("profile_id values must be unique.")
}
if (anyDuplicated(genes$gene_id)) {
    stop("gene_id values must be unique.")
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
counts <- methods::as(
    Matrix::t(counts_profiles_by_genes),
    "CsparseMatrix"
)
rownames(counts) <- genes$gene_id
colnames(counts) <- profiles$profile_id
observed_library_sizes <- as.numeric(Matrix::colSums(counts))
if (!isTRUE(all.equal(
    unname(observed_library_sizes),
    as.numeric(profiles$library_size),
    tolerance = 0
))) {
    stop("Profile library sizes do not match matrix column sums.")
}

dir.create(dirname(combined_out), recursive = TRUE, showWarnings = FALSE)
dir.create(celltype_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_out), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(plot_out, width = 8, height = 6, onefile = TRUE)

results <- list()
diagnostics <- list()
cell_types <- unique(profiles$cell_type)
for (index in seq_along(cell_types)) {
    cell_type <- cell_types[[index]]
    selected <- which(profiles$cell_type == cell_type)
    metadata <- profiles[selected, , drop = FALSE]
    metadata$sex <- droplevels(factor(metadata$sex))
    metadata$batch <- droplevels(factor(metadata$batch))
    counts_cell_type <- as.matrix(counts[, selected, drop = FALSE])
    n_profiles <- nrow(metadata)
    age_span <- diff(range(metadata$age))

    support_reason <- ""
    if (n_profiles < min_replicates) {
        support_reason <- "insufficient_replicates"
    } else if (!is.finite(age_span) || age_span < min_age_span) {
        support_reason <- "insufficient_age_span"
    }
    if (nzchar(support_reason)) {
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
            cell_type = cell_type,
            analysis_tier = metadata$analysis_tier[[1L]],
            model_scope = "not_tested",
            n_replicates = n_profiles,
            age_span_years = age_span,
            design_columns = NA_integer_,
            design_rank = NA_integer_,
            residual_df = NA_integer_,
            n_genes_total = nrow(counts_cell_type),
            n_genes_tested = 0L,
            status = "excluded",
            reason = support_reason,
            stringsAsFactors = FALSE
        )
        next
    }

    primary <- design_diagnostics(metadata, primary_formula)
    model_scope <- "adjusted_primary"
    model_formula <- primary_formula
    selected_design <- primary
    if (
        !primary$estimable ||
        is.na(primary$residual_df) ||
        primary$residual_df < min_residual_df
    ) {
        model_scope <- "age_only_sensitivity"
        model_formula <- "~ age_decade"
        selected_design <- design_diagnostics(metadata, model_formula)
    }

    if (
        !selected_design$estimable ||
        is.na(selected_design$residual_df) ||
        selected_design$residual_df < min_residual_df
    ) {
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
            cell_type = cell_type,
            analysis_tier = metadata$analysis_tier[[1L]],
            model_scope = model_scope,
            n_replicates = n_profiles,
            age_span_years = age_span,
            design_columns = selected_design$columns,
            design_rank = selected_design$rank,
            residual_df = selected_design$residual_df,
            n_genes_total = nrow(counts_cell_type),
            n_genes_tested = 0L,
            status = "excluded",
            reason = if (
                nzchar(selected_design$reason)
            ) selected_design$reason else "insufficient_residual_df",
            stringsAsFactors = FALSE
        )
        next
    }

    age_column <- which(
        colnames(selected_design$design) == "age_decade"
    )
    model <- tryCatch(
        run_model(
            counts_cell_type,
            metadata,
            selected_design$design,
            age_column
        ),
        error = function(error) error
    )
    if (inherits(model, "error")) {
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
            cell_type = cell_type,
            analysis_tier = metadata$analysis_tier[[1L]],
            model_scope = model_scope,
            n_replicates = n_profiles,
            age_span_years = age_span,
            design_columns = selected_design$columns,
            design_rank = selected_design$rank,
            residual_df = selected_design$residual_df,
            n_genes_total = nrow(counts_cell_type),
            n_genes_tested = 0L,
            status = "model_failed",
            reason = conditionMessage(model),
            stringsAsFactors = FALSE
        )
        next
    }

    table <- model$table
    table$gene_id <- rownames(table)
    result <- data.frame(
        gene_id = table$gene_id,
        cell_type = cell_type,
        analysis_tier = metadata$analysis_tier[[1L]],
        model_scope = model_scope,
        model_formula = model_formula,
        log2_fc_per_10_years = table$logFC,
        average_log_cpm = table$logCPM,
        ql_f_statistic = table$F,
        p_value = table$PValue,
        fdr_within_celltype = table$FDR,
        n_replicates = n_profiles,
        age_span_years = age_span,
        stringsAsFactors = FALSE
    )
    results[[length(results) + 1L]] <- result
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        cell_type = cell_type,
        analysis_tier = metadata$analysis_tier[[1L]],
        model_scope = model_scope,
        n_replicates = n_profiles,
        age_span_years = age_span,
        design_columns = selected_design$columns,
        design_rank = selected_design$rank,
        residual_df = selected_design$residual_df,
        n_genes_total = nrow(counts_cell_type),
        n_genes_tested = model$n_tested,
        status = "completed",
        reason = "",
        stringsAsFactors = FALSE
    )

    limma::plotMD(
        model$test,
        main = sprintf(
            "%s (%s; n=%d)",
            cell_type,
            metadata$analysis_tier[[1L]],
            n_profiles
        )
    )
    graphics::abline(h = 0, col = "grey40", lty = 2)
}
grDevices::dev.off()

diagnostics_table <- do.call(rbind, diagnostics)
write_csv(diagnostics_table, diagnostics_out)
if (!length(results)) {
    stop("No cell type produced an edgeR result table.")
}

combined <- do.call(rbind, results)
combined$fdr_global <- stats::p.adjust(combined$p_value, method = "BH")
combined <- combined[
    order(
        combined$analysis_tier,
        combined$cell_type,
        combined$p_value,
        combined$gene_id
    ),
]
write_csv(combined, combined_out)

manifest <- list()
for (index in seq_along(cell_types)) {
    cell_type <- cell_types[[index]]
    cell_result <- combined[combined$cell_type == cell_type, , drop = FALSE]
    if (!nrow(cell_result)) {
        next
    }
    output_filename <- sprintf(
        "%02d_%s.csv.gz",
        index,
        safe_label(cell_type)
    )
    output_path <- file.path(celltype_dir, output_filename)
    write_csv(cell_result, output_path)
    manifest[[length(manifest) + 1L]] <- data.frame(
        cell_type = cell_type,
        analysis_tier = cell_result$analysis_tier[[1L]],
        model_scope = cell_result$model_scope[[1L]],
        n_genes_tested = nrow(cell_result),
        result_path = file.path("by_cell_type", output_filename),
        stringsAsFactors = FALSE
    )
}
write_csv(do.call(rbind, manifest), manifest_out)

dir.create(dirname(session_info_out), recursive = TRUE, showWarnings = FALSE)
session_connection <- file(session_info_out, open = "wt")
sink(session_connection)
print(sessionInfo())
sink()
close(session_connection)

message(
    sprintf(
        "[pseudobulk_edger] completed %d cell types and %d gene-level tests",
        length(unique(combined$cell_type)),
        nrow(combined)
    )
)
