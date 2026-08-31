#!/usr/bin/env Rscript

parse_named_args <- function(values) {
    if (length(values) %% 2L != 0L) stop("Arguments must be --name value pairs.")
    keys <- values[seq.int(1L, length(values), by = 2L)]
    if (any(!startsWith(keys, "--"))) stop("Argument names must start with --.")
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

design_diagnostics <- function(metadata, formula_text, require_repeated = TRUE) {
    formula <- tryCatch(stats::as.formula(formula_text), error = function(error) error)
    if (inherits(formula, "error")) {
        return(list(estimable = FALSE, reason = conditionMessage(formula)))
    }
    fixed_formula <- reformulas::nobars(formula)
    design <- tryCatch(
        stats::model.matrix(fixed_formula, data = metadata),
        error = function(error) error
    )
    if (inherits(design, "error")) {
        return(list(estimable = FALSE, reason = conditionMessage(design)))
    }
    rank <- qr(design)$rank
    columns <- ncol(design)
    age_column <- which(colnames(design) == "age_decade")
    residual_df <- nrow(design) - rank
    subject_counts <- table(metadata$subject_id)
    repeated_support <- !require_repeated || any(subject_counts > 1L)
    estimable <- (
        length(age_column) == 1L && rank == columns && residual_df >= 0L &&
        repeated_support
    )
    list(
        estimable = estimable,
        reason = if (estimable) "" else "rank_deficient_age_absent_or_no_repeated_subjects",
        formula = formula,
        fixed_design = design,
        rank = rank,
        columns = columns,
        residual_df = residual_df,
        n_subjects = length(subject_counts),
        n_subjects_repeated = sum(subject_counts > 1L),
        max_samples_per_subject = max(subject_counts)
    )
}

validate_result_table <- function(table, context, require_z = FALSE) {
    required <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
    if (require_z) required <- c(required, "z.std")
    missing <- setdiff(required, colnames(table))
    if (length(missing)) {
        stop(context, " result table lacks columns: ", paste(missing, collapse = ", "))
    }
    nonfinite <- required[vapply(
        table[required], function(values) any(!is.finite(values)), logical(1L)
    )]
    if (length(nonfinite)) {
        stop(
            context, " produced non-finite inferential values in: ",
            paste(nonfinite, collapse = ", ")
        )
    }
    if (any(table$P.Value < 0 | table$P.Value > 1)) {
        stop(context, " produced p-values outside [0, 1].")
    }
    invisible(table)
}

label_diagnostic_page <- function(cell_type) {
    graphics::mtext(
        paste("Population:", cell_type),
        side = 3,
        line = -1.1,
        adj = 0.98,
        cex = 0.72,
        col = "#444444"
    )
}

run_model <- function(counts, metadata, formula, fixed_design, cell_type) {
    y <- edgeR::DGEList(counts = counts)
    keep <- edgeR::filterByExpr(y, design = fixed_design)
    if (!any(keep)) stop("No genes passed design-aware filterByExpr.")
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::normLibSizes(y, method = "TMM")
    parameter <- BiocParallel::SerialParam()
    vobj <- variancePartition::voomWithDreamWeights(
        y, formula, metadata, plot = TRUE, BPPARAM = parameter
    )
    label_diagnostic_page(cell_type)
    fit <- variancePartition::dream(
        vobj, formula, metadata, BPPARAM = parameter
    )
    fit <- variancePartition::eBayes(fit)
    if (!("age_decade" %in% colnames(fit))) {
        stop("The fitted dream model does not expose age_decade.")
    }
    table <- variancePartition::topTable(
        fit, coef = "age_decade", number = Inf, sort.by = "none"
    )
    validate_result_table(table, "dream", require_z = TRUE)
    list(y = y, fit = fit, table = table, n_tested = sum(keep))
}

run_voom_diagnostic <- function(counts, metadata, formula, fixed_design, cell_type) {
    y <- edgeR::DGEList(counts = counts)
    keep <- edgeR::filterByExpr(y, design = fixed_design)
    if (!any(keep)) stop("No genes passed design-aware filterByExpr.")
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::normLibSizes(y, method = "TMM")
    variancePartition::voomWithDreamWeights(
        y,
        formula,
        metadata,
        plot = TRUE,
        BPPARAM = BiocParallel::SerialParam()
    )
    label_diagnostic_page(cell_type)
    invisible(sum(keep))
}

run_fixed_model <- function(counts, metadata, formula, fixed_design) {
    y <- edgeR::DGEList(counts = counts)
    keep <- edgeR::filterByExpr(y, design = fixed_design)
    if (!any(keep)) stop("No genes passed design-aware filterByExpr.")
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::normLibSizes(y, method = "TMM")
    vobj <- limma::voom(y, design = fixed_design, plot = FALSE)
    fit <- limma::lmFit(vobj, design = fixed_design)
    fit <- limma::eBayes(fit)
    if (!("age_decade" %in% colnames(fit))) {
        stop("The fixed sensitivity model does not expose age_decade.")
    }
    table <- limma::topTable(
        fit, coef = "age_decade", number = Inf, sort.by = "none"
    )
    validate_result_table(table, "fixed sensitivity")
    list(y = y, fit = fit, table = table, n_tested = sum(keep))
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
matrix_path <- require_arg(args, "matrix")
profiles_path <- require_arg(args, "profiles")
genes_path <- require_arg(args, "genes")
combined_out <- require_arg(args, "combined_out")
sensitivity_out <- require_arg(args, "sensitivity_out")
celltype_dir <- require_arg(args, "celltype_dir")
manifest_out <- require_arg(args, "manifest_out")
diagnostics_out <- require_arg(args, "diagnostics_out")
plot_out <- require_arg(args, "plot_out")
runtime_versions_out <- require_arg(args, "runtime_versions_out")
session_info_out <- require_arg(args, "session_info_out")
primary_formula <- require_arg(args, "primary_formula")
sensitivity_sample_rule <- require_arg(args, "sensitivity_sample_rule")
min_sample_units <- as.integer(require_arg(args, "min_sample_units"))
min_subjects <- as.integer(require_arg(args, "min_subjects"))
min_age_span <- as.numeric(require_arg(args, "min_age_span"))
min_residual_df <- as.integer(require_arg(args, "min_residual_df"))
expected_r_version <- require_arg(args, "expected_r_version")
expected_bioconductor_version <- require_arg(args, "expected_bioconductor_version")
expected_edger_version <- require_arg(args, "expected_edger_version")
expected_variance_partition_version <- require_arg(
    args, "expected_variance_partition_version"
)
diagnostics_only <- identical(tolower(args[["diagnostics_only"]]), "true")
if (sensitivity_sample_rule != "earliest_age_then_sample_id") {
    stop("Unsupported one-sample-per-subject selection rule: ", sensitivity_sample_rule)
}

suppressPackageStartupMessages({
    library(BiocManager)
    library(BiocParallel)
    library(edgeR)
    library(limma)
    library(Matrix)
    library(variancePartition)
})

runtime_versions <- data.frame(
    component = c("R", "Bioconductor", "edgeR", "variancePartition", "limma", "Matrix"),
    expected = c(
        expected_r_version, expected_bioconductor_version,
        expected_edger_version, expected_variance_partition_version, "", ""
    ),
    observed = c(
        as.character(getRversion()), as.character(BiocManager::version()),
        as.character(utils::packageVersion("edgeR")),
        as.character(utils::packageVersion("variancePartition")),
        as.character(utils::packageVersion("limma")),
        as.character(utils::packageVersion("Matrix"))
    ),
    stringsAsFactors = FALSE
)
write_csv(runtime_versions, runtime_versions_out)
if (!all(c(
    startsWith(as.character(getRversion()), expected_r_version),
    as.character(BiocManager::version()) == expected_bioconductor_version,
    as.character(utils::packageVersion("edgeR")) == expected_edger_version,
    as.character(utils::packageVersion("variancePartition")) ==
        expected_variance_partition_version
))) {
    stop("The repeated-measures R runtime does not match the approved version contract.")
}

profiles <- utils::read.csv(profiles_path, stringsAsFactors = FALSE, check.names = FALSE)
genes <- utils::read.csv(genes_path, stringsAsFactors = FALSE, check.names = FALSE)
required_profile_columns <- c(
    "profile_id", "sample_unit_id", "subject_id", "cell_type", "age",
    "age_decade", "sex", "batch", "analysis_tier", "n_cells", "library_size"
)
missing_profile_columns <- setdiff(required_profile_columns, colnames(profiles))
if (length(missing_profile_columns)) {
    stop("Profile metadata lacks required columns: ", paste(missing_profile_columns, collapse = ", "))
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
) stop("Matrix dimensions do not match profile and gene metadata.")
counts <- methods::as(Matrix::t(counts_profiles_by_genes), "CsparseMatrix")
rownames(counts) <- genes$gene_id
colnames(counts) <- profiles$profile_id
if (!isTRUE(all.equal(
    unname(as.numeric(Matrix::colSums(counts))),
    as.numeric(profiles$library_size), tolerance = 0
))) stop("Profile library sizes do not match matrix column sums.")

dir.create(dirname(combined_out), recursive = TRUE, showWarnings = FALSE)
dir.create(celltype_dir, recursive = TRUE, showWarnings = FALSE)
plot_temp <- tempfile(pattern = "dream_diagnostics_", fileext = ".pdf")
grDevices::pdf(plot_temp, width = 8, height = 6, onefile = TRUE)
results <- list()
sensitivity_results <- list()
diagnostics <- list()
diagnostic_pages <- 0L
cell_types <- unique(profiles$cell_type)
for (index in seq_along(cell_types)) {
    cell_type <- cell_types[[index]]
    selected <- which(profiles$cell_type == cell_type)
    metadata <- profiles[selected, , drop = FALSE]
    metadata$sex <- droplevels(factor(metadata$sex))
    metadata$batch <- droplevels(factor(metadata$batch))
    metadata$subject_id <- droplevels(factor(metadata$subject_id))
    rownames(metadata) <- metadata$profile_id
    counts_cell_type <- as.matrix(counts[, selected, drop = FALSE])
    n_profiles <- nrow(metadata)
    n_subjects <- length(unique(metadata$subject_id))
    age_span <- diff(range(metadata$age))

    support_reason <- ""
    if (n_profiles < min_sample_units) support_reason <- "insufficient_sample_units"
    if (n_subjects < min_subjects) support_reason <- "insufficient_subjects"
    if (!is.finite(age_span) || age_span < min_age_span) support_reason <- "insufficient_age_span"
    design <- design_diagnostics(metadata, primary_formula)
    if (!design$estimable) support_reason <- design$reason
    if (isTRUE(design$estimable) && design$residual_df < min_residual_df) {
        support_reason <- "insufficient_fixed_effect_residual_df"
    }
    if (nzchar(support_reason)) {
        if (diagnostics_only) {
            stop("Diagnostic plot is not estimable for ", cell_type, ": ", support_reason)
        }
        message(sprintf("[%s] primary excluded: %s", cell_type, support_reason))
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
            cell_type = cell_type, analysis_tier = metadata$analysis_tier[[1L]],
            model_scope = "not_tested", n_sample_units = n_profiles,
            n_subjects = n_subjects, age_span_years = age_span,
            design_columns = ifelse(is.null(design$columns), NA, design$columns),
            design_rank = ifelse(is.null(design$rank), NA, design$rank),
            residual_df = ifelse(is.null(design$residual_df), NA, design$residual_df),
            n_genes_total = nrow(counts_cell_type), n_genes_tested = 0L,
            status = "excluded", reason = support_reason, stringsAsFactors = FALSE
        )
        next
    }

    if (diagnostics_only) {
        diagnostic_result <- tryCatch(
            run_voom_diagnostic(
                counts_cell_type,
                metadata,
                design$formula,
                design$fixed_design,
                cell_type
            ),
            error = function(error) error
        )
        if (inherits(diagnostic_result, "error")) {
            stop(
                "Diagnostic plot failed for ",
                cell_type,
                ": ",
                conditionMessage(diagnostic_result)
            )
        }
        diagnostic_pages <- diagnostic_pages + 1L
        message(sprintf(
            "[%s] diagnostic plot completed for %d genes",
            cell_type,
            diagnostic_result
        ))
        next
    }

    model <- tryCatch(
        run_model(
            counts_cell_type,
            metadata,
            design$formula,
            design$fixed_design,
            cell_type
        ),
        error = function(error) error
    )
    if (inherits(model, "error")) {
        message(sprintf("[%s] primary model failed: %s", cell_type, conditionMessage(model)))
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
            cell_type = cell_type, analysis_tier = metadata$analysis_tier[[1L]],
            model_scope = "adjusted_repeated_measures", n_sample_units = n_profiles,
            n_subjects = n_subjects, age_span_years = age_span,
            design_columns = design$columns, design_rank = design$rank,
            residual_df = design$residual_df, n_genes_total = nrow(counts_cell_type),
            n_genes_tested = 0L, status = "model_failed",
            reason = conditionMessage(model), stringsAsFactors = FALSE
        )
        next
    }

    table <- model$table
    result <- data.frame(
        gene_id = rownames(table), cell_type = cell_type,
        analysis_tier = metadata$analysis_tier[[1L]],
        model_scope = "adjusted_repeated_measures",
        model_formula = primary_formula,
        log2_fc_per_10_years = table$logFC,
        average_log_cpm = table$AveExpr,
        moderated_t_statistic = table$t,
        z_standardized = if ("z.std" %in% colnames(table)) table$z.std else NA_real_,
        p_value = table$P.Value,
        fdr_within_celltype = table$adj.P.Val,
        n_sample_units = n_profiles, n_subjects = n_subjects,
        age_span_years = age_span, stringsAsFactors = FALSE
    )
    results[[length(results) + 1L]] <- result
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        cell_type = cell_type, analysis_tier = metadata$analysis_tier[[1L]],
        model_scope = "adjusted_repeated_measures", n_sample_units = n_profiles,
        n_subjects = n_subjects, age_span_years = age_span,
        design_columns = design$columns, design_rank = design$rank,
        residual_df = design$residual_df, n_genes_total = nrow(counts_cell_type),
        n_genes_tested = model$n_tested, status = "completed", reason = "",
        stringsAsFactors = FALSE
    )

    sensitivity_order <- order(
        as.character(metadata$subject_id),
        metadata$age,
        as.character(metadata$sample_unit_id)
    )
    first_subject_row <- !duplicated(
        as.character(metadata$subject_id[sensitivity_order])
    )
    sensitivity_indices <- sensitivity_order[first_subject_row]
    sensitivity_metadata <- metadata[sensitivity_indices, , drop = FALSE]
    sensitivity_metadata$sex <- droplevels(sensitivity_metadata$sex)
    sensitivity_metadata$batch <- droplevels(sensitivity_metadata$batch)
    sensitivity_metadata$subject_id <- droplevels(sensitivity_metadata$subject_id)
    sensitivity_counts <- counts_cell_type[, sensitivity_indices, drop = FALSE]
    sensitivity_formula <- reformulas::nobars(design$formula)
    sensitivity_design <- design_diagnostics(
        sensitivity_metadata,
        sensitivity_formula,
        require_repeated = FALSE
    )
    sensitivity_age_span <- diff(range(sensitivity_metadata$age))
    sensitivity_reason <- ""
    if (nrow(sensitivity_metadata) < min_subjects) {
        sensitivity_reason <- "insufficient_subjects"
    }
    if (!is.finite(sensitivity_age_span) || sensitivity_age_span < min_age_span) {
        sensitivity_reason <- "insufficient_age_span"
    }
    if (!sensitivity_design$estimable) sensitivity_reason <- sensitivity_design$reason
    if (
        isTRUE(sensitivity_design$estimable) &&
        sensitivity_design$residual_df < min_residual_df
    ) {
        sensitivity_reason <- "insufficient_fixed_effect_residual_df"
    }
    if (nzchar(sensitivity_reason)) {
        message(sprintf("[%s] sensitivity excluded: %s", cell_type, sensitivity_reason))
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
            cell_type = cell_type, analysis_tier = metadata$analysis_tier[[1L]],
            model_scope = "one_sample_per_subject", n_sample_units = nrow(sensitivity_metadata),
            n_subjects = nrow(sensitivity_metadata), age_span_years = sensitivity_age_span,
            design_columns = ifelse(
                is.null(sensitivity_design$columns), NA, sensitivity_design$columns
            ),
            design_rank = ifelse(
                is.null(sensitivity_design$rank), NA, sensitivity_design$rank
            ),
            residual_df = ifelse(
                is.null(sensitivity_design$residual_df), NA, sensitivity_design$residual_df
            ),
            n_genes_total = nrow(sensitivity_counts), n_genes_tested = 0L,
            status = "excluded", reason = sensitivity_reason, stringsAsFactors = FALSE
        )
    } else {
        sensitivity_model <- tryCatch(
            run_fixed_model(
                sensitivity_counts,
                sensitivity_metadata,
                sensitivity_formula,
                sensitivity_design$fixed_design
            ),
            error = function(error) error
        )
        if (inherits(sensitivity_model, "error")) {
            message(sprintf(
                "[%s] sensitivity model failed: %s",
                cell_type,
                conditionMessage(sensitivity_model)
            ))
            diagnostics[[length(diagnostics) + 1L]] <- data.frame(
                cell_type = cell_type, analysis_tier = metadata$analysis_tier[[1L]],
                model_scope = "one_sample_per_subject",
                n_sample_units = nrow(sensitivity_metadata),
                n_subjects = nrow(sensitivity_metadata),
                age_span_years = sensitivity_age_span,
                design_columns = sensitivity_design$columns,
                design_rank = sensitivity_design$rank,
                residual_df = sensitivity_design$residual_df,
                n_genes_total = nrow(sensitivity_counts), n_genes_tested = 0L,
                status = "model_failed", reason = conditionMessage(sensitivity_model),
                stringsAsFactors = FALSE
            )
        } else {
            sensitivity_table <- sensitivity_model$table
            sensitivity_result <- data.frame(
                gene_id = rownames(sensitivity_table), cell_type = cell_type,
                analysis_tier = metadata$analysis_tier[[1L]],
                model_scope = "one_sample_per_subject",
                model_formula = paste(deparse(sensitivity_formula), collapse = " "),
                selection_rule = sensitivity_sample_rule,
                log2_fc_per_10_years = sensitivity_table$logFC,
                average_log_cpm = sensitivity_table$AveExpr,
                moderated_t_statistic = sensitivity_table$t,
                p_value = sensitivity_table$P.Value,
                fdr_within_celltype = sensitivity_table$adj.P.Val,
                n_sample_units = nrow(sensitivity_metadata),
                n_subjects = nrow(sensitivity_metadata),
                age_span_years = sensitivity_age_span,
                stringsAsFactors = FALSE
            )
            sensitivity_results[[length(sensitivity_results) + 1L]] <- sensitivity_result
            diagnostics[[length(diagnostics) + 1L]] <- data.frame(
                cell_type = cell_type, analysis_tier = metadata$analysis_tier[[1L]],
                model_scope = "one_sample_per_subject",
                n_sample_units = nrow(sensitivity_metadata),
                n_subjects = nrow(sensitivity_metadata),
                age_span_years = sensitivity_age_span,
                design_columns = sensitivity_design$columns,
                design_rank = sensitivity_design$rank,
                residual_df = sensitivity_design$residual_df,
                n_genes_total = nrow(sensitivity_counts),
                n_genes_tested = sensitivity_model$n_tested,
                status = "completed", reason = "", stringsAsFactors = FALSE
            )
        }
    }
}
grDevices::dev.off()
dir.create(dirname(plot_out), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(plot_temp, plot_out, overwrite = TRUE, copy.mode = FALSE)) {
    stop("Failed to copy the diagnostic PDF from container-local temporary storage.")
}
if (!identical(unname(tools::md5sum(plot_temp)), unname(tools::md5sum(plot_out)))) {
    stop("Diagnostic PDF checksum changed while copying to the configured output path.")
}
unlink(plot_temp)

if (diagnostics_only) {
    if (diagnostic_pages != length(cell_types)) {
        stop("Diagnostic-only run did not produce one plot for every eligible cell type.")
    }
    message(sprintf(
        "[pseudobulk_dream] regenerated %d diagnostic pages",
        diagnostic_pages
    ))
    quit(save = "no", status = 0L)
}

diagnostics_table <- do.call(rbind, diagnostics)
write_csv(diagnostics_table, diagnostics_out)
if (!length(results)) stop("No cell type produced a repeated-measures result table.")
combined <- do.call(rbind, results)
combined$fdr_global <- stats::p.adjust(combined$p_value, method = "BH")
combined <- combined[order(combined$analysis_tier, combined$cell_type, combined$p_value, combined$gene_id), ]
write_csv(combined, combined_out)

if (!length(sensitivity_results)) {
    stop("No cell type produced a one-sample-per-subject sensitivity result table.")
}
sensitivity_combined <- do.call(rbind, sensitivity_results)
sensitivity_combined$fdr_global <- stats::p.adjust(
    sensitivity_combined$p_value, method = "BH"
)
sensitivity_combined <- sensitivity_combined[
    order(
        sensitivity_combined$analysis_tier,
        sensitivity_combined$cell_type,
        sensitivity_combined$p_value,
        sensitivity_combined$gene_id
    ),
]
write_csv(sensitivity_combined, sensitivity_out)

manifest <- list()
for (index in seq_along(cell_types)) {
    cell_type <- cell_types[[index]]
    cell_result <- combined[combined$cell_type == cell_type, , drop = FALSE]
    if (!nrow(cell_result)) next
    output_filename <- sprintf("%02d_%s.csv.gz", index, safe_label(cell_type))
    write_csv(cell_result, file.path(celltype_dir, output_filename))
    manifest[[length(manifest) + 1L]] <- data.frame(
        cell_type = cell_type, analysis_tier = cell_result$analysis_tier[[1L]],
        model_scope = cell_result$model_scope[[1L]], n_genes_tested = nrow(cell_result),
        result_path = file.path("by_cell_type", output_filename), stringsAsFactors = FALSE
    )
}
write_csv(do.call(rbind, manifest), manifest_out)

dir.create(dirname(session_info_out), recursive = TRUE, showWarnings = FALSE)
session_connection <- file(session_info_out, open = "wt")
sink(session_connection)
print(sessionInfo())
sink()
close(session_connection)
message(sprintf(
    "[pseudobulk_dream] completed %d cell types and %d gene-level tests",
    length(unique(combined$cell_type)), nrow(combined)
))
