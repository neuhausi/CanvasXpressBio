#' CanvasXpress visualizations for Bioconductor data structures
#'
#' Bridges the CanvasXpress interactive visualization library to core
#' Bioconductor classes. [cxplot()] renders a
#' [SummarizedExperiment::SummarizedExperiment] (and any class that extends it,
#' such as `SingleCellExperiment` or `DESeqDataSet`) as an interactive
#' CanvasXpress heatmap carrying its row and column annotations; [cxvolcano()]
#' draws a volcano plot from a differential-expression table; [cxsurvival()]
#' draws a Kaplan-Meier curve from `colData`; and [cxboxplot()] draws a paired
#' boxplot of one feature with per-subject connections.
#'
#' @name CanvasXpressBio-package
#' @keywords internal
#' @importFrom methods setGeneric setMethod is
#' @importFrom stats var rnorm runif rexp rbinom t.test setNames
#' @importFrom canvasXpress canvasXpress
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
"_PACKAGE"


#' Visualize a Bioconductor object with CanvasXpress
#'
#' Renders a Bioconductor data object as an interactive CanvasXpress chart. The
#' assay matrix becomes the chart data, with `colData` supplied as sample
#' (column) annotations and `rowData` as variable (row) annotations, so the
#' annotations are available for interactive grouping, coloring and sorting.
#'
#' @param object A Bioconductor data object. Currently a
#'   [SummarizedExperiment::SummarizedExperiment] or any class extending it
#'   (e.g. `SingleCellExperiment`, `DESeqDataSet`).
#' @param ... Additional CanvasXpress configuration passed through to
#'   [canvasXpress::canvasXpress()] (for example `title` or `colorScheme`).
#'
#' @return A `canvasXpress` htmlwidget.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays  = list(counts = matrix(rnorm(20), nrow = 5,
#'                    dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))),
#'     colData = data.frame(group = rep(c("A", "B"), 2),
#'                          row.names = paste0("s", 1:4)))
#' cxplot(se, title = "Example heatmap")
#'
#' @export
setGeneric("cxplot", function(object, ...) standardGeneric("cxplot"))


# Build the overlay / clustering config for cxplot. `overlays` and
# `var_overlays` may be TRUE (all annotation columns) or a character vector of
# column names; `cluster` toggles row and column dendrograms.
.cxOverlayArgs <- function(overlays, var_overlays, cluster,
    smp_annot, var_annot, use_var_annot) {
    pick <- function(request, annot) {
        if (isTRUE(request)) {
            return(colnames(annot))
        }
        if (is.character(request)) {
            return(intersect(request, colnames(annot)))
        }
        NULL
    }
    smp <- pick(overlays, smp_annot)
    var <- if (use_var_annot) pick(var_overlays, var_annot) else NULL

    extra <- list()
    if (length(smp)) extra$smpOverlays <- as.list(smp)
    if (length(var)) extra$varOverlays <- as.list(var)
    if (isTRUE(cluster)) {
        extra$samplesClustered <- TRUE
        extra$variablesClustered <- TRUE
    }
    extra
}


#' @rdname cxplot
#'
#' @param assay Assay to plot, given as a name or an index. Defaults to the
#'   first assay.
#' @param graphType CanvasXpress graph type. Defaults to `"Heatmap"`.
#' @param n Optional integer. When supplied, only the `n` most variable
#'   features (rows) are kept before plotting — a common convenience for large
#'   genomics matrices. Defaults to `NULL`, which keeps all features.
#' @param overlays Sample (`colData`) annotations to draw as overlay tracks
#'   beside the heatmap. `TRUE` shows every `colData` column; a character
#'   vector names specific ones; `NULL` (the default) shows none.
#' @param var_overlays Feature (`rowData`) annotations to draw as overlay
#'   tracks, accepting the same values as `overlays`.
#' @param cluster Logical; when `TRUE` both samples and features are
#'   hierarchically clustered (adding dendrograms). Defaults to `FALSE`.
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @export
setMethod(
    "cxplot", "SummarizedExperiment",
    function(object, assay = 1L, graphType = "Heatmap", n = NULL,
        overlays = NULL, var_overlays = NULL, cluster = FALSE, ...) {
        mat <- as.matrix(SummarizedExperiment::assay(object, assay))
        if (is.null(rownames(mat))) {
            rownames(mat) <- paste0("V", seq_len(nrow(mat)))
        }
        if (is.null(colnames(mat))) {
            colnames(mat) <- paste0("S", seq_len(ncol(mat)))
        }

        var_annot <- as.data.frame(SummarizedExperiment::rowData(object))
        smp_annot <- as.data.frame(SummarizedExperiment::colData(object))
        rownames(smp_annot) <- colnames(mat)
        if (nrow(var_annot) == nrow(mat)) {
            rownames(var_annot) <- rownames(mat)
        }

        # Keep the n most variable features when requested.
        if (!is.null(n) && is.finite(n) && n < nrow(mat)) {
            vars_by_row <- apply(mat, 1L, stats::var)
            keep <- sort(order(vars_by_row, decreasing = TRUE)[seq_len(n)])
            if (nrow(var_annot) == length(vars_by_row)) {
                var_annot <- var_annot[keep, , drop = FALSE]
            }
            mat <- mat[keep, , drop = FALSE]
        }

        use_var_annot <- nrow(var_annot) == nrow(mat) && ncol(var_annot)
        extra <- .cxOverlayArgs(overlays, var_overlays, cluster,
            smp_annot, var_annot, use_var_annot)

        do.call(canvasXpress::canvasXpress, c(
            list(
                data = mat,
                smpAnnot = if (ncol(smp_annot)) smp_annot else NULL,
                varAnnot = if (use_var_annot) var_annot else NULL,
                graphType = graphType),
            extra,
            list(...)))
    })


# Build the volcano's visual elements: the Regulation annotation, the
# red / blue / grey colour key, and the threshold-gate lines (vertical at the
# fold-change cut-offs, horizontal at the p-value cut-off).
.cxVolcanoVisuals <- function(lfc, nlp, label_fc, label_p, ids) {
    regulation <- ifelse(
        nlp > label_p & lfc > label_fc, "Up-regulated",
        ifelse(nlp > label_p & lfc < -label_fc,
            "Down-regulated", "Not significant"))
    annot <- data.frame(
        Regulation = regulation,
        row.names = ids,
        check.names = FALSE,
        stringsAsFactors = FALSE)
    key <- list(Regulation = list(
        "Up-regulated" = "rgba(197,27,38,0.75)",
        "Down-regulated" = "rgba(33,102,172,0.75)",
        "Not significant" = "rgba(150,150,150,0.35)"))
    gate <- "rgba(120,120,120,0.8)"
    decorations <- list(line = list(
        list(color = gate, width = 1, x = label_fc),
        list(color = gate, width = 1, x = -label_fc),
        list(color = gate, width = 1, y = label_p)))
    list(annot = annot, key = key, decorations = decorations)
}


# Build the labelBy / labelSelect arguments for a volcano. Returns the config
# list plus the optional set of genes to flag for an explicit-label selection.
.cxVolcanoLabels <- function(label, label_genes, label_fc, label_p) {
    if (!isTRUE(label)) {
        return(list(args = list(), highlight = NULL))
    }
    if (!is.null(label_genes)) {
        select <- list(list("Highlight", "==", "yes"))
        return(list(
            args = list(
                labelBy = "vars",
                labelSelect = select,
                optimizeTextPosition = TRUE),
            highlight = label_genes))
    }
    select <- list(
        "AND",
        list("y", ">", label_p),
        list("OR",
            list("x", "<", -label_fc),
            list("x", ">", label_fc)))
    list(
        args = list(
            labelBy = "vars",
            labelSelect = select,
            optimizeTextPosition = TRUE),
        highlight = NULL)
}


#' Volcano plot of differential-expression results with CanvasXpress
#'
#' Draws an interactive volcano plot (log fold change against
#' \eqn{-\log_{10}} p-value) from a table of differential-expression results.
#' The input is coerced with [as.data.frame()], so a `DESeqResults` object, a
#' `limma::topTable()` data frame, or any comparable table works.
#'
#' @param x A data-frame-like table of results.
#' @param logfc,pval Column names holding the log fold change and the
#'   (adjusted) p-value. The defaults suit `DESeq2` output.
#' @param label Logical; when `TRUE` (the default) genes are labelled with
#'   their row names. Set to `FALSE` to draw an unlabelled volcano.
#' @param label_genes Optional character vector of gene (row) names to label.
#'   When supplied, only these genes are labelled — useful for calling out a
#'   named set on a dense plot. When `NULL` (the default) every significant
#'   gene is labelled instead.
#' @param label_fc,label_p Thresholds that decide which genes are labelled: a
#'   gene is labelled when its absolute log fold change exceeds `label_fc` and
#'   its \eqn{-\log_{10}} p-value exceeds `label_p`. Defaults of `1` and `2`
#'   label points outside a two-fold change and below p \eqn{\approx} 0.01.
#' @param ... Additional CanvasXpress configuration passed to
#'   [canvasXpress::canvasXpress()].
#'
#' @return A `canvasXpress` htmlwidget.
#'
#' @examples
#' res <- data.frame(
#'     log2FoldChange = rnorm(100),
#'     padj           = runif(100),
#'     row.names      = paste0("gene", 1:100))
#' cxvolcano(res, title = "Volcano")
#'
#' @export
cxvolcano <- function(x, logfc = "log2FoldChange", pval = "padj",
    label = TRUE, label_genes = NULL, label_fc = 1, label_p = 2, ...) {
    df <- as.data.frame(x)
    if (!all(c(logfc, pval) %in% colnames(df))) {
        stop("Columns '", logfc, "' and '", pval,
            "' must both be present in 'x'.", call. = FALSE)
    }
    ok <- is.finite(df[[logfc]]) & is.finite(df[[pval]]) & df[[pval]] > 0
    df <- df[ok, , drop = FALSE]
    if (is.null(rownames(df))) {
        rownames(df) <- paste0("F", seq_len(nrow(df)))
    }

    # Features (genes) are the rows (`vars`); the fold change and p-value are
    # the columns (`smps`) that drive the axes. Do NOT transpose — that would
    # collapse the whole table to two plotted points.
    lfc <- df[[logfc]]
    nlp <- -log10(df[[pval]])
    ids <- rownames(df)
    mat <- as.matrix(data.frame(
        "log2FC" = lfc,
        "negLogP" = nlp,
        row.names = ids,
        check.names = FALSE))

    vis <- .cxVolcanoVisuals(lfc, nlp, label_fc, label_p, ids)
    lab <- .cxVolcanoLabels(label, label_genes, label_fc, label_p)
    annot <- vis$annot
    if (!is.null(lab$highlight)) {
        annot$Highlight <- ifelse(ids %in% lab$highlight, "yes", "no")
    }

    do.call(canvasXpress::canvasXpress, c(
        list(
            data = mat,
            varAnnot = annot,
            graphType = "Scatter2D",
            xAxis = list("log2FC"),
            yAxis = list("negLogP"),
            xAxisTitle = "log2 fold change",
            yAxisTitle = "-log10 p-value",
            colorBy = "Regulation",
            colorKey = vis$key,
            showDecorations = TRUE,
            decorations = vis$decorations),
        lab$args,
        list(...)))
}


#' Kaplan-Meier survival plot with CanvasXpress
#'
#' Draws an interactive Kaplan-Meier survival curve from per-subject survival
#' data, optionally stratified into groups. Survival columns are read from the
#' `colData` of a [SummarizedExperiment::SummarizedExperiment] (so the same
#' object that feeds [cxplot()] can produce the survival plot) or from a plain
#' data frame.
#'
#' @param object A `SummarizedExperiment` (its `colData` is used) or a
#'   data-frame-like table with one row per subject.
#' @param time,status Column names holding the follow-up time and the event
#'   status (`1`/`TRUE` = event, `0`/`FALSE` = censored).
#' @param group Optional column name used to stratify and colour the curves
#'   (for example a treatment arm or cohort). `NULL` draws a single curve.
#' @param tooltip Optional character vector of extra `colData` column names to
#'   carry with each subject, so hovering a step or a censored observation
#'   shows that subject's metadata (for example age, sex or disease stage).
#'   `NULL` (the default) carries only the grouping column.
#' @param ... Additional CanvasXpress configuration passed to
#'   [canvasXpress::canvasXpress()].
#'
#' @return A `canvasXpress` htmlwidget.
#'
#' @examples
#' cd <- data.frame(
#'     time   = c(5, 12, 18, 24, 7, 15, 22, 30),
#'     status = c(1, 1, 0, 1, 1, 0, 1, 0),
#'     arm    = rep(c("A", "B"), each = 4))
#' cxsurvival(cd, time = "time", status = "status", group = "arm")
#'
#' @importFrom SummarizedExperiment colData
#' @export
cxsurvival <- function(object, time = "time", status = "status",
    group = NULL, tooltip = NULL, ...) {
    cd <- if (methods::is(object, "SummarizedExperiment")) {
        as.data.frame(SummarizedExperiment::colData(object))
    } else {
        as.data.frame(object)
    }
    needed <- c(time, status, group, tooltip)
    if (!all(needed %in% colnames(cd))) {
        missing <- setdiff(needed, colnames(cd))
        stop("Columns ", paste(sQuote(missing), collapse = ", "),
            " are missing from the survival data.", call. = FALSE)
    }

    # CanvasXpress reads a Kaplan-Meier curve from a plot whose variables
    # (rows) are the subjects and whose `time` / `status` columns drive the
    # axes. Each subject is one row.
    subjects <- rownames(cd)
    if (is.null(subjects)) subjects <- paste0("S", seq_len(nrow(cd)))
    mat <- as.matrix(data.frame(
        time = as.numeric(cd[[time]]),
        status = as.numeric(cd[[status]]),
        row.names = subjects,
        check.names = FALSE))

    # Carry the grouping column and any requested tooltip columns as per-subject
    # annotations: CanvasXpress shows them in the hover tooltip for each step /
    # censored point, and the grouping column also colours the curves.
    annot_cols <- unique(c(group, tooltip))
    group_args <- list()
    if (length(annot_cols)) {
        var_annot <- data.frame(
            row.names = subjects,
            check.names = FALSE,
            stringsAsFactors = FALSE)
        for (col in annot_cols) var_annot[[col]] <- as.character(cd[[col]])
        group_args <- list(varAnnot = var_annot)
        if (!is.null(group)) group_args$colorBy <- group
    }

    do.call(canvasXpress::canvasXpress, c(
        list(
            data = mat,
            graphType = "KaplanMeier",
            xAxis = list("time"),
            yAxis = list("status"),
            xAxisTitle = "Time",
            yAxisTitle = "Survival probability",
            kmRiskTable = TRUE,
            showKMConfidenceIntervals = TRUE,
            showKMMedianSurvivalTime = TRUE),
        group_args,
        list(...)))
}


#' Paired boxplot with per-subject connections using CanvasXpress
#'
#' Draws a grouped boxplot of one measurement across conditions, overlaying the
#' individual points and drawing a line between the repeated measurements of
#' each subject. The measurement is one feature (row) of a
#' [SummarizedExperiment::SummarizedExperiment] assay, and the grouping,
#' faceting, connecting and colouring are taken from its `colData`.
#'
#' @param object A `SummarizedExperiment`.
#' @param feature The assay feature (row name) whose values are plotted.
#' @param group `colData` column that defines the box groups along the x-axis
#'   (for example a time point).
#' @param segregate Optional `colData` column that splits the plot into panels
#'   (for example a cohort).
#' @param connect Optional `colData` column identifying the subject, so repeated
#'   measurements are joined by a line.
#' @param color Optional `colData` column used to colour the points and
#'   connecting lines.
#' @param points How the individual observations are drawn over the boxes:
#'   `"all"` (the default) places every point, `"sina"` packs them into a
#'   beeswarm that stays legible with many points, `"jitter"` spreads them
#'   randomly, and `"none"` draws boxes only. The beeswarm and jitter layouts
#'   scale to hundreds or thousands of points without overplotting.
#' @param assay Assay to read, given as a name or index. Defaults to the first.
#' @param ... Additional CanvasXpress configuration passed to
#'   [canvasXpress::canvasXpress()].
#'
#' @return A `canvasXpress` htmlwidget.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays  = list(expr = matrix(rnorm(8), nrow = 1,
#'                    dimnames = list("IL6", paste0("s", 1:8)))),
#'     colData = data.frame(
#'         Subject   = rep(paste0("P", 1:4), each = 2),
#'         Condition = rep(c("Baseline", "Week 12"), 4),
#'         row.names = paste0("s", 1:8)))
#' cxboxplot(se, feature = "IL6", group = "Condition", connect = "Subject")
#'
#' @importFrom SummarizedExperiment assay colData
#' @export
cxboxplot <- function(object, feature, group, segregate = NULL,
    connect = NULL, color = NULL, points = c("all", "sina", "jitter", "none"),
    assay = 1L, ...) {
    points <- match.arg(points)
    if (!methods::is(object, "SummarizedExperiment")) {
        stop("'object' must be a SummarizedExperiment.", call. = FALSE)
    }
    mat_all <- as.matrix(SummarizedExperiment::assay(object, assay))
    if (!feature %in% rownames(mat_all)) {
        stop("Feature ", sQuote(feature), " is not a row of the assay.",
            call. = FALSE)
    }
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    needed <- c(group, segregate, connect, color)
    if (!all(needed %in% colnames(cd))) {
        missing <- setdiff(needed, colnames(cd))
        stop("Columns ", paste(sQuote(missing), collapse = ", "),
            " are missing from colData.", call. = FALSE)
    }

    # One measurement row (the chosen feature) across all samples; the sample
    # annotations carry the grouping / faceting / connection / colour factors.
    samples <- colnames(mat_all)
    if (is.null(samples)) samples <- paste0("S", seq_len(ncol(mat_all)))
    mat <- matrix(mat_all[feature, ], nrow = 1,
        dimnames = list(feature, samples))

    smp_annot <- data.frame(
        row.names = samples,
        check.names = FALSE,
        stringsAsFactors = FALSE)
    for (col in needed) smp_annot[[col]] <- as.character(cd[[col]])

    box_args <- list(
        graphType = "Boxplot",
        graphOrientation = "vertical",
        groupingFactors = list(group),
        xAxisTitle = feature)
    # Point overlay style. "sina" packs the observations into a beeswarm that
    # stays readable at scale; "jitter" spreads them randomly; "all" places
    # every point; "none" leaves the boxes bare.
    if (points == "sina") {
        box_args$showBoxplotOriginalData <- TRUE
        box_args$sina <- TRUE
    } else if (points == "jitter") {
        box_args$showBoxplotOriginalData <- TRUE
        box_args$jitter <- TRUE
    } else if (points == "all") {
        box_args$showBoxplotOriginalData <- TRUE
    }
    if (!is.null(segregate)) box_args$segregateSamplesBy <- list(segregate)
    if (!is.null(connect)) {
        box_args$connectBy <- connect
        box_args$connectByPointColor <- TRUE
        box_args$connectByWidth <- 1
    }
    if (!is.null(color)) box_args$colorBy <- color

    do.call(canvasXpress::canvasXpress, c(
        list(data = mat, smpAnnot = smp_annot),
        box_args,
        list(...)))
}


#' Build the demonstration SummarizedExperiment used in the vignette
#'
#' Assembles a small, self-contained `SummarizedExperiment` so the vignette can
#' show every figure from one object with a single function call instead of an
#' inline data-construction block. The study follows 24 patients in two cohorts,
#' each measured at baseline and at week 12 over a panel of twelve marker genes.
#' `colData` records the cohort, time point, subject, a follow-up time with an
#' event indicator, and the direction of each subject's biomarker change;
#' `rowData` records whether each gene is an up- or down-regulated marker.
#'
#' The derived analyses that the other figures need are attached to the object's
#' [S4Vectors::metadata()] so the vignette needs no further wrangling:
#' `deResults` is a differential-expression table (the twelve markers among a
#' synthetic genome-wide background) for [cxvolcano()], and `albPaired` holds
#' the per-cohort paired t-test p-values for the `ALB` marker shown by
#' [cxboxplot()]. The random draws are seeded reproducibly with
#' [withr::with_seed()], so the object is identical on every call.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] with `deResults` and
#'   `albPaired` in its `metadata()`.
#'
#' @examples
#' se <- makeDemoSE()
#' se
#' S4Vectors::metadata(se)$albPaired
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay colData rowData
#' @export
makeDemoSE <- function() {
    d <- withr::with_seed(1, .demoData())

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(expr = d$expr),
        colData = S4Vectors::DataFrame(d$coldata),
        rowData = S4Vectors::DataFrame(
            Direction = ifelse(rownames(d$expr) %in% d$up_genes, "Up", "Down"),
            row.names = rownames(d$expr)))

    S4Vectors::metadata(se) <- list(
        deResults = d$de,
        albPaired = c("Cohort A" = d$pA, "Cohort B" = d$pB))
    se
}


# Draw the demonstration data (expression matrix, sample table, DE table and
# the per-cohort paired p-values). Kept separate so [makeDemoSE()] can run it
# under a fixed seed via withr::with_seed() rather than calling set.seed() in
# package code.
.demoData <- function() {
    up_genes   <- c("IL6", "CXCL8", "MMP9", "TNF", "CCL2", "VEGFA")
    down_genes <- c("ALB", "CDH1", "FOXO1", "PPARG", "ADIPOQ", "LEP")
    genes      <- c(up_genes, down_genes)

    patients   <- c(sprintf("PA%02d", seq_len(12)),
        sprintf("PB%02d", seq_len(12)))
    cohort     <- rep(c("Cohort A", "Cohort B"), each = 12)
    names(cohort) <- patients

    # Sample table: each patient measured at Baseline and Week 12.
    coldata <- data.frame(
        Patient   = rep(patients, each = 2),
        Condition = rep(c("Baseline", "Week 12"), times = length(patients)),
        stringsAsFactors = FALSE)
    coldata$Cohort <- cohort[coldata$Patient]
    rownames(coldata) <- paste(coldata$Patient,
        ifelse(coldata$Condition == "Baseline", "BL", "W12"),
        sep = "_")

    # Expression: baseline around 8; at week 12 the up genes rise and the down
    # genes fall. Cohort A responds consistently; cohort B is heterogeneous, so
    # the biomarker moves both ways there.
    baseline <- matrix(rnorm(length(genes) * length(patients), 8, 0.8),
        nrow = length(genes), dimnames = list(genes, patients))
    effect   <- ifelse(genes %in% up_genes, 1.6, -1.6)
    strength <- setNames(rep(1.0, length(patients)), patients)
    strength[cohort == "Cohort B"] <-
        runif(sum(cohort == "Cohort B"), -0.7, 0.8)

    expr <- matrix(0, nrow = length(genes), ncol = nrow(coldata),
        dimnames = list(genes, rownames(coldata)))
    for (j in seq_len(nrow(coldata))) {
        p <- coldata$Patient[j]
        delta <- if (coldata$Condition[j] == "Week 12")
            effect * strength[[p]] + rnorm(length(genes), 0, 0.4) else 0
        expr[, j] <- baseline[, p] + delta
    }

    # Per-patient survival (cohort A does better) and the direction of the ALB
    # change from baseline to week 12 (ALB is the biomarker in the boxplot).
    alb_bl  <- expr["ALB", coldata$Condition == "Baseline"]
    alb_w12 <- expr["ALB", coldata$Condition == "Week 12"]
    updown  <- ifelse(alb_w12 >= alb_bl, "up", "down")
    surv_time <- round(ifelse(cohort == "Cohort A",
        rexp(length(patients), 1 / 28),
        rexp(length(patients), 1 / 16)) + 2)
    surv_status <- rbinom(length(patients), 1,
        ifelse(cohort == "Cohort A", .5, .75))

    coldata$Time   <- surv_time[match(coldata$Patient, patients)]
    coldata$Status <- surv_status[match(coldata$Patient, patients)]
    coldata$updown <- updown[match(coldata$Patient, patients)]

    # Differential-expression table for the volcano: the observed marker effects
    # placed among a synthetic genome-wide background.
    marker_fc <- rowMeans(expr[, coldata$Condition == "Week 12"]) -
        rowMeans(expr[, coldata$Condition == "Baseline"])
    marker_p  <- 10^(-abs(marker_fc) * runif(length(genes), 2.5, 3.5))
    n_bg  <- 1988
    bg_fc <- rnorm(n_bg, 0, 0.45)
    bg_p  <- pmin(10^(-(abs(bg_fc) * runif(n_bg, 1.0, 3.0))), 1)
    de <- data.frame(
        log2FoldChange = c(marker_fc, bg_fc),
        padj           = c(marker_p, bg_p),
        row.names      = c(genes, paste0("bg", seq_len(n_bg))))

    # Per-cohort paired t-test p-values for ALB (shown on the boxplot panels).
    pA <- t.test(alb_w12[cohort == "Cohort A"], alb_bl[cohort == "Cohort A"],
        paired = TRUE)$p.value
    pB <- t.test(alb_w12[cohort == "Cohort B"], alb_bl[cohort == "Cohort B"],
        paired = TRUE)$p.value

    list(expr = expr, coldata = coldata, up_genes = up_genes,
        de = de, pA = pA, pB = pB)
}


#' Build a dense single-feature SummarizedExperiment for a beeswarm demo
#'
#' Returns a `SummarizedExperiment` with one feature measured across several
#' groups and many observations per group, so a boxplot with a sina / beeswarm
#' overlay ([cxboxplot()] with `points = "sina"`) can be shown to stay legible
#' at scale rather than collapsing into an overplotted mass.
#'
#' @param n_per_group Integer; observations in each group. Defaults to `500`.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] with a single
#'   `"Expression"` feature and a `Group` column in its `colData`.
#'
#' @examples
#' bse <- makeBeeswarmDemo(n_per_group = 100)
#' cxboxplot(bse, feature = "Expression", group = "Group", points = "sina")
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
makeBeeswarmDemo <- function(n_per_group = 500) {
    groups  <- c("Ctrl", "Low dose", "Mid dose", "High dose", "Combo")
    centers <- c(8, 8.6, 9.4, 10.2, 9.0)

    samples <- character(0)
    values  <- numeric(0)
    grp     <- character(0)
    for (gi in seq_along(groups)) {
        ids <- sprintf("%s_%04d", gsub(" ", "", groups[gi]),
            seq_len(n_per_group))
        samples <- c(samples, ids)
        grp     <- c(grp, rep(groups[gi], n_per_group))
        # Seed per group so the result is reproducible without set.seed() in
        # package code (withr restores the RNG state afterwards).
        values  <- c(values,
            withr::with_seed(gi, rnorm(n_per_group, centers[gi], 0.9)))
    }

    SummarizedExperiment::SummarizedExperiment(
        assays  = list(expr = matrix(values, nrow = 1,
            dimnames = list("Expression", samples))),
        colData = S4Vectors::DataFrame(Group = grp, row.names = samples))
}
