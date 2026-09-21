## ==========================================================================
## scriptname:  "neuroblastoma_efs.R"
## title:       "Elk-1 as a Potential Transcriptional Regulator of Cell
##              Surface Markers CD44, CD133, and CD114 in Neuroblastoma"
## author:      "CoRD-OMICS - CimenLab - Simay Akpinar"
## date:        "21.09.2026"
## license:     "MIT License (see LICENSE file for details)"
## description: "Expression and survival analysis of the ETS transcription
##              factors ELK1 and the PEA3 subfamily (ETV1, ETV4, ETV5) in
##              neuroblastoma, using public RNA-seq data from cell lines
##              (GSE89413) and patients (GSE62564). Part A: ELK1 expression
##              across cell lines, its comparison with the surface markers
##              CSF3R (CD114), CD44 and PROM1 (CD133), and a heatmap of the
##              most variable surface markers. Part B: stratification of
##              patients into ELK1-High and ELK1-Low groups (optimal cutpoint,
##              median or tertile), Kaplan-Meier analysis and time-dependent
##              ROC for event-free survival (EFS), correlation of ELK1 with
##              the focus markers and a wider marker panel, violin plots, and
##              MYCN expression across cell lines. Part C: the same pipeline
##              run independently for each PEA3 gene (ETV1, ETV4, ETV5),
##              without averaging."
##
## Run top to bottom. Set the paths in the config sections before running.
## ==========================================================================


## ---- packages ------------------------------------------------------------------
cran_pkgs <- c("ggplot2","dplyr","tidyr","tibble","stringr","readxl","ggrepel",
               "pheatmap","ggpubr","writexl","igraph","scales",
               "survival","survminer","timeROC","cmprsk")
to_inst <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_inst)) install.packages(to_inst, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("GEOquery","limma","org.Hs.eg.db","Biobase")
to_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_bioc)) BiocManager::install(to_bioc, ask = FALSE, update = FALSE)
suppressMessages(invisible(lapply(c("ggplot2","dplyr","tidyr","tibble","stringr","readxl",
                                    "ggrepel","pheatmap","ggpubr","writexl","survival","survminer",
                                    "GEOquery","limma","org.Hs.eg.db","Biobase","scales"),
                                  library, character.only = TRUE)))


## ################################################################################
## PART A - Cell lines (GSE89413): ELK1
## ################################################################################

## ---- config ------------------------------------------------------------------
## Set the paths below to your local directories.
DATA_DIR <- "path/to/your/data"
OUT_DIR  <- "path/to/your/output"

outdir  <- normalizePath(file.path(OUT_DIR, "GSE89413_outputs"), mustWork = FALSE)
FIG_DIR <- file.path(outdir, "Figures")
RES_DIR <- file.path(outdir, "Results")
for (d in c(outdir, FIG_DIR, RES_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
cat("Figure folder :", FIG_DIR, "\n"); cat("Result folder :", RES_DIR, "\n")

CELLLINE_FILE <- file.path(DATA_DIR, "GSE89413_2016-10-30-NBL-cell-line-STAR-fpkm.txt")

FOCUS_CELLS <- c("SKNBE2", "SHSY5Y", "KELLY")
FOCUS3      <- c("CSF3R", "CD44", "PROM1")

PAL_FOCUSCELL <- c(SKNBE2 = "#d6604d", SHSY5Y = "#4393c3", KELLY = "#1a9850")
hm_pal  <- colorRampPalette(c("navy","white","firebrick3"))(100)
GRP_PAL <- c(`Neuro+Stem`="#9970ab", Neuroblastoma="#74c476", Stemness="#fd8d3c", Other="grey80")
FOCUS_PAL<- c(focus = "#d6604d", other = "grey85")

theme_pub <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.35),
                   axis.text  = ggplot2::element_text(colour = "black", size = base_size - 1),
                   axis.title = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 1),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
                   legend.key = ggplot2::element_blank())
}

open_pdf <- function(path, w, h) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); grDevices::pdf(path, width = w, height = h) }
open_png <- function(path, w, h) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); grDevices::png(path, width = round(w*300), height = round(h*300), res = 300) }
FIGS <- list()
.draw_to <- function(opener, draw) {
  ok <- tryCatch({ opener(); TRUE }, error = function(e) { message("device could not be opened: ", conditionMessage(e)); FALSE })
  if (!ok) return(FALSE); on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  isTRUE(tryCatch({ draw(); TRUE }, error = function(e) { message("drawing error: ", conditionMessage(e)); FALSE }))
}
add_fig <- function(name, draw, w = 8, h = 6) {
  fp <- file.path(FIG_DIR, paste0(name, ".pdf")); pp <- file.path(FIG_DIR, paste0(name, ".png"))
  .draw_to(function() open_pdf(fp, w, h), draw); .draw_to(function() open_png(pp, w, h), draw)
  got <- c(if (file.exists(fp)) "pdf", if (file.exists(pp)) "png")
  cat("  ", name, " -> ", paste(got, collapse = "+"), if (length(got) < 2) "  (MISSING!)" else "", "\n", sep = "")
  FIGS[[name]] <<- list(w = w, h = h); invisible(name)
}
fig_pdf  <- function(plot, name, w = 8, h = 6) add_fig(name, function() print(plot), w, h)
note_fig <- function(name, w = 8, h = 6) { FIGS[[name]] <<- list(w = w, h = h); invisible(name) }
res_xlsx <- function(x, name) {
  f <- file.path(RES_DIR, paste0(name, ".xlsx")); writexl::write_xlsx(as.data.frame(x), f)
  if (file.exists(f)) cat("  -> table:", f, "\n"); invisible(f)
}



## ---- load_tpm ------------------------------------------------------------------
.norm_cell <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
read_tpm <- function(path) {
  stopifnot(file.exists(path))
  l1 <- readLines(path, n = 1, warn = FALSE)
  delim <- if (grepl("\t", l1)) "\t" else if (grepl(";", l1)) ";" else if (grepl(",", l1)) "," else ""
  df <- read.delim(path, sep = delim, header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  is_num_col <- vapply(df, function(c) mean(!is.na(suppressWarnings(as.numeric(c)))) > 0.9, logical(1))
  id_cols <- which(!is_num_col); num_cols <- which(is_num_col)
  if (length(num_cols) < 3) stop("Numeric column not found - sep might be wrong. l1: ", l1)
  sym_idx <- id_cols[ vapply(df[id_cols], function(c) any(grepl("[A-Za-z]", c)), logical(1)) ]
  sym_idx <- if (length(sym_idx)) sym_idx[1] else id_cols[1]
  sym <- toupper(trimws(as.character(df[[sym_idx]])))
  M <- as.matrix(df[, num_cols, drop = FALSE]); mode(M) <- "numeric"; rownames(M) <- sym
  if (any(duplicated(sym))) {
    ord <- order(rowMeans(M, na.rm = TRUE), decreasing = TRUE)
    M <- M[ord, , drop = FALSE]; M <- M[!duplicated(rownames(M)), , drop = FALSE]
  }
  M[rownames(M) != "" & !is.na(rownames(M)), , drop = FALSE]
}

tpm_raw <- read_tpm(CELLLINE_FILE)
rng <- range(tpm_raw, na.rm = TRUE)
cat(sprintf("Raw matrix: %d genes x %d cell lines | value range: %.2f .. %.2f\n",
            nrow(tpm_raw), ncol(tpm_raw), rng[1], rng[2]))
is_linear <- rng[2] > 30 || any(tpm_raw > 30, na.rm = TRUE)
expr <- if (is_linear) log2(tpm_raw + 1) else tpm_raw
cat(sprintf("Scale: %s -> expression matrix = %s\n", ifelse(is_linear,"LINEAR (FPKM)","already log-like"),
            ifelse(is_linear,"log2(FPKM+1)","as is")))
cat("\nCell lines:\n"); print(colnames(expr))
cat("\nAre focus genes present in the matrix?\n")
for (g in c("ELK1", FOCUS3)) cat(sprintf("  [%s] %s\n", g, ifelse(g %in% rownames(expr),"PRESENT","ABSENT")))
cell_lookup <- setNames(colnames(expr), .norm_cell(colnames(expr)))
focus_cells_real <- cell_lookup[.norm_cell(FOCUS_CELLS)]
cat("\nDid focus cells match?\n")
for (i in seq_along(FOCUS_CELLS))
  cat(sprintf("  %-8s -> %s\n", FOCUS_CELLS[i], ifelse(is.na(focus_cells_real[i]),"NOT FOUND",focus_cells_real[i])))
focus_cells_real <- as.character(na.omit(focus_cells_real))


## ---- gene_lists ------------------------------------------------------------------
read_markers <- function(file) {
  if (!file.exists(file)) return(character(0))
  m <- read.csv(file, stringsAsFactors = FALSE)
  gn <- grep("Gene.?Names", names(m), value = TRUE)[1]; if (is.na(gn)) gn <- names(m)[1]
  m[[gn]] %>% stringr::str_split("\\s+") %>% unlist() %>% toupper() %>% unique() %>% .[. != ""]
}
cd_markers   <- read_markers(file.path(DATA_DIR, "cd_antigen_markers.csv"))
surf_markers <- read_markers(file.path(DATA_DIR, "surface_marker_candidates.csv"))

## Neuroblastoma genes are read from neuroblastoma_markers.csv (Marker + Alias columns).
## Extra literature genes can be appended to neuro_genes below.
neuro_genes  <- character(0)
nmk_f <- file.path(DATA_DIR, "neuroblastoma_markers.csv")
if (file.exists(nmk_f)) {
  neuro_mk <- read.csv(nmk_f, stringsAsFactors = FALSE)
  neuro_genes <- c(neuro_mk$Marker,
                   stringr::str_split(neuro_mk$Alias, "[,;/ ]+", simplify = FALSE) %>% unlist()) %>%
    toupper() %>% stringr::str_trim() %>% unique() %>% .[. != "" & !is.na(.)]
}

## Stemness signatures (gene sets used for the Stemness / Neuro+Stem grouping).
## Edit freely: genes from the literature can be added, removed or replaced here,
## and whole signatures can be added. Use official HGNC gene symbols.
stem_sigs <- list(
  Pluripotency = toupper(c("POU5F1","SOX2","NANOG","KLF4","MYC","LIN28A","LIN28B",
                           "SALL4","ZFP42","FOXO1","DPPA4","UTF1","PRDM14","DNMT3L")),
  Neural_Stem  = toupper(c("NES","NOTCH1","NOTCH2","SOX2","PAX6","VIM","GFAP",
                           "PROM1","OLIG2","ASCL1","HES1","HES5","ID1","ID3")),
  Cancer_Stem  = toupper(c("CD44","PROM1","ALDH1A1","ALDH1A3","ITGA6","ABCG2",
                           "BMI1","EZH2","HMGA2","IGF2BP1","IGF2BP2","IGF2BP3","CD24","PCNA","MKI67")),
  MYCN_Pathway = toupper(c("MYCN","ODC1","LDHA","ENO1","CCND1","CDK4","MDM2",
                           "ID2","MAX","MXD1","E2F1")))
stem_genes  <- unique(unlist(stem_sigs))

group_of <- function(g) {
  inN <- g %in% neuro_genes; inS <- g %in% stem_genes
  ifelse(inN & inS, "Neuro+Stem", ifelse(inN,"Neuroblastoma", ifelse(inS,"Stemness","Other")))
}
col_ann   <- data.frame(Focus = ifelse(colnames(expr) %in% focus_cells_real,"focus","other"),
                        row.names = colnames(expr))
ann_colors<- list(Focus = FOCUS_PAL, Group = GRP_PAL)
cat("neuro:", length(neuro_genes), "| stem:", length(stem_genes),
    "| CD:", length(cd_markers), "| surface:", length(surf_markers), "\n")

## ---- elk1_across_cells ------------------------------------------------------------------
stopifnot("ELK1" %in% rownames(expr))
elk1_df <- data.frame(Cell = colnames(expr), ELK1 = as.numeric(expr["ELK1", ])) %>%
  dplyr::arrange(ELK1) %>%
  dplyr::mutate(Cell = factor(Cell, levels = Cell),
                Focus = ifelse(as.character(Cell) %in% focus_cells_real, "focus", "other"))
p_elk1 <- ggplot(elk1_df, aes(ELK1, Cell, fill = Focus)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = FOCUS_PAL, name = NULL) +
  labs(title = "ELK1 expression - neuroblastoma cell lines",
       subtitle = "log2(FPKM+1) | red = SKNBE2 / SHSY5Y / KELLY",
       x = "ELK1 (log2 FPKM)", y = NULL) + theme_pub(11)
print(p_elk1); fig_pdf(p_elk1, "02_ELK1_across_celllines", 7, max(6, ncol(expr)*0.22))

res_xlsx(elk1_df %>% dplyr::select(Cell, ELK1), "ELK1_expression_per_cellline")

## ---- focus_combined ------------------------------------------------------------------
combo_genes <- intersect(c("ELK1", FOCUS3), rownames(expr))
cb_long <- as.data.frame(expr[combo_genes, , drop = FALSE]) %>% tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Cell", values_to = "Expr") %>%
  dplyr::mutate(Gene = factor(Gene, levels = combo_genes),
                Focus = ifelse(Cell %in% focus_cells_real, "focus", "other"),
                Cell  = factor(Cell, levels = levels(elk1_df$Cell)))
p_cb <- ggplot(cb_long, aes(Expr, Cell, fill = Focus)) +
  geom_col(width = 0.78, colour = "black", linewidth = 0.12) +
  facet_wrap(~ Gene, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = FOCUS_PAL, name = NULL) +
  labs(title = "ELK1 + CSF3R / CD44 / PROM1 - cell lines (ordered by ELK1)",
       subtitle = "log2(FPKM+1) | red = 3 focus lines", x = "log2 FPKM", y = NULL) +
  theme_pub(10)
print(p_cb); fig_pdf(p_cb, "03_ELK1_focus_across_cells", 11, max(6, ncol(expr)*0.22))

## ---- marker_heatmaps -------------------------------------------------------------
## Compact heatmap: fixed small size, row names hidden (500 genes are unreadable anyway).
draw_marker_hm <- function(genes, title, fname, add_elk1 = TRUE, w = 7, h = 6, compact = TRUE) {
  g <- intersect(genes, rownames(expr))
  if (add_elk1) g <- unique(c("ELK1", g))
  
  s <- expr[g, , drop = FALSE]
  s <- s[apply(s, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]
  
  cat(title, "- genes:", nrow(s), "\n")
  if (nrow(s) < 3) {
    message(title, ": <3 genes - skipped")
    return(invisible())
  }
  
  if (compact) {
    # Compact mode: fixed small size, row names hidden (for large gene lists)
    safe_h <- h
    show_rn <- FALSE
  } else {
    calc_h <- max(6, nrow(s)*0.13 + 2)
    safe_h <- min(calc_h, 45)
    show_rn <- nrow(s) <= 60
  }
  
  pheatmap(s, annotation_col = col_ann, annotation_colors = list(Focus = FOCUS_PAL),
           scale = "row", color = hm_pal, border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE,
           clustering_method = "ward.D2", show_rownames = show_rn, fontsize_row = 6, fontsize_col = 8,
           main = title, filename = file.path(FIG_DIR, paste0(fname, ".png")),
           width = w, height = safe_h)
  
  note_fig(fname, w, safe_h)
  res_xlsx(data.frame(Gene = rownames(s), s, check.names = FALSE), fname)
}

vg_surf <- apply(expr[intersect(surf_markers, rownames(expr)), ], 1, var, na.rm = TRUE)
top_surf <- names(sort(vg_surf, decreasing = TRUE))[1:min(500, length(vg_surf))]

draw_marker_hm(top_surf, "Surface Markers Top 500 (z-score) - all lines", "05_Surface_heatmap",
               w = 7, h = 6, compact = TRUE)


## ################################################################################
## PART B - Patients (GSE62564): ELK1
## ################################################################################

## ---- config ------------------------------------------------------------------
## Set the paths below to your local directories.
DATA_DIR <- "path/to/your/data"
OUT_DIR  <- "path/to/your/output"

outdir  <- normalizePath(file.path(OUT_DIR, "GSE62564_outputs"), mustWork = FALSE)
FIG_DIR <- file.path(outdir, "Figures")
RES_DIR <- file.path(outdir, "Results")
CACHE   <- file.path(outdir, "geo_cache")
SUPP    <- file.path(outdir, "GSE62564_Supp")
for (d in c(outdir, FIG_DIR, RES_DIR, CACHE, SUPP)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
cat("Figure folder :", FIG_DIR, "\n"); cat("Result folder :", RES_DIR, "\n")

PRIMARY_COMPARISON <- "ELK1"
ELK1_SPLIT         <- "optimal"
FOCUS3             <- c("CSF3R", "CD44", "PROM1")
ELK1_REFSEQS       <- c("NM_005229","NM_001114123","NM_001283069")

PAL_ELK1 <- c("ELK1 High" = "#d6604d", "ELK1 Low" = "#4393c3")
PAL_COMP <- c("High Risk"="#d6604d","Low Risk"="#4393c3","MYCN Amp"="#b2182b",
              "MYCN Non-Amp"="#4dac26","Stage 4"="#b2182b","Stage 4S"="#4393c3")
hm_pal <- colorRampPalette(c("navy","white","firebrick3"))(100)

theme_pub <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.35),
                   axis.text  = ggplot2::element_text(colour = "black", size = base_size - 1),
                   axis.title = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 1),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
                   legend.key = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey75"),
                   strip.text = ggplot2::element_text(face = "bold", size = base_size - 1))
}

open_pdf <- function(path, w, h) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); grDevices::pdf(path, width = w, height = h) }
open_png <- function(path, w, h) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); grDevices::png(path, width = round(w*300), height = round(h*300), res = 300) }
FIGS <- list()
.draw_to <- function(opener, draw) {
  ok <- tryCatch({ opener(); TRUE }, error = function(e) { message("device could not be opened: ", conditionMessage(e)); FALSE })
  if (!ok) return(FALSE)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  isTRUE(tryCatch({ draw(); TRUE }, error = function(e) { message("drawing error: ", conditionMessage(e)); FALSE }))
}
add_fig <- function(name, draw, w = 8, h = 6) {
  fp <- file.path(FIG_DIR, paste0(name, ".pdf")); pp <- file.path(FIG_DIR, paste0(name, ".png"))
  .draw_to(function() open_pdf(fp, w, h), draw); .draw_to(function() open_png(pp, w, h), draw)
  got <- c(if (file.exists(fp)) "pdf", if (file.exists(pp)) "png")
  cat("  ", name, " -> ", paste(got, collapse = "+"), if (length(got) < 2) "  (MISSING!)" else "", "\n", sep = "")
  FIGS[[name]] <<- list(w = w, h = h); invisible(name)
}
fig_pdf  <- function(plot, name, w = 8, h = 6) add_fig(name, function() print(plot), w, h)
res_xlsx <- function(x, name) {
  f <- file.path(RES_DIR, paste0(name, ".xlsx")); writexl::write_xlsx(as.data.frame(x), f)
  if (file.exists(f)) cat("  -> table:", f, "\n"); invisible(f)
}

## ---- download_diag ------------------------------------------------------------------
cat("[1] Downloading GSE62564 clinical data (pData)...\n")
gse <- GEOquery::getGEO("GSE62564", GSEMatrix = TRUE, getGPL = FALSE, destdir = CACHE)
if (is.list(gse)) gse <- gse[[1]]
pdata <- Biobase::pData(gse)
cat(sprintf("Clinical: %d samples\n", nrow(pdata)))
cat("\npData columns (first 40):\n"); print(head(names(pdata), 40))

cat("\n[1] Downloading supplementary RNA-seq expression matrix (~15 MB)...\n")

supp_paths <- GEOquery::getGEOSuppFiles("GSE62564", baseDir = SUPP, makeDirectory = FALSE)
expr_file  <- rownames(supp_paths)[grepl("txt.gz|tsv.gz|csv.gz", rownames(supp_paths), ignore.case=TRUE)][1]
cat("Supplementary file:", basename(expr_file), "\n")

l1    <- readLines(expr_file, n = 1, warn = FALSE)
delim <- if (grepl("\t", l1)) "\t" else if (grepl(";", l1)) ";" else ","
raw   <- read.delim(expr_file, sep = delim, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
id_raw <- as.character(raw[[1]])
M      <- as.matrix(raw[, -1, drop = FALSE]); mode(M) <- "numeric"; rownames(M) <- id_raw

is_refseq <- mean(grepl("^N[MR]_", id_raw)) > 0.5
cat(sprintf("\nExpression matrix (raw): %d rows x %d samples | ID type: %s\n",
            nrow(M), ncol(M), ifelse(is_refseq, "RefSeq (NM_/NR_)", "symbol/other")))
cat("Sample column names (first 5):\n"); print(head(colnames(M), 5))

rng <- range(M, na.rm = TRUE)
is_linear <- rng[2] > 30 || any(M > 30, na.rm = TRUE)
cat(sprintf("Value range: %.2f..%.2f -> %s\n", rng[1], rng[2],
            ifelse(is_linear, "LINEAR -> log2(x+1)", "already log-like")))

## ---- gene_matrix ------------------------------------------------------------------
if (is_refseq) {
  acc <- sub("\\.\\d+$", "", rownames(M))
  map <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=unique(acc), columns="SYMBOL", keytype="REFSEQ"),
                  error=function(e) NULL)
  if (is.null(map)) {
    map <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=unique(acc), columns="SYMBOL", keytype="ACCNUM"),
                    error=function(e) data.frame(REFSEQ=unique(acc), SYMBOL=NA_character_))
    names(map)[1] <- "REFSEQ"
  }
  map <- map[!duplicated(map$REFSEQ), ]
  sym <- toupper(map$SYMBOL[match(acc, map$REFSEQ)])
} else {
  sym <- toupper(trimws(rownames(M)))
}
valid <- !is.na(sym) & sym != "" & sym != "NA"
M2 <- M[valid, , drop=FALSE]; sym <- sym[valid]
expr_lin <- M2
M2 <- if (is_linear) log2(M2 + 1) else M2

ord  <- order(rowMeans(M2, na.rm=TRUE), decreasing=TRUE); keep <- !duplicated(sym[ord])
df.expr_ann <- M2[ord, , drop=FALSE][keep, , drop=FALSE]; rownames(df.expr_ann) <- sym[ord][keep]
df.expr_ann <- df.expr_ann[order(rownames(df.expr_ann)), , drop=FALSE]
cat(sprintf("Gene-level matrix: %d genes x %d samples\n", nrow(df.expr_ann), ncol(df.expr_ann)))

if (is_refseq) {
  elk1_rows <- which(sub("\\.\\d+$","",rownames(M)) %in% ELK1_REFSEQS)
  if (length(elk1_rows)==0) elk1_rows <- which(sym[match(rownames(M2),rownames(M2))] == "ELK1")
  if (length(elk1_rows) >= 1) {
    e_lin <- if (length(elk1_rows)>1) colMeans(M[elk1_rows, , drop=FALSE], na.rm=TRUE) else as.numeric(M[elk1_rows, ])
    elk1_expr <- if (is_linear) log2(e_lin + 1) else e_lin
  } else elk1_expr <- as.numeric(df.expr_ann["ELK1", ])
} else {
  elk1_expr <- as.numeric(df.expr_ann["ELK1", ])
}
names(elk1_expr) <- colnames(df.expr_ann)
stopifnot("ELK1" %in% rownames(df.expr_ann) || all(is.finite(elk1_expr)))
cat(sprintf("ELK1 vector: range %.2f..%.2f\n", min(elk1_expr,na.rm=TRUE), max(elk1_expr,na.rm=TRUE)))
for (g in c("ELK1", FOCUS3)) cat(sprintf("  [%s] in gene matrix: %s\n", g, ifelse(g %in% rownames(df.expr_ann),"PRESENT","ABSENT")))

## ---- clinical ------------------------------------------------------------------
gc <- function(pat) {
  col <- grep(pat, names(pdata), ignore.case=TRUE, value=TRUE)[1]
  if (is.na(col)) return(rep(NA_character_, nrow(pdata)))
  gsub("^[^:]+:\\s*", "", trimws(as.character(pdata[[col]])))
}
clin <- data.frame(
  MatchID   = gsub(" \\[.*?\\]", "", as.character(pdata$title)),
  GEO       = rownames(pdata),
  os_time   = suppressWarnings(as.numeric(gc("os.?day")))/30.44,
  os_event  = suppressWarnings(as.numeric(gc("os.?bin"))),
  efs_time  = suppressWarnings(as.numeric(gc("efs.?day")))/30.44,
  efs_event = suppressWarnings(as.numeric(gc("efs.?bin"))),
  mycn      = gc("mycn"),
  inss      = gc("inss|stage"),
  high_risk = gc("high.?risk|risk"),
  age_days  = suppressWarnings(as.numeric(gc("age.?at.?diag|age.?diag|age.?day|^age"))),
  death_cause = gc("death.?from.?disease|cause.?of.?death|death.?cause"),
  stringsAsFactors = FALSE)
clin$MYCN      <- ifelse(grepl("amp",clin$mycn,ignore.case=TRUE) & !grepl("non|not|no.?amp",clin$mycn,ignore.case=TRUE),"Amplified",
                         ifelse(grepl("non|not|no.?amp|0",clin$mycn,ignore.case=TRUE),"Not Amplified",NA))
clin$Stage_bin <- ifelse(grepl("4",clin$inss) & !grepl("4s",clin$inss,ignore.case=TRUE),"Stage 4","Non-Stage 4")
.med_age  <- suppressWarnings(median(clin$age_days, na.rm=TRUE))
.age_unit <- if (!is.finite(.med_age)) "unknown" else if (.med_age>400) "days" else if (.med_age>30) "months" else "years"
clin$age_days  <- dplyr::case_when(.age_unit=="days"~clin$age_days, .age_unit=="months"~clin$age_days*30.44,
                                   .age_unit=="years"~clin$age_days*365.25, TRUE~NA_real_)
clin$Age_group <- ifelse(is.na(clin$age_days), NA, ifelse(clin$age_days>=547,">=18mo","<18mo"))
cat(sprintf("Age: unit=%s | raw median=%.1f | >=18mo=%d, <18mo=%d, NA=%d\n", .age_unit, .med_age,
            sum(clin$Age_group==">=18mo",na.rm=TRUE), sum(clin$Age_group=="<18mo",na.rm=TRUE), sum(is.na(clin$Age_group))))
clin$Risk      <- ifelse(grepl("high|1",clin$high_risk,ignore.case=TRUE),"High Risk",
                         ifelse(grepl("low|0",clin$high_risk,ignore.case=TRUE),"Low Risk",NA))

common <- intersect(clin$MatchID, colnames(df.expr_ann))
cat(sprintf("Matched samples (clinical <-> expression): %d / clinical %d / expression %d\n",
            length(common), nrow(clin), ncol(df.expr_ann)))
if (length(common) < 50) {
  cat(">>> WARNING: low match rate. clin$MatchID first 3:", paste(head(clin$MatchID,3),collapse=", "),
      "| expression column first 3:", paste(head(colnames(df.expr_ann),3),collapse=", "), "\n")
}
rownames(clin) <- clin$MatchID
clin <- clin[clin$MatchID %in% common, ]
df.expr_ann <- df.expr_ann[, common, drop=FALSE]
elk1_expr   <- elk1_expr[common]
clin <- clin[common, ]
cat(sprintf("OS available: %d | EFS available: %d / %d\n",
            sum(!is.na(clin$os_time)&!is.na(clin$os_event)), sum(!is.na(clin$efs_time)&!is.na(clin$efs_event)), nrow(clin)))
write.csv(clin, file.path(RES_DIR,"Clinical_data_clean.csv"), row.names=TRUE)

## ---- elk1_groups ------------------------------------------------------------------
samples <- colnames(df.expr_ann); elk1v <- elk1_expr[samples]
if (ELK1_SPLIT=="optimal" && requireNamespace("survminer",quietly=TRUE) &&
    sum(!is.na(clin$efs_time) & !is.na(clin$efs_event)) > 20) {
  cut_df <- data.frame(t=clin[samples,"efs_time"], e=clin[samples,"efs_event"], ELK1=elk1v)
  cut_df <- cut_df[!is.na(cut_df$t)&!is.na(cut_df$e), ]
  sc <- survminer::surv_cutpoint(cut_df, time="t", event="e", variables="ELK1", minprop=0.2)
  cutval <- sc$cutpoint$cutpoint
  elk1_grp <- setNames(ifelse(elk1v>=cutval,"ELK1 High","ELK1 Low"), samples)
  add_fig("00_ELK1_optimal_cutpoint", function() print(plot(sc,"ELK1")), 7, 5)
  cat(sprintf("ELK1 optimal cutpoint = %.3f\n", cutval))
} else if (ELK1_SPLIT=="tertile") {
  qs <- quantile(elk1v, c(1/3,2/3), na.rm=TRUE)
  elk1_grp <- setNames(ifelse(elk1v>=qs[2],"ELK1 High",ifelse(elk1v<=qs[1],"ELK1 Low",NA)), samples)
} else {
  med <- median(elk1v, na.rm=TRUE); elk1_grp <- setNames(ifelse(elk1v>=med,"ELK1 High","ELK1 Low"), samples)
  cat(sprintf("ELK1 median = %.3f\n", med))
}
clin$ELK1_expr <- elk1v[rownames(clin)]; clin$ELK1_group <- elk1_grp[rownames(clin)]

p_dist <- ggplot(data.frame(ELK1=elk1v, Group=factor(elk1_grp[samples],levels=names(PAL_ELK1))), aes(ELK1, fill=Group)) +
  geom_histogram(bins=40, colour="black", linewidth=0.2, alpha=0.85) +
  scale_fill_manual(values=PAL_ELK1, na.value="grey80") + theme_pub() +
  labs(title="ELK1 expression distribution - GSE62564 (RNA-seq)",
       subtitle=sprintf("method: %s - High=%d, Low=%d", ELK1_SPLIT, sum(elk1_grp=="ELK1 High",na.rm=TRUE), sum(elk1_grp=="ELK1 Low",na.rm=TRUE)),
       x="ELK1 (log2)", y="patient count", fill=NULL)
print(p_dist); fig_pdf(p_dist,"01_ELK1_distribution",7,4.5)
cat("\nELK1 group distribution:\n"); print(table(clin$ELK1_group, useNA="ifany"))
res_xlsx(data.frame(PatientID=samples, ELK1_expr=round(elk1v,3), ELK1_group=elk1_grp[samples],
                    MYCN=clin[samples,"MYCN"], Risk=clin[samples,"Risk"], Stage=clin[samples,"Stage_bin"]),
         "ELK1_patient_classification")

if (PRIMARY_COMPARISON=="ELK1") { clin$Group <- clin$ELK1_group; GROUP_LEVELS <- c("ELK1 Low","ELK1 High"); GRP_PAL <- PAL_ELK1
} else if (PRIMARY_COMPARISON=="high_risk") { clin$Group <- clin$Risk; GROUP_LEVELS <- c("Low Risk","High Risk"); GRP_PAL <- PAL_COMP[GROUP_LEVELS]
} else if (PRIMARY_COMPARISON=="mycn") { clin$Group <- ifelse(clin$MYCN=="Amplified","MYCN Amp",ifelse(clin$MYCN=="Not Amplified","MYCN Non-Amp",NA)); GROUP_LEVELS <- c("MYCN Non-Amp","MYCN Amp"); GRP_PAL <- PAL_COMP[GROUP_LEVELS]
} else if (PRIMARY_COMPARISON=="stage4vs4S") { clin$Group <- ifelse(grepl("4s",clin$inss,ignore.case=TRUE),"Stage 4S",ifelse(grepl("4",clin$inss),"Stage 4",NA)); GROUP_LEVELS <- c("Stage 4S","Stage 4"); GRP_PAL <- PAL_COMP[GROUP_LEVELS] }

keep_s <- samples[!is.na(clin[samples,"Group"])]
sample_order <- c(keep_s[clin[keep_s,"Group"]==GROUP_LEVELS[1]], keep_s[clin[keep_s,"Group"]==GROUP_LEVELS[2]])
clin_ord <- clin[sample_order, ]; df.expr_ann <- df.expr_ann[, sample_order, drop=FALSE]
clin_ord$Group <- factor(clin_ord$Group, levels=GROUP_LEVELS)
ann_col_df <- data.frame(Group=clin_ord$Group, row.names=sample_order)
if (!all(is.na(clin_ord$MYCN))) ann_col_df$MYCN <- clin_ord$MYCN
ann_colors_base <- list(Group = GRP_PAL[GROUP_LEVELS])
cat(sprintf("\nComparison: %s | n=%d\n", PRIMARY_COMPARISON, length(sample_order)))

## ---- survival ------------------------------------------------------------------
have_efs <- sum(!is.na(clin$efs_time) & !is.na(clin$efs_event)) > 10

if (have_efs) {
  # Kaplan-Meier (EFS)
  kme <- clin %>% dplyr::filter(!is.na(efs_time), !is.na(efs_event), !is.na(ELK1_group))
  fit_efs <- survminer::surv_fit(survival::Surv(efs_time, efs_event) ~ ELK1_group, data = kme)
  
  p_efs <- survminer::ggsurvplot(
    fit_efs, data = kme, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    palette = unname(PAL_ELK1), title = "Event-Free Survival by ELK1 - GSE62564",
    xlab = "Time (months)", ylab = "EFS probability", ggtheme = theme_pub(11)
  )
  add_fig("03_KM_EFS_ELK1", function() print(p_efs), 8, 8)
  
  # Time-dependent ROC (EFS)
  roc_df <- kme %>% dplyr::filter(!is.na(ELK1_expr))
  if (nrow(roc_df) > 30) {
    tmax <- suppressWarnings(max(roc_df$efs_time[roc_df$efs_event == 1], na.rm = TRUE))
    tt_pts <- c(12, 36, 60); tt_pts <- tt_pts[is.finite(tmax) & tt_pts < tmax]
    
    if (length(tt_pts) >= 1) {
      roc <- tryCatch(
        timeROC::timeROC(
          T = roc_df$efs_time, 
          delta = roc_df$efs_event, 
          marker = roc_df$ELK1_expr, 
          cause = 1, 
          times = tt_pts, 
          iid = FALSE
        ),
        error = function(e) { message("timeROC error: ", conditionMessage(e)); NULL }
      )
      
      if (!is.null(roc)) {
        tcols <- paste0("t=", tt_pts)
        cols  <- c("#4393c3", "#d6604d", "#1a9850")[seq_along(tt_pts)]
        add_fig("07_TimeROC", function() {
          plot(
            roc$FP[, tcols[1]], roc$TP[, tcols[1]], type = "l", col = cols[1], lwd = 2,
            xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
            main = "ELK1 time-dependent ROC (EFS) - GSE62564"
          )
          if (length(tt_pts) > 1) {
            for (i in 2:length(tt_pts)) lines(roc$FP[, tcols[i]], roc$TP[, tcols[i]], col = cols[i], lwd = 2)
          }
          abline(0, 1, lty = 3, col = "grey60")
          legend("bottomright", bty = "n", col = cols, lwd = 2, legend = sprintf("%d mo: AUC=%.2f", tt_pts, roc$AUC[tcols]))
        }, 6, 6)
        cat("Time-ROC (EFS) AUC:\n"); print(round(roc$AUC, 3))
      }
    } else {
      cat("Time points exceed follow-up duration - ROC skipped.\n")
    }
  }
} else {
  cat("EFS data insufficient/missing - survival analysis skipped.\n")
}

## ---- markers_setup ------------------------------------------------------------------
read_markers <- function(file) {
  if (!file.exists(file)) return(character(0))
  m <- read.csv(file, stringsAsFactors=FALSE); gn <- grep("Gene.?Names", names(m), value=TRUE)[1]
  if (is.na(gn)) gn <- names(m)[1]
  m[[gn]] %>% stringr::str_split("\\s+") %>% unlist() %>% toupper() %>% unique() %>% .[. != ""]
}
cd_markers   <- read_markers(file.path(DATA_DIR,"cd_antigen_markers.csv"))
surf_markers <- read_markers(file.path(DATA_DIR,"surface_marker_candidates.csv"))
neuro_genes  <- character(0)

nmk_f <- file.path(DATA_DIR, "neuroblastoma_markers.csv")
if (file.exists(nmk_f)) {
  neuro_mk <- read.csv(nmk_f, stringsAsFactors = FALSE)
  neuro_genes <- c(neuro_mk$Marker, stringr::str_split(neuro_mk$Alias, "[,;/ ]+", simplify = FALSE) %>% unlist()) %>%
    toupper() %>% stringr::str_trim() %>% unique() %>% .[. != "" & !is.na(.)]
}

## Additional neuroblastoma-related genes, merged with neuroblastoma_markers.csv.
## Edit freely: genes from the literature can be added or removed here.
neuro_genes <- unique(c(neuro_genes, toupper(c(
  "PHOX2B","PHOX2A","HAND2","GATA3","TH","DDC","DBH","SLC6A2","PNMT","VIM","FN1","CDH2",
  "TWIST1","ZEB1","ZEB2","SNAI1","SNAI2","NES","SOX10","SOX2","PROM1","NOTCH1","NOTCH2",
  "MYCN","MYC","ODC1","CCND1","CDK4","E2F1","MDM2","BCL2","MCL1","BIRC5","NTRK1","NTRK2",
  "NTRK3","RET","ALK","NGF","BDNF","CD44","CD24","ALDH1A1","ABCG2","BMI1","EZH2","LMO1","ISL1")), FOCUS3))
## Stemness signatures (gene sets used for the Stemness / Neuro+Stem grouping).
## Edit freely: genes from the literature can be added, removed or replaced here,
## and whole signatures can be added. Use official HGNC gene symbols.
stem_sigs <- list(
  Pluripotency = toupper(c("POU5F1","SOX2","NANOG","KLF4","MYC","LIN28A","LIN28B","SALL4","ZFP42","FOXO1","DPPA4","UTF1","PRDM14","DNMT3L")),
  Neural_Stem  = toupper(c("NES","NOTCH1","NOTCH2","SOX2","PAX6","VIM","GFAP","PROM1","OLIG2","ASCL1","HES1","HES5","ID1","ID3")),
  Cancer_Stem  = toupper(c("CD44","PROM1","ALDH1A1","ALDH1A3","ITGA6","ABCG2","BMI1","EZH2","HMGA2","IGF2BP1","IGF2BP2","IGF2BP3","CD24","PCNA","MKI67")),
  MYCN_Pathway = toupper(c("MYCN","ODC1","LDHA","ENO1","CCND1","CDK4","MDM2","ID2","MAX","MXD1","E2F1")))
stem_genes <- unique(unlist(stem_sigs))
group_of <- function(g){ inN<-g%in%neuro_genes; inS<-g%in%stem_genes
ifelse(inN&inS,"Neuro+Stem",ifelse(inN,"Neuroblastoma",ifelse(inS,"Stemness","Other"))) }

## ---- correlation ------------------------------------------------------------------
elk1_full <- as.numeric(df.expr_ann["ELK1", sample_order])
scg <- intersect(FOCUS3, rownames(df.expr_ann))
if (length(scg)>=1) {
  sc <- lapply(scg, function(g) data.frame(Gene=g, ELK1=elk1_full,
                                           Marker=as.numeric(df.expr_ann[g, sample_order]), Group=clin_ord$Group)) %>% dplyr::bind_rows()
  p_sc <- ggplot(sc, aes(ELK1, Marker)) +
    geom_smooth(method="lm", se=TRUE, colour="black", fill="grey85", alpha=0.4, linewidth=0.8) +
    geom_point(aes(colour=Group), size=1.6, alpha=0.6) +
    ggpubr::stat_cor(method="pearson", label.x.npc=0.04, label.y.npc=0.97, size=3.2) +
    facet_wrap(~Gene, nrow=1, scales="free_y") + scale_colour_manual(values=ann_colors_base$Group) +
    theme_pub() + labs(title="ELK1 - CSF3R/CD44/PROM1 - GSE62564", x="ELK1 (log2)", y="Marker (log2)", colour=NULL)
  print(p_sc); fig_pdf(p_sc, "12_ELK1_focus_scatter", 11, 4)
}
fg <- intersect(c("ELK1", FOCUS3), rownames(df.expr_ann))
if (length(fg)>=2) {
  cg <- cor(t(df.expr_ann[fg, sample_order]), method="pearson")
  pheatmap(cg, cluster_rows=TRUE, cluster_cols=TRUE, display_numbers=TRUE, number_format="%.2f",
           color=colorRampPalette(c("#9ecae1","white","#fcbba1"))(100), border_color="grey80",
           main="Gene-gene correlation (ELK1/CSF3R/CD44/PROM1) - GSE62564",
           filename=file.path(FIG_DIR,"13_GeneCorr.png"), width=5, height=5)
  FIGS[["13_GeneCorr"]] <- list(w=5,h=5)
}
panel <- setdiff(intersect(unique(c(FOCUS3, neuro_genes, stem_genes, cd_markers, surf_markers)), rownames(df.expr_ann)), "ELK1")
panel <- panel[apply(df.expr_ann[panel, sample_order, drop=FALSE], 1, function(x) sd(x,na.rm=TRUE)>0)]
cors  <- sapply(panel, function(g) suppressWarnings(cor(elk1_full, as.numeric(df.expr_ann[g, sample_order]), method="pearson")))
corr_tbl <- data.frame(SYMBOL=panel, ELK1_corr=round(cors,3), Group=sapply(panel, group_of)) %>% dplyr::arrange(dplyr::desc(ELK1_corr))
res_xlsx(corr_tbl, "ELK1_correlation_panel")
cat("\nELK1 - focus markers:\n"); print(corr_tbl %>% dplyr::filter(SYMBOL %in% FOCUS3))
top_corr <- corr_tbl %>%
  dplyr::slice(unique(c(head(order(-corr_tbl$ELK1_corr),25), head(order(corr_tbl$ELK1_corr),25)))) %>%
  dplyr::arrange(dplyr::desc(ELK1_corr))
row1  <- matrix(top_corr$ELK1_corr, nrow=1, dimnames=list("ELK1_r", top_corr$SYMBOL))
.lim  <- max(abs(top_corr$ELK1_corr), na.rm=TRUE)
.brks <- seq(-.lim, .lim, length.out = length(hm_pal) + 1)
pheatmap(row1, cluster_rows=FALSE, cluster_cols=FALSE, color=hm_pal, breaks=.brks, border_color=NA,
         fontsize_col=7, angle_col=90, main="ELK1 Pearson r - markers (GSE62564)",
         filename=file.path(FIG_DIR,"14_ELK1_corr_singlerow.png"), width=12, height=2.4)
FIGS[["14_ELK1_corr_singlerow"]] <- list(w=12,h=2.4)

## ---- ELK1 violin plots ------------------------------------------------------
star_fn <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "****", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}
vio_pal <- if (exists("PAL_ELK1")) PAL_ELK1 else setNames(c("#d6604d","#4393c3"), rev(GROUP_LEVELS))

parse_gene_names <- function(file, gene_col = "Gene Names") {
  if (!file.exists(file)) { warning("File not found: ", file); return(character(0)) }
  df  <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  col <- grep(gene_col, names(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(col)) { warning("Column not found: ", gene_col); return(character(0)) }
  tokens <- unlist(strsplit(df[[col]], "\\s+"))
  unique(toupper(trimws(tokens[nzchar(tokens)])))
}
parse_nb_markers <- function(file) {
  if (!file.exists(file)) { warning("File not found: ", file); return(character(0)) }
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  syms <- character(0)
  if ("Alias" %in% names(df)) {
    tok <- unlist(strsplit(df$Alias, "[,;/ ]+")); tok <- toupper(trimws(tok))
    syms <- c(syms, tok[grepl("^[A-Z][A-Z0-9]{1,9}$", tok)])
  }
  if ("Marker" %in% names(df)) syms <- c(syms, toupper(trimws(df$Marker)))
  unique(syms[nzchar(syms)])
}

f_cd   <- file.path(DATA_DIR, "cd_antigen_markers.csv")
f_nb   <- file.path(DATA_DIR, "neuroblastoma_markers.csv")
f_surf <- file.path(DATA_DIR, "surface_marker_candidates.csv")

genes_cd   <- intersect(parse_gene_names(f_cd,   "Gene Names"), rownames(df.expr_ann))
genes_surf <- intersect(parse_gene_names(f_surf, "Gene Names"), rownames(df.expr_ann))
genes_nb   <- intersect(parse_nb_markers(f_nb),               rownames(df.expr_ann))

make_single_violin <- function(gene, grp_var, grp_levels, expr_mat,
                               sample_ord, clin_df, pal) {
  if (!gene %in% rownames(expr_mat)) return(NULL)
  df_g <- data.frame(Sample = sample_ord, Expr = as.numeric(expr_mat[gene, sample_ord]),
                     Group = factor(as.character(clin_df[sample_ord, grp_var]), levels = grp_levels),
                     stringsAsFactors = FALSE) %>% dplyr::filter(!is.na(Group), is.finite(Expr))
  n_tab <- table(df_g$Group)
  n_str <- paste(sapply(grp_levels, function(g) sprintf("%s (n=%d)", g, if (g %in% names(n_tab)) n_tab[[g]] else 0)), collapse = " vs ")
  g1 <- df_g$Expr[df_g$Group == grp_levels[1]]; g2 <- df_g$Expr[df_g$Group == grp_levels[2]]
  p_wilcox <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) NA_real_)
  p_anova  <- tryCatch(summary(aov(Expr ~ Group, data = df_g))[[1]][1,"Pr(>F)"], error = function(e) NA_real_)
  star     <- star_fn(p_wilcox)
  y_max  <- max(df_g$Expr, na.rm = TRUE); y_rng  <- diff(range(df_g$Expr, na.rm = TRUE))
  y_star <- y_max + 0.08 * y_rng; y_seg  <- y_max + 0.04 * y_rng
  ggplot(df_g, aes(x = Group, y = Expr, fill = Group, colour = Group)) +
    geom_violin(trim = TRUE, alpha = 0.35, linewidth = 0.4, scale = "width", colour = NA) +
    geom_jitter(width = 0.18, size = 1.5, alpha = 0.55, shape = 16, stroke = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.35, linewidth = 0.6, colour = "black", fatten = 0) +
    annotate("segment", x = 1, xend = 2, y = y_seg, yend = y_seg, linewidth = 0.5, colour = "grey30") +
    annotate("text", x = 1.5, y = y_star, label = star, size = 5.5, fontface = "bold", colour = "black", vjust = 0) +
    scale_fill_manual(values = pal, breaks = grp_levels) +
    scale_colour_manual(values = pal, breaks = grp_levels) +
    scale_x_discrete(limits = grp_levels,
                     labels = sapply(grp_levels, function(g) sprintf("%s\n(n=%d)", g, if (g %in% names(n_tab)) n_tab[[g]] else 0))) +
    coord_cartesian(ylim = c(min(df_g$Expr, na.rm=TRUE) - 0.05*y_rng, y_star + 0.18 * y_rng)) +
    theme_pub(12) +
    labs(title = sprintf("%s  |  %d patients", gene, nrow(df_g)),
         subtitle = sprintf("%s\nWilcoxon p = %s  |  ANOVA p = %s", n_str,
                            ifelse(is.na(p_wilcox), "NA", format(p_wilcox, digits = 3, scientific = TRUE)),
                            ifelse(is.na(p_anova),  "NA", format(p_anova,  digits = 3, scientific = TRUE))),
         x = PRIMARY_COMPARISON, y = sprintf("log2 of %s", gene), fill = NULL, colour = NULL) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey35"),
          axis.title.x = element_text(face = "bold"))
}

save_violin_pdf <- function(gene_vec, pdf_name, source_label) {
  valid_genes <- intersect(gene_vec, rownames(df.expr_ann))
  valid_samples <- intersect(sample_order, colnames(df.expr_ann))
  if (length(valid_genes) == 0) { message("Violin skipped (no matching genes): ", source_label); return(invisible(NULL)) }
  avail <- valid_genes[apply(df.expr_ann[valid_genes, valid_samples, drop = FALSE], 1, function(x) sd(x, na.rm = TRUE) > 0)]
  if (length(avail) == 0) { message("Violin skipped (no variance): ", source_label); return(invisible(NULL)) }
  cat(sprintf("\n[%s] %d genes -> generating PDF...\n", source_label, length(avail)))
  pdf_path <- file.path(FIG_DIR, paste0(pdf_name, ".pdf"))
  grDevices::pdf(pdf_path, width = 5, height = 6)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  stat_rows <- list()
  for (g in avail) {
    p <- make_single_violin(g, "Group", GROUP_LEVELS, df.expr_ann, valid_samples, clin_ord, vio_pal)
    if (!is.null(p)) {
      print(p)
      df_g <- data.frame(Expr = as.numeric(df.expr_ann[g, valid_samples]),
                         Group = factor(as.character(clin_ord[valid_samples, "Group"]), levels = GROUP_LEVELS)) %>%
        dplyr::filter(!is.na(Group), is.finite(Expr))
      g1 <- df_g$Expr[df_g$Group == GROUP_LEVELS[1]]; g2 <- df_g$Expr[df_g$Group == GROUP_LEVELS[2]]
      p_w <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) NA_real_)
      p_a <- tryCatch(summary(aov(Expr ~ Group, data = df_g))[[1]][1,"Pr(>F)"], error = function(e) NA_real_)
      stat_rows[[g]] <- data.frame(Gene = g, n_total = nrow(df_g), n_group1 = length(g1), n_group2 = length(g2),
                                   median_group1 = round(median(g1, na.rm=TRUE), 3), median_group2 = round(median(g2, na.rm=TRUE), 3),
                                   p_wilcoxon = p_w, p_anova = p_a, star = star_fn(p_w))
    }
  }
  try(grDevices::dev.off(), silent = TRUE)
  cat(sprintf("  -> %s (%d pages)\n", basename(pdf_path), length(avail)))
  if (length(stat_rows) > 0) {
    stat_df <- dplyr::bind_rows(stat_rows) %>% dplyr::arrange(p_wilcoxon)
    res_xlsx(stat_df, paste0(pdf_name, "_stats"))
  }
  invisible(pdf_path)
}

save_violin_pdf(genes_cd,   "25_CD_antigen_violin_per_gene",    "cd_antigen_markers")
save_violin_pdf(genes_surf, "26_Surface_markers_violin_per_gene", "surface_marker_candidates")
save_violin_pdf(genes_nb,   "27_Neuroblastoma_CSV_violin_per_gene", "neuroblastoma_markers")
genes_all <- unique(c(genes_cd, genes_surf, genes_nb))
save_violin_pdf(genes_all,  "28_All_CSV_violin_per_gene", "all CSV combined")

preview_genes <- intersect(c("PROM1","CD44","CSF3R"), rownames(df.expr_ann))
if (length(preview_genes) > 0) {
  for (g in preview_genes) {
    p <- make_single_violin(g, "Group", GROUP_LEVELS, df.expr_ann, sample_order, clin_ord, vio_pal)
    if (!is.null(p)) { print(p); fig_pdf(p, paste0("29_preview_", g), w = 5, h = 6) }
  }
}


## ################################################################################
## SUPPLEMENTARY - MYCN expression across cell lines (MYCN status only, no TF overlay)
## ################################################################################

DATA_DIR_CELL <- "path/to/your/cell_line_data"
MYCN_FILE <- file.path(DATA_DIR_CELL, "2016-11-17-CellLine-MYCN-status.txt")
PAL_MYCN  <- c(Amplified = "#1a7a1a", Nonamplified = "#8ecae6")

stopifnot(exists("tpm_raw"), "MYCN" %in% rownames(tpm_raw))

mycn_status <- read.delim(MYCN_FILE, stringsAsFactors = FALSE, check.names = FALSE)
names(mycn_status) <- c("CellLine", "Status")
mycn_status$CellLine <- toupper(gsub("[^A-Za-z0-9]", "", mycn_status$CellLine))
.norm_cell2 <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

mycn_fpkm <- data.frame(
  CellLineNorm = .norm_cell2(colnames(tpm_raw)),
  CellLine     = colnames(tpm_raw),
  MYCN         = as.numeric(tpm_raw["MYCN", ]),
  stringsAsFactors = FALSE
)
mycn_fpkm <- merge(mycn_fpkm, mycn_status, by.x = "CellLineNorm", by.y = "CellLine", all.x = TRUE)
mycn_fpkm$Status[is.na(mycn_fpkm$Status)] <- "Unknown"
mycn_fpkm$Status <- factor(mycn_fpkm$Status, levels = c("Amplified", "Nonamplified", "Unknown"))
mycn_fpkm <- mycn_fpkm %>% dplyr::arrange(MYCN) %>% dplyr::mutate(CellLine = factor(CellLine, levels = CellLine))
pal_mycn_full <- c(PAL_MYCN, Unknown = "grey70")

p_mycn <- ggplot(mycn_fpkm, aes(x = CellLine, y = MYCN, fill = Status)) +
  geom_col(colour = NA, width = 0.75) +
  scale_fill_manual(values = pal_mycn_full, name = "Status") +
  labs(title = "MYCN Expression Across Cell Lines", x = "Cell Line", y = "FPKM") +
  theme_pub(11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
print(p_mycn)
fig_pdf(p_mycn, "30_MYCN_status", 14, 7)
res_xlsx(mycn_fpkm %>% dplyr::select(CellLine, MYCN, Status), "MYCN_by_status")


############################################################################
# PEA3 FAMILY (ETV1, ETV4, ETV5) - each gene analysed independently
# (no averaged score): classification, EFS, ROC, correlation and violin
# plots are produced separately for every gene.
############################################################################

## ---- 0. SETTINGS ------------------------------------------------------------
## Set the paths below to your local directories.
suppressWarnings(suppressMessages({
  
  DATA_DIR_CELL <- "path/to/your/pea3_data"
  DATA_DIR_PAT  <- "path/to/your/pea3_data"
  OUT_DIR       <- "path/to/your/pea3_output"
  
  FIG_DIR <- file.path(OUT_DIR, "Figures")
  RES_DIR <- file.path(OUT_DIR, "Results")
  CACHE   <- file.path(OUT_DIR, "geo_cache")
  SUPP    <- file.path(OUT_DIR, "GSE62564_Supp")
  for (d in c(OUT_DIR, FIG_DIR, RES_DIR, CACHE, SUPP)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  
  CELLLINE_FILE <- file.path(DATA_DIR_CELL, "GSE89413_2016-10-30-NBL-cell-line-STAR-fpkm.txt")
  MYCN_FILE     <- file.path(DATA_DIR_CELL, "2016-11-17-CellLine-MYCN-status.txt")
  
  PEA3_GENES <- c("ETV1", "ETV4", "ETV5")  # PEA3 sub-family, each analysed separately
  FOCUS3     <- c("CSF3R", "CD44", "PROM1")
  ELK1_SPLIT <- "optimal"                  # split method: "optimal" | "median" | "tertile" (also applied to each PEA3 gene)
  
  PAL_MYCN <- c(Amplified = "#1a7a1a", Nonamplified = "#8ecae6")
  hm_pal   <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  
  # Per-gene palette (High = red, Low = blue)
  pal_for_gene <- function(gene) setNames(c("#d6604d", "#4393c3"), c(paste(gene, "High"), paste(gene, "Low")))
  
  theme_pub <- function(base_size = 13) {
    ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.35),
        axis.text  = ggplot2::element_text(colour = "black", size = base_size - 1),
        axis.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 1),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
        legend.key = ggplot2::element_blank())
  }
  
  open_pdf <- function(path, w, h) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); grDevices::pdf(path, width = w, height = h) }
  open_png <- function(path, w, h) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); grDevices::png(path, width = round(w*300), height = round(h*300), res = 300) }
  .draw_to <- function(opener, draw) {
    ok <- tryCatch({ opener(); TRUE }, error = function(e) { message("device could not be opened: ", conditionMessage(e)); FALSE })
    if (!ok) return(FALSE)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    isTRUE(tryCatch({ draw(); TRUE }, error = function(e) { message("drawing error: ", conditionMessage(e)); FALSE }))
  }
  add_fig <- function(name, draw, w = 8, h = 6) {
    fp <- file.path(FIG_DIR, paste0(name, ".pdf")); pp <- file.path(FIG_DIR, paste0(name, ".png"))
    .draw_to(function() open_pdf(fp, w, h), draw); .draw_to(function() open_png(pp, w, h), draw)
    got <- c(if (file.exists(fp)) "pdf", if (file.exists(pp)) "png")
    cat("  ", name, " -> ", paste(got, collapse = "+"), if (length(got) < 2) "  (MISSING!)" else "", "\n", sep = "")
    invisible(name)
  }
  fig_pdf <- function(plot, name, w = 8, h = 6) add_fig(name, function() print(plot), w, h)
  res_xlsx <- function(x, name) {
    f <- file.path(RES_DIR, paste0(name, ".xlsx")); writexl::write_xlsx(as.data.frame(x), f)
    if (file.exists(f)) cat("  -> table:", f, "\n"); invisible(f)
  }
  star_fn <- function(p) {
    dplyr::case_when(is.na(p) ~ "ns", p < 0.0001 ~ "****", p < 0.001 ~ "***",
                     p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
  }
  read_markers <- function(file) {
    if (!file.exists(file)) { message("Marker file not found, skipping: ", file); return(character(0)) }
    m <- read.csv(file, stringsAsFactors = FALSE)
    gn <- grep("Gene.?Names", names(m), value = TRUE)[1]; if (is.na(gn)) gn <- names(m)[1]
    m[[gn]] %>% stringr::str_split("\\s+") %>% unlist() %>% toupper() %>% unique() %>% .[. != ""]
  }
  
}))



## ############################################################################
## PEA3 - CELL LINES (GSE89413): PEA3 genes, surface heatmap, MYCN bar
## ############################################################################

## A1 - Load FPKM matrix --------------------------------------------------
read_tpm <- function(path) {
  stopifnot(file.exists(path))
  l1 <- readLines(path, n = 1, warn = FALSE)
  delim <- if (grepl("\t", l1)) "\t" else if (grepl(";", l1)) ";" else if (grepl(",", l1)) "," else ""
  df <- read.delim(path, sep = delim, header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  is_num_col <- vapply(df, function(c) mean(!is.na(suppressWarnings(as.numeric(c)))) > 0.9, logical(1))
  id_cols <- which(!is_num_col); num_cols <- which(is_num_col)
  if (length(num_cols) < 3) stop("Numeric column not found - sep might be wrong. l1: ", l1)
  sym_idx <- id_cols[ vapply(df[id_cols], function(c) any(grepl("[A-Za-z]", c)), logical(1)) ]
  sym_idx <- if (length(sym_idx)) sym_idx[1] else id_cols[1]
  sym <- toupper(trimws(as.character(df[[sym_idx]])))
  Mm <- as.matrix(df[, num_cols, drop = FALSE]); mode(Mm) <- "numeric"; rownames(Mm) <- sym
  if (any(duplicated(sym))) {
    ord <- order(rowMeans(Mm, na.rm = TRUE), decreasing = TRUE)
    Mm <- Mm[ord, , drop = FALSE]; Mm <- Mm[!duplicated(rownames(Mm)), , drop = FALSE]
  }
  Mm[rownames(Mm) != "" & !is.na(rownames(Mm)), , drop = FALSE]
}

tpm_raw <- read_tpm(CELLLINE_FILE)
rng <- range(tpm_raw, na.rm = TRUE)
is_linear <- rng[2] > 30 || any(tpm_raw > 30, na.rm = TRUE)
expr <- if (is_linear) log2(tpm_raw + 1) else tpm_raw
cat(sprintf("Cell line matrix: %d genes x %d lines | scale: %s\n",
            nrow(expr), ncol(expr), ifelse(is_linear, "log2(FPKM+1)", "already log")))
for (g in PEA3_GENES) if (!g %in% rownames(expr)) message("WARNING: ", g, " not in matrix.")

## A2 - Surface marker list ----------------------------------------------
surf_markers <- read_markers(file.path(DATA_DIR_CELL, "surface_marker_candidates.csv"))

## A3 - Cell line bar plot: one plot per PEA3 gene (no averaging) -----------
for (gene in PEA3_GENES) {
  if (!gene %in% rownames(expr)) { message("Skipped (not in matrix): ", gene); next }
  df_g <- data.frame(Cell = colnames(expr), Expr = as.numeric(expr[gene, ])) %>%
    dplyr::arrange(Expr) %>%
    dplyr::mutate(Cell = factor(Cell, levels = Cell))
  p_g <- ggplot(df_g, aes(Expr, Cell)) +
    geom_col(width = 0.75, fill = "#4393c3", colour = "black", linewidth = 0.2) +
    labs(title = sprintf("%s expression - neuroblastoma cell lines", gene),
         subtitle = "log2(FPKM+1)", x = sprintf("%s (log2 FPKM)", gene), y = NULL) +
    theme_pub(11)
  print(p_g)
  fig_pdf(p_g, sprintf("01_%s_across_celllines", gene), 7, max(6, ncol(expr) * 0.22))
  res_xlsx(df_g, sprintf("%s_expression_per_cellline", gene))
}

## A4 - Surface marker heatmap (compact: row names hidden, fixed 7x6 in) ---
surf_in_expr <- intersect(surf_markers, rownames(expr))
if (length(surf_in_expr) >= 3) {
  vg_surf  <- apply(expr[surf_in_expr, , drop = FALSE], 1, var, na.rm = TRUE)
  top_surf <- names(sort(vg_surf, decreasing = TRUE))[seq_len(min(500, length(vg_surf)))]
  g_hm <- unique(c(intersect(PEA3_GENES, rownames(expr)), top_surf))
  s <- expr[g_hm, , drop = FALSE]
  s <- s[apply(s, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]
  pheatmap(s, scale = "row", color = hm_pal, border_color = NA,
           cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D2",
           show_rownames = FALSE, fontsize_col = 8,
           main = "Surface Markers - top 500 most variable (z-score) + ETV1/ETV4/ETV5",
           filename = file.path(FIG_DIR, "02_Surface_heatmap.png"), width = 7, height = 6)
  cat("  02_Surface_heatmap -> png (compact)\n")
  res_xlsx(data.frame(Gene = rownames(s), s, check.names = FALSE), "Surface_heatmap_matrix")
} else {
  message("Surface marker list not found/too few - heatmap skipped. Expected file: ",
          file.path(DATA_DIR_CELL, "surface_marker_candidates.csv"))
}

## A5 - PEA3 genes (one facet each) + focus markers ------------------------
combo_genes <- intersect(c(PEA3_GENES, FOCUS3), rownames(expr))
cb_long <- as.data.frame(expr[combo_genes, , drop = FALSE]) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Cell", values_to = "Expr") %>%
  dplyr::mutate(Gene = factor(Gene, levels = combo_genes))

p_cb <- ggplot(cb_long, aes(Expr, Cell)) +
  geom_col(width = 0.78, fill = "#74a9cf", colour = "black", linewidth = 0.12) +
  facet_wrap(~ Gene, nrow = 1, scales = "free_x") +
  labs(title = "ETV1 / ETV4 / ETV5 (separately) + CSF3R / CD44 / PROM1 - cell lines",
       subtitle = "log2(FPKM+1)", x = "log2 FPKM", y = NULL) +
  theme_pub(10)
print(p_cb)
fig_pdf(p_cb, "03_PEA3_focus_across_cells", 13, max(6, ncol(expr) * 0.22))

## A6 - MYCN bar plot (MYCN expression only, no overlay) --------------------
mycn_status <- read.delim(MYCN_FILE, stringsAsFactors = FALSE, check.names = FALSE)
names(mycn_status) <- c("CellLine", "Status")
mycn_status$CellLine <- toupper(gsub("[^A-Za-z0-9]", "", mycn_status$CellLine))
.norm_cell <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

stopifnot("MYCN" %in% rownames(tpm_raw))

mycn_fpkm <- data.frame(
  CellLineNorm = .norm_cell(colnames(tpm_raw)),
  CellLine     = colnames(tpm_raw),
  MYCN         = as.numeric(tpm_raw["MYCN", ]),
  stringsAsFactors = FALSE
)
mycn_fpkm <- merge(mycn_fpkm, mycn_status, by.x = "CellLineNorm", by.y = "CellLine", all.x = TRUE)
mycn_fpkm$Status[is.na(mycn_fpkm$Status)] <- "Unknown"
mycn_fpkm$Status <- factor(mycn_fpkm$Status, levels = c("Amplified", "Nonamplified", "Unknown"))
mycn_fpkm <- mycn_fpkm %>% dplyr::arrange(MYCN) %>% dplyr::mutate(CellLine = factor(CellLine, levels = CellLine))
pal_mycn_full <- c(PAL_MYCN, Unknown = "grey70")

p_mycn <- ggplot(mycn_fpkm, aes(x = CellLine, y = MYCN, fill = Status)) +
  geom_col(colour = NA, width = 0.75) +
  scale_fill_manual(values = pal_mycn_full, name = "Status") +
  labs(title = "MYCN Expression Across Cell Lines", x = "Cell Line", y = "FPKM") +
  theme_pub(11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
print(p_mycn)
fig_pdf(p_mycn, "04_MYCN_status", 14, 7)
res_xlsx(mycn_fpkm %>% dplyr::select(CellLine, MYCN, Status), "MYCN_by_status_PEA3section")


## ############################################################################
## PEA3 - PATIENTS (GSE62564, n~498): independent classification + EFS + ROC +
## correlation + violin for each PEA3 gene
## ############################################################################

## B1 - Download from GEO (clinical data + supplementary RNA-seq matrix) ---
cat("[GSE62564] downloading clinical data...\n")
gse <- GEOquery::getGEO("GSE62564", GSEMatrix = TRUE, getGPL = FALSE, destdir = CACHE)
if (is.list(gse)) gse <- gse[[1]]
pdata <- Biobase::pData(gse)

cat("[GSE62564] downloading supplementary RNA-seq expression matrix...\n")
supp_paths <- GEOquery::getGEOSuppFiles("GSE62564", baseDir = SUPP, makeDirectory = FALSE)
expr_file  <- rownames(supp_paths)[grepl("txt.gz|tsv.gz|csv.gz", rownames(supp_paths), ignore.case = TRUE)][1]

l1    <- readLines(expr_file, n = 1, warn = FALSE)
delim <- if (grepl("\t", l1)) "\t" else if (grepl(";", l1)) ";" else ","
raw   <- read.delim(expr_file, sep = delim, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
id_raw <- as.character(raw[[1]])
Mpat   <- as.matrix(raw[, -1, drop = FALSE]); mode(Mpat) <- "numeric"; rownames(Mpat) <- id_raw
is_refseq_pat <- mean(grepl("^N[MR]_", id_raw)) > 0.5
rng2 <- range(Mpat, na.rm = TRUE)
is_linear_pat <- rng2[2] > 30 || any(Mpat > 30, na.rm = TRUE)

## B2 - RefSeq -> gene symbol (gene-level matrix) ---------------------------
if (is_refseq_pat) {
  acc <- sub("\\.\\d+$", "", rownames(Mpat))
  map <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = unique(acc), columns = "SYMBOL", keytype = "REFSEQ"),
                  error = function(e) NULL)
  if (is.null(map)) {
    map <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = unique(acc), columns = "SYMBOL", keytype = "ACCNUM"),
                    error = function(e) data.frame(REFSEQ = unique(acc), SYMBOL = NA_character_))
    names(map)[1] <- "REFSEQ"
  }
  map <- map[!duplicated(map$REFSEQ), ]
  sym <- toupper(map$SYMBOL[match(acc, map$REFSEQ)])
} else {
  sym <- toupper(trimws(rownames(Mpat)))
}
valid <- !is.na(sym) & sym != "" & sym != "NA"
M2 <- Mpat[valid, , drop = FALSE]; sym <- sym[valid]
M2 <- if (is_linear_pat) log2(M2 + 1) else M2
ord  <- order(rowMeans(M2, na.rm = TRUE), decreasing = TRUE); keep <- !duplicated(sym[ord])
df.expr_ann <- M2[ord, , drop = FALSE][keep, , drop = FALSE]; rownames(df.expr_ann) <- sym[ord][keep]
df.expr_ann <- df.expr_ann[order(rownames(df.expr_ann)), , drop = FALSE]
for (g in PEA3_GENES) if (!g %in% rownames(df.expr_ann)) message("WARNING: ", g, " not in patient matrix.")
cat(sprintf("[GSE62564] Gene-level matrix: %d genes x %d patients\n", nrow(df.expr_ann), ncol(df.expr_ann)))

## B3 - Clean clinical data and match to expression samples ------------------
gc_ <- function(pat) {
  col <- grep(pat, names(pdata), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(col)) return(rep(NA_character_, nrow(pdata)))
  gsub("^[^:]+:\\s*", "", trimws(as.character(pdata[[col]])))
}
clin <- data.frame(
  MatchID   = gsub(" \\[.*?\\]", "", as.character(pdata$title)),
  GEO       = rownames(pdata),
  os_time   = suppressWarnings(as.numeric(gc_("os.?day"))) / 30.44,
  os_event  = suppressWarnings(as.numeric(gc_("os.?bin"))),
  efs_time  = suppressWarnings(as.numeric(gc_("efs.?day"))) / 30.44,
  efs_event = suppressWarnings(as.numeric(gc_("efs.?bin"))),
  mycn      = gc_("mycn"),
  inss      = gc_("inss|stage"),
  high_risk = gc_("high.?risk|risk"),
  stringsAsFactors = FALSE
)
clin$MYCN <- ifelse(grepl("amp", clin$mycn, ignore.case = TRUE) & !grepl("non|not|no.?amp", clin$mycn, ignore.case = TRUE),
                    "Amplified", ifelse(grepl("non|not|no.?amp|0", clin$mycn, ignore.case = TRUE), "Not Amplified", NA))
rownames(clin) <- clin$MatchID
common <- intersect(clin$MatchID, colnames(df.expr_ann))
cat(sprintf("Matched samples (clinical <-> expression): %d\n", length(common)))
clin <- clin[common, ]
df.expr_ann <- df.expr_ann[, common, drop = FALSE]

## B4 - Violin plot function (one gene per call) ----------------------------
make_single_violin <- function(gene, grp_var, grp_levels, expr_mat, sample_ord, clin_df, pal) {
  if (!gene %in% rownames(expr_mat)) return(NULL)
  df_g <- data.frame(
    Sample = sample_ord, Expr = as.numeric(expr_mat[gene, sample_ord]),
    Group = factor(as.character(clin_df[sample_ord, grp_var]), levels = grp_levels),
    stringsAsFactors = FALSE
  ) %>% dplyr::filter(!is.na(Group), is.finite(Expr))
  n_tab <- table(df_g$Group)
  g1 <- df_g$Expr[df_g$Group == grp_levels[1]]; g2 <- df_g$Expr[df_g$Group == grp_levels[2]]
  p_wilcox <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) NA_real_)
  star <- star_fn(p_wilcox)
  y_max <- max(df_g$Expr, na.rm = TRUE); y_rng <- diff(range(df_g$Expr, na.rm = TRUE))
  y_star <- y_max + 0.08 * y_rng; y_seg <- y_max + 0.04 * y_rng
  ggplot(df_g, aes(x = Group, y = Expr, fill = Group, colour = Group)) +
    geom_violin(trim = TRUE, alpha = 0.35, linewidth = 0.4, scale = "width", colour = NA) +
    geom_jitter(width = 0.18, size = 1.5, alpha = 0.55, shape = 16, stroke = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.35, linewidth = 0.6, colour = "black", fatten = 0) +
    annotate("segment", x = 1, xend = 2, y = y_seg, yend = y_seg, linewidth = 0.5, colour = "grey30") +
    annotate("text", x = 1.5, y = y_star, label = star, size = 5.5, fontface = "bold", colour = "black", vjust = 0) +
    scale_fill_manual(values = pal, breaks = grp_levels) +
    scale_colour_manual(values = pal, breaks = grp_levels) +
    scale_x_discrete(limits = grp_levels,
                     labels = sapply(grp_levels, function(g) sprintf("%s\n(n=%d)", g, if (g %in% names(n_tab)) n_tab[[g]] else 0))) +
    coord_cartesian(ylim = c(min(df_g$Expr, na.rm = TRUE) - 0.05 * y_rng, y_star + 0.18 * y_rng)) +
    theme_pub(12) +
    labs(title = sprintf("%s  |  %d patients", gene, nrow(df_g)),
         subtitle = sprintf("Wilcoxon p = %s", ifelse(is.na(p_wilcox), "NA", format(p_wilcox, digits = 3, scientific = TRUE))),
         x = sprintf("%s Group", gene), y = sprintf("log2 of %s", gene), fill = NULL, colour = NULL) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey35"))
}

## Writes one PNG file per gene.
save_violin_png <- function(gene_vec, base_name, source_label, samp_order, grp_levels, clin_df, pal) {
  valid_genes <- intersect(gene_vec, rownames(df.expr_ann))
  if (length(valid_genes) == 0) { message("Violin skipped (no matching genes): ", source_label); return(invisible(NULL)) }
  avail <- valid_genes[apply(df.expr_ann[valid_genes, samp_order, drop = FALSE], 1, function(x) sd(x, na.rm = TRUE) > 0)]
  if (length(avail) == 0) { message("Violin skipped (no variance): ", source_label); return(invisible(NULL)) }
  cat(sprintf("\n[%s] %d genes -> generating PNG...\n", source_label, length(avail)))
  stat_rows <- list()
  for (g in avail) {
    p <- make_single_violin(g, "Group", grp_levels, df.expr_ann, samp_order, clin_df, pal)
    if (!is.null(p)) {
      png_path <- file.path(FIG_DIR, paste0(base_name, "_", g, ".png"))
      dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
      grDevices::png(png_path, width = round(5 * 300), height = round(6 * 300), res = 300)
      print(p)
      try(grDevices::dev.off(), silent = TRUE)
      cat(sprintf("  -> %s\n", basename(png_path)))
      g1 <- as.numeric(df.expr_ann[g, samp_order])[clin_df[samp_order, "Group"] == grp_levels[1]]
      g2 <- as.numeric(df.expr_ann[g, samp_order])[clin_df[samp_order, "Group"] == grp_levels[2]]
      p_w <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) NA_real_)
      stat_rows[[g]] <- data.frame(Gene = g, n_low = length(g1), n_high = length(g2),
                                   median_low = round(median(g1, na.rm = TRUE), 3),
                                   median_high = round(median(g2, na.rm = TRUE), 3),
                                   p_wilcoxon = p_w, star = star_fn(p_w))
    }
  }
  if (length(stat_rows) > 0) res_xlsx(dplyr::bind_rows(stat_rows) %>% dplyr::arrange(p_wilcoxon), paste0(base_name, "_stats"))
}

surf_markers_pat <- read_markers(file.path(DATA_DIR_PAT, "surface_marker_candidates.csv"))
cd_markers_pat   <- read_markers(file.path(DATA_DIR_PAT, "cd_antigen_markers.csv"))

## B5 - Independent full analysis per PEA3 gene
## (classification / EFS / ROC / correlation / violin)
run_pea3_gene_analysis <- function(gene) {
  
  cat(sprintf("\n=========================== %s ===========================\n", gene))
  if (!gene %in% rownames(df.expr_ann)) { message(gene, " not in patient matrix - skipped."); return(invisible(NULL)) }
  
  samples <- colnames(df.expr_ann)
  gene_v  <- as.numeric(df.expr_ann[gene, samples]); names(gene_v) <- samples
  pal_g   <- pal_for_gene(gene)
  grp_levels <- c(paste(gene, "Low"), paste(gene, "High"))
  
  ## -- classification (by the gene's own expression, no averaging) --
  if (ELK1_SPLIT == "optimal" && sum(!is.na(clin$efs_time) & !is.na(clin$efs_event)) > 20) {
    cut_df <- data.frame(t = clin[samples, "efs_time"], e = clin[samples, "efs_event"], G = gene_v)
    cut_df <- cut_df[!is.na(cut_df$t) & !is.na(cut_df$e), ]
    sc <- tryCatch(survminer::surv_cutpoint(cut_df, time = "t", event = "e", variables = "G", minprop = 0.2),
                   error = function(e) NULL)
    if (!is.null(sc)) {
      cutval <- sc$cutpoint$cutpoint
      gene_grp <- setNames(ifelse(gene_v >= cutval, grp_levels[2], grp_levels[1]), samples)
      add_fig(sprintf("05_%s_optimal_cutpoint", gene), function() print(plot(sc, "G")), 7, 5)
      cat(sprintf("%s optimal cutpoint = %.3f\n", gene, cutval))
    } else {
      med <- median(gene_v, na.rm = TRUE)
      gene_grp <- setNames(ifelse(gene_v >= med, grp_levels[2], grp_levels[1]), samples)
    }
  } else if (ELK1_SPLIT == "tertile") {
    qs <- quantile(gene_v, c(1/3, 2/3), na.rm = TRUE)
    gene_grp <- setNames(ifelse(gene_v >= qs[2], grp_levels[2], ifelse(gene_v <= qs[1], grp_levels[1], NA)), samples)
  } else {
    med <- median(gene_v, na.rm = TRUE)
    gene_grp <- setNames(ifelse(gene_v >= med, grp_levels[2], grp_levels[1]), samples)
    cat(sprintf("%s median = %.3f\n", gene, med))
  }
  
  clin_g <- clin
  clin_g[[paste0(gene, "_expr")]]  <- gene_v[rownames(clin_g)]
  clin_g[[paste0(gene, "_group")]] <- gene_grp[rownames(clin_g)]
  
  p_dist <- ggplot(data.frame(Expr = gene_v, Group = factor(gene_grp[samples], levels = grp_levels)),
                   aes(Expr, fill = Group)) +
    geom_histogram(bins = 40, colour = "black", linewidth = 0.2, alpha = 0.85) +
    scale_fill_manual(values = pal_g, na.value = "grey80") + theme_pub() +
    labs(title = sprintf("%s expression distribution - GSE62564 (RNA-seq)", gene),
         subtitle = sprintf("method: %s | High=%d, Low=%d", ELK1_SPLIT,
                            sum(gene_grp == grp_levels[2], na.rm = TRUE), sum(gene_grp == grp_levels[1], na.rm = TRUE)),
         x = sprintf("%s (log2)", gene), y = "patient count", fill = NULL)
  print(p_dist)
  fig_pdf(p_dist, sprintf("06_%s_distribution", gene), 7, 4.5)
  cat(sprintf("%s group distribution:\n", gene)); print(table(gene_grp, useNA = "ifany"))
  res_xlsx(data.frame(PatientID = samples, Expr = round(gene_v, 3), Group = gene_grp[samples],
                      MYCN = clin_g[samples, "MYCN"]), sprintf("%s_patient_classification", gene))
  
  keep_s <- samples[!is.na(clin_g[samples, paste0(gene, "_group")])]
  samp_order_g <- c(keep_s[clin_g[keep_s, paste0(gene, "_group")] == grp_levels[1]],
                    keep_s[clin_g[keep_s, paste0(gene, "_group")] == grp_levels[2]])
  clin_ord_g <- clin_g[samp_order_g, ]
  clin_ord_g$Group <- factor(clin_ord_g[[paste0(gene, "_group")]], levels = grp_levels)
  cat(sprintf("Comparison: %s High vs Low | n=%d\n", gene, length(samp_order_g)))
  
  ## -- Event-free survival (Kaplan-Meier) --
  if (sum(!is.na(clin_g$efs_time) & !is.na(clin_g$efs_event)) > 10) {
    kme <- clin_g %>% dplyr::filter(!is.na(efs_time), !is.na(efs_event), !is.na(.data[[paste0(gene,"_group")]]))
    kme$Group <- factor(kme[[paste0(gene, "_group")]], levels = grp_levels)
    fit_efs <- survminer::surv_fit(survival::Surv(efs_time, efs_event) ~ Group, data = kme)
    p_efs <- survminer::ggsurvplot(fit_efs, data = kme, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                                   palette = unname(pal_g[grp_levels]),
                                   title = sprintf("Event-Free Survival by %s - GSE62564", gene),
                                   xlab = "Time (months)", ylab = "EFS probability", ggtheme = theme_pub(11))
    print(p_efs)
    add_fig(sprintf("07_KM_EFS_%s", gene), function() print(p_efs), 8, 8)
  } else {
    message("EFS data insufficient - KM skipped (", gene, ").")
  }
  
  ## -- Time-dependent ROC (EFS) --
  roc_df <- clin_g %>% dplyr::filter(!is.na(.data[[paste0(gene,"_expr")]]), !is.na(efs_time), !is.na(efs_event))
  if (nrow(roc_df) > 30) {
    tmax <- suppressWarnings(max(roc_df$efs_time[roc_df$efs_event == 1], na.rm = TRUE))
    tt_pts <- c(12, 36, 60); tt_pts <- tt_pts[is.finite(tmax) & tt_pts < tmax]
    if (length(tt_pts) >= 1) {
      roc <- tryCatch(
        timeROC::timeROC(T = roc_df$efs_time, delta = roc_df$efs_event, marker = roc_df[[paste0(gene,"_expr")]],
                         cause = 1, times = tt_pts, iid = FALSE),
        error = function(e) { message("timeROC error: ", conditionMessage(e)); NULL })
      if (!is.null(roc)) {
        tcols <- paste0("t=", tt_pts)
        cols  <- c("#4393c3", "#d6604d", "#1a9850")[seq_along(tt_pts)]
        add_fig(sprintf("08_TimeROC_EFS_%s", gene), function() {
          plot(roc$FP[, tcols[1]], roc$TP[, tcols[1]], type = "l", col = cols[1], lwd = 2,
               xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
               main = sprintf("%s time-dependent ROC (EFS) - GSE62564", gene))
          if (length(tt_pts) > 1) for (i in 2:length(tt_pts)) lines(roc$FP[, tcols[i]], roc$TP[, tcols[i]], col = cols[i], lwd = 2)
          abline(0, 1, lty = 3, col = "grey60")
          legend("bottomright", bty = "n", col = cols, lwd = 2, legend = sprintf("%d mo: AUC=%.2f", tt_pts, roc$AUC[tcols]))
        }, 6, 6)
        cat(sprintf("Time-ROC (EFS) AUC - %s:\n", gene)); print(round(roc$AUC, 3))
      }
    } else message("Time points exceed follow-up duration - ROC skipped (", gene, ").")
  } else message("EFS + ", gene, " data insufficient - ROC skipped.")
  
  ## -- Correlation (gene vs CSF3R/CD44/PROM1 + wider marker panel) --
  gene_full <- as.numeric(df.expr_ann[gene, samp_order_g])
  scg <- intersect(FOCUS3, rownames(df.expr_ann))
  if (length(scg) >= 1) {
    sc <- lapply(scg, function(mk) data.frame(Gene = mk, TF = gene_full,
                                              Marker = as.numeric(df.expr_ann[mk, samp_order_g]), Group = clin_ord_g$Group)) %>% dplyr::bind_rows()
    p_sc <- ggplot(sc, aes(TF, Marker)) +
      geom_smooth(method = "lm", se = TRUE, colour = "black", fill = "grey85", linewidth = 0.8) +
      geom_point(aes(colour = Group), size = 1.6, alpha = 0.6) +
      ggpubr::stat_cor(method = "pearson", label.x.npc = 0.04, label.y.npc = 0.97, size = 3.2) +
      facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
      scale_colour_manual(values = pal_g) +
      theme_pub() + labs(title = sprintf("%s - CSF3R/CD44/PROM1 - GSE62564", gene),
                         x = sprintf("%s (log2)", gene), y = "Marker (log2)", colour = NULL)
    print(p_sc)
    fig_pdf(p_sc, sprintf("09_%s_focus_scatter", gene), 11, 4)
  }
  fg <- intersect(c(gene, FOCUS3), rownames(df.expr_ann))
  if (length(fg) >= 2) {
    cg <- cor(t(df.expr_ann[fg, samp_order_g]), method = "pearson")
    pheatmap(cg, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE, number_format = "%.2f",
             color = colorRampPalette(c("#9ecae1", "white", "#fcbba1"))(100), border_color = "grey80",
             main = sprintf("Gene-gene correlation (%s/CSF3R/CD44/PROM1) - GSE62564", gene),
             filename = file.path(FIG_DIR, sprintf("10_GeneCorr_%s.png", gene)), width = 6, height = 6)
    cat(sprintf("  10_GeneCorr_%s -> png\n", gene))
  }
  panel <- setdiff(intersect(unique(c(FOCUS3, surf_markers_pat, cd_markers_pat)), rownames(df.expr_ann)), gene)
  if (length(panel) >= 3) {
    panel <- panel[apply(df.expr_ann[panel, samp_order_g, drop = FALSE], 1, function(x) sd(x, na.rm = TRUE) > 0)]
    cors  <- sapply(panel, function(mk) suppressWarnings(cor(gene_full, as.numeric(df.expr_ann[mk, samp_order_g]), method = "pearson")))
    corr_tbl <- data.frame(SYMBOL = panel, Corr = round(cors, 3)) %>% dplyr::arrange(dplyr::desc(Corr))
    names(corr_tbl)[2] <- paste0(gene, "_corr")
    res_xlsx(corr_tbl, sprintf("%s_correlation_panel", gene))
    cat(sprintf("%s - focus markers:\n", gene)); print(corr_tbl %>% dplyr::filter(SYMBOL %in% FOCUS3))
  }
  
  ## -- Violin plots (FOCUS3 only, saved as PNG) --
  save_violin_png(FOCUS3, sprintf("11_%s_Focus3_violin", gene), sprintf("%s / CSF3R-CD44-PROM1", gene),
                  samp_order_g, grp_levels, clin_ord_g, pal_g)
  
  for (mk in intersect(FOCUS3, rownames(df.expr_ann))) {
    p <- make_single_violin(mk, "Group", grp_levels, df.expr_ann, samp_order_g, clin_ord_g, pal_g)
    if (!is.null(p)) { print(p); fig_pdf(p, sprintf("14_preview_%s_%s", gene, mk), 5, 6) }
  }
  
  invisible(list(gene = gene, clin = clin_g, sample_order = samp_order_g))
}

## B6 - Run the pipeline separately for each PEA3 gene -----------------------
pea3_results <- lapply(PEA3_GENES, run_pea3_gene_analysis)
names(pea3_results) <- PEA3_GENES

cat("\n================================================================\n")
cat("COMPLETED (PEA3 = ETV1, ETV4, ETV5 - EACH ANALYZED SEPARATELY).\n")
cat("Figures :", FIG_DIR, "\nResults :", RES_DIR, "\n")
cat("================================================================\n")