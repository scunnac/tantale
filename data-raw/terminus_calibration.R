# Calibrate what counts as a canonical TALE terminus (ledger §59).
#
# NTERM/CTERM mark a terminus whose protein match to the TALE N- or
# C-terminal profile is full-length or close to it. "Close" is learned here
# from the termini of curated TALEs: the TALEs of tale_annotations (ten
# published X. oryzae genomes), minus the truncTALEs and the TALEs with an
# unusual_feature, whose two termini are left out (maintainer, Q174).
#
# Steps:
#   1. download the ten assemblies from NCBI (GenBank accessions in
#      tale_annotations$genome_id)
#   2. run tell_tales() on each
#   3. match its arrays to the curated TALEs by RVD string, within a strain
#   4. tabulate the terminus match features (array_report.tsv) by label:
#      canonical, truncTALE, unusual (an unusual_feature), unannotated (an
#      array with no curated TALE of that RVD string)
#
# Large files go to CALIB_DIR (default: ../tantale_calibration, outside the
# repository); rerunning skips what is already there.
#
# Run with:  Rscript data-raw/terminus_calibration.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

CALIB_DIR <- Sys.getenv("TANTALE_CALIB_DIR", normalizePath(file.path("..", "tantale_calibration"),
                                                           mustWork = FALSE))
genome_dir <- file.path(CALIB_DIR, "genomes")
dir.create(genome_dir, recursive = TRUE, showWarnings = FALSE)

genomes <- unique(tale_annotations[c("strain", "genome_id")])

#### 1. Assemblies ####
genome_path <- function(strain) file.path(genome_dir, paste0(strain, ".fa"))
for (i in seq_len(nrow(genomes))) {
  out <- genome_path(genomes$strain[i])
  if (file.exists(out)) next
  zip <- tempfile(fileext = ".zip")
  url <- paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/",
                genomes$genome_id[i], "/download?include_annotation_type=GENOME_FASTA")
  utils::download.file(url, zip, mode = "wb", quiet = TRUE)
  fna <- grep("_genomic\\.fna$", utils::unzip(zip, list = TRUE)$Name, value = TRUE)
  stopifnot(length(fna) == 1L)
  utils::unzip(zip, files = fna, exdir = tempdir(), junkpaths = TRUE)
  file.rename(file.path(tempdir(), basename(fna)), out)
  unlink(zip)
}
genomes$path <- genome_path(genomes$strain)
print(transform(genomes, size_mb = round(file.size(path) / 1e6, 1))[c("strain", "genome_id", "size_mb")])

#### 2. tell_tales() ####
run_dir <- file.path(CALIB_DIR, "telltales")
for (strain in genomes$strain) {
  out <- file.path(run_dir, strain)
  if (file.exists(file.path(out, "array_report.tsv"))) next
  unlink(out, recursive = TRUE)
  message("tell_tales() on ", strain)
  suppressWarnings(suppressMessages(tell_tales(genome_path(strain), output_dir = out)))
}

#### 3. Arrays matched to the curated TALEs ####
bare_rvds <- function(x) {
  x <- toupper(x)
  x <- gsub("^(NTERM|XXXXX)-|-(CTERM|XXXXX)$", "", x)
  x
}
reports <- lapply(stats::setNames(nm = genomes$strain), function(strain) {
  r <- readr::read_tsv(file.path(run_dir, strain, "array_report.tsv"), show_col_types = FALSE)
  r$strain <- strain
  r
})
reports <- dplyr::bind_rows(reports)
reports$rvds <- bare_rvds(reports$rvd_string)

curated <- tale_annotations
curated$rvds <- toupper(curated$rvd_seq)
curated$label <- dplyr::case_when(curated$truncTALE ~ "truncTALE",
                                  !is.na(curated$unusual_feature) ~ "unusual",
                                  TRUE ~ "canonical")
# one label per strain and RVD string; identical TALEs with different
# labels would make the match ambiguous
labels <- curated %>%
  dplyr::group_by(strain, rvds) %>%
  dplyr::summarise(label = if (dplyr::n_distinct(label) == 1L) label[1] else "ambiguous",
                   tal_names = paste(stats::na.omit(tal_name), collapse = "/"),
                   .groups = "drop")
reports <- dplyr::left_join(reports, labels, by = c("strain", "rvds"))
reports$label[is.na(reports$label)] <- "unannotated"
reports$label[reports$rvds == ""] <- "no RVDs"

cat("\nArrays by label:\n"); print(table(reports$label))
missed <- dplyr::anti_join(labels, reports, by = c("strain", "rvds"))
cat("\nCurated TALEs with no array of the same RVD string:", nrow(missed), "\n")
print(as.data.frame(missed))

#### 4. Terminus features ####
features <- dplyr::bind_rows(lapply(c(nterm = "N-terminus", cterm = "C-terminus"), function(end) {
  pre <- if (end == "N-terminus") "nterm_" else "cterm_"
  cols <- c("dna_hit", "dna_score", "dna_cover", "dna_pieces", "aa_evalue", "aa_score",
            "aa_profile_gap", "aa_far_gap", "aa_cover", "aa_domains", "aa_hit", "aa_length")
  x <- reports[c("strain", "array_id", "label", "tal_names", paste0(pre, cols))]
  names(x) <- c("strain", "array_id", "label", "tal_names", cols)
  x$terminus <- end
  x
}))
features <- features[c("strain", "array_id", "tal_names", "label", "terminus",
                       setdiff(names(features), c("strain", "array_id", "tal_names", "label", "terminus")))]
readr::write_tsv(features, file.path(CALIB_DIR, "terminus_features.tsv"))

summarise_features <- function(x) {
  x %>%
    dplyr::group_by(terminus, label) %>%
    dplyr::summarise(n = dplyr::n(),
                     no_segment = sum(is.na(aa_length)),
                     aa_cover_min = min(aa_cover, na.rm = TRUE),
                     aa_far_gap_max = max(aa_far_gap, na.rm = TRUE),
                     aa_profile_gap_max = max(aa_profile_gap, na.rm = TRUE),
                     aa_domains_max = max(aa_domains, na.rm = TRUE),
                     aa_length_min = min(aa_length, na.rm = TRUE),
                     aa_length_max = max(aa_length, na.rm = TRUE),
                     dna_cover_min = min(dna_cover, na.rm = TRUE),
                     .groups = "drop")
}
cat("\nTerminus features by label:\n")
print(as.data.frame(suppressWarnings(summarise_features(features))))
