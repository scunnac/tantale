# Build inst/extdata/annotaleExampleOutput, read by tales_from_annotale()'s
# example and tests: AnnoTALE predict + analyze on the shipped BAI3 sample
# (4 TALEs), keeping analyze's three files and predict's GFF3 (for the
# contig names). predict's GenBank and sequence files are left out for size.
#
# Rerun after updating the AnnoTALE jar.
#
# Run with:  Rscript data-raw/make_annotale_example_output.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

subject <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                       package = "tantale", mustWork = TRUE)
run <- file.path(tempdir(), "annotale_example")
unlink(run, recursive = TRUE)
run_annotale_predict(subject, output_dir = run)

to <- "inst/extdata/annotaleExampleOutput"
unlink(to, recursive = TRUE)
dir.create(file.path(to, "Analyze"), recursive = TRUE)
dir.create(file.path(to, "Predict"))
file.copy(file.path(run, "Analyze", c("TALE_Protein_parts.fasta", "TALE_DNA_parts.fasta",
                                      "TALE_RVDs.fasta")),
          file.path(to, "Analyze"))
# renamed: AnnoTALE's name has parentheses, which R CMD check calls
# non-portable; tales_from_annotale() looks for "GFF__*.gff3"
file.copy(list.files(file.path(run, "Predict"), "^GFF__.*\\.gff3$", full.names = TRUE),
          file.path(to, "Predict", "GFF__TALE_predictions_bai3_sample.gff3"))
