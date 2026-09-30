
# Functions that should be decommissioned

- tale_parts_to_rvd
- repeat_to_rvd_map_distalr: output can be trivially obtained from a tales object
- repeat_to_rvd_map: this was needed when run_distal was in effect to map the annoTale RDV sequences to the domain codes. But now that there is the tales_compare function, I do not see a convincing use case.
- .rvd_to_repeat_align (never called)

# Still needed?
.pairwise_distances_rename_legacy
.tales_rename_legacy


# Renaming

- tales_width() should be renamed to tales_msa_width() or should be converted to a method for the width() generic.

- tales_from_telltale should be renamed to tales_from_telltales because it takes results from tell_tales

- repeat_to_rvd_align -> .repeat_to_rvd_align because it is internal.

# Functions that I inspected to determine what action needs to be taken

- .domain_to_cluster_align -> keep

- .domain_to_sim_align -> keep

- .rvd_to_match_align -> keep

- repeat_to_rvd_align (never called) -> keep because it builds fixture data in test_plot_tales_msa.R.

- .rvds_from_annotale_file -> keep just in case. It be used on the run_annotale_predict() output, but we may as well want to write a tales_from_annotale() function, similar to tales_from_telltale(). Actually this MUST be done and advertised in website.

# Redundant functions?
 - .tidy_biostrings_msa
 Not redundant but problem with @param msa multiple sequence alignment file -> reading from file is not implemented

# Code refactoring

in plot.tales_msa the content of both 'if (!is.null(rvd_align)) {' blocks should be merged to improve legibility.


# A proposal to bring package size in line with rOpenSci requirements.

Place jars and fixtures strictly necessary for the website in an archive.
Insert that archive in a 'assets' branch.
Update `tantale_setup()` to have it download the archive from git, unpack and move the jars to a more standard jar folder and the fixtures to inst/extdata.

In preparation, we need to write a short manuscript that would fit in the Journal of Open Software and deposit it in BiorXiv.

The package needs to have a CRAN or OSI accepted license. The R packages book includes a helpful section on licenses.

For testing your functions creating plots, we suggest using vdiffr, an extension of the testthat package that relies on testthat snapshot tests.

https://devguide.ropensci.org/pkg_ci.html#whichci -> this is going to be a major challenge for me.


# Misc.

plot.tales() should have a parameter specifying the variable used for the ggplot2 wrap.

note that if you name the main vignette of your package “pkg-name.Rmd”, it’ll be accessible from the navbar as a Get started link instead of via Articles > Vignette Title.




# The many issues of tell_tales and tales_from_telltales

The log of tell_tale is displayed weirdly in RStudio: as if the line breaks where converted to space.

**tell_tale should clearly indicate that its array report analyses hmmer hits**. In the doc and messages there should be no ambiguity as to whether we are refering to dna hits or AnnoTALE output:
"has_all_domains" columns should be renamed to something like "all_domains_in_dna_hits".. This is too long but we need to find something less missleading. 'Domain' evoke protein but the check is done on hmmer dna profiles. The same comment is valid for "n_domain_hits".

When running:
```{r discover_bai311_raw, message = FALSE}
invisible(tell_tales(subject_file = bai311_fa, output_dir = bai311_raw_dir,
                     cterm_min_score = 300))
```
These warnings should be accounted for:
```
Warning messages:
1: In GenomeInfoDb::renameSeqlevels(gr, value = seqlevels) :
  invalid seqlevels 'seq2' ignored
2: Some HMMER hits overlap, so the inferred RVD sequences may carry artefactual insertions.
ℹ Check these regions: "ROI_00001", "ROI_00002", "ROI_00003", "ROI_00005", "ROI_00006", "ROI_00007", "ROI_00008", and "ROI_00009" 
3: Annotale failed to parse TALE domains for ROI_00003. 
4: Annotale failed to parse TALE domains for ROI_00005. 
```

Reviewing the tale_mining article of website made me realise that we need to audit how termini codes are added to arrays in tell_tales and tales_from_telltale because I suspect **something is flawed or misleading** (see my pain explaining what happened to the weird ROIs in my review of the site):

There are inconsistencies in the way tale parts are reconciliated between tell_talse dna hits and AnnoTALE.
tell_tale can identify all domains at the dna level but because of one or several frame shifts, the longest orf occupies only a subregion of the array sequence. AnnoTALE takes the predicted longest orf on the predicted 'array' DNA span and identify tale 'domains' on both orf and protein product. As far as I understand it, there is no guarantee that the terminus domains are canonical (true? yes). Furthermore, there is no guarantee that both tell_tales and AnnoTALE sequences correctly aligns.
Then the question is what are the inconsistencies between these two visions?
A the core of my indecision is that I would like to check if a termini reported by AnnoTALE is truly supported by an underlying dna domain hit or not. The problem is that I do not know how to do that.
What I could do is take the AnnoTALE termini AA or dna and ask if that fits with the corresponding hmm. What I could also ask is if in this array, tell_tales detected a hmmmer for this termini irrespective of where it is.
This would prevent plainly erroneous cases where AnnoTALE N- or C-term sequences have nothing to do with a genuine TALE terminal domain but are considered as such because a corresponding nhmmer hit has been found earlier...
Where shall this be implemented? Probably in .tale_parts or in .tale_parts_from_file (may be more complicated)



- In .tale_parts rvds should be obtained from AnnoTALE "TALE_RVDs.fasta" files with .rvds_from_annotale_file and not from the tell_tales output (which is AnnoTALE rvds plus termini if nhmmer hits found on dna): .tale_parts parse RVDs from tell_tale "rvd_sequences.fas" which is obtained by using AnnoTALE and then adds terminus codes (if requested) using this code:
```r
  if (extremity_codes) {
    # This is necessary for other tantale utilities that can operate on 'full'
    # domains sequences, ie downstream of distal, for TALE domains sequences
    # alignments.
    anchors <- tales_anchor_codes()   # NTERM, CTERM, XXXXX
    for (s in names(rvds)) {
      present <- as.character(by_array[[s]]$query_name)
      rvds[s] <- paste(ifelse(hmm$nterm %in% present, anchors[1], anchors[3]),
                       rvds[s], sep = rvd_sep)
      rvds[s] <- paste(rvds[s],
                       ifelse(hmm$cterm %in% present, anchors[2], anchors[3]),
                       sep = rvd_sep)
    }
  }
```
Which translate in "If a hhmm dna domain is found for the corresponding terminus, append it, if not,
append "XXXXX"". This could be recycled to have the info about nhmmer hits in .tale_parts but other than that should be discounted everywhere.

- In addition regarding tale_parts_from_file, some critical points:
  1. I am not sure .tale_parts_from_file fires warnings if termini are not found by AnnoTALE.
  2. **.tale_parts_from_file adds NA rows in the table for missing "N-terminus" or "C-terminus".**
  3. .tale_parts_from_file always return arrays tibble with both terminus. If not found by AnnoTALE, string is NA

- Instead of rvd_string in "array_report.tsv" should report found_N-term_nhmmer_hit and found_C-term_nhmmer_hit 








# What kind of shape conversion do we need? What do we have?

 -> NOTHING TO DO HERE, just to be kept as a note.

  - vector of seq (rvd, dom_code) with separator -> list of vectors
      * .split_list

  - vector of seq with separator -> list of vectors
      * .tidy_biostrings_msa (sep = "")

  - algn matrix -> align matrix with other values
      * .domain_to_cluster_align
      * .domain_to_sim_align
      * .rvd_to_match_align
      * .rvd_to_repeat_align (never called)
      * repeat_to_rvd_align (never called)

  - align matrix -> to long (ie ~tales_msa)
      * .matrix_to_long()

  - tales_msa -> align matrix
      * as.matrix.tales_msa
  
  - tales: vector of tales 'values' (rvd, dom_code, aa, dna) with a specified separator ("-", "", ...) -> biostring set
      * tales_coded_strings
      * tales_rvd_strings
      * tales_get_dna_seq
      * tales_get_protein_seq






