
# Functions that should be decommissioned

- tale_parts_to_rvd
- repeat_to_rvd_map_distalr: output can be trivially obtained from a tales object
- repeat_to_rvd_map: this was needed when run_distal was in effect to map the annoTale RDV sequences to the domain codes. But now that there is the tales_compare function, I do not see a convincing use case.
- .rvd_to_repeat_align (never called)

# Still needed?
.pairwise_distances_rename_legacy
.tales_rename_legacy


# Renaming

tales_width() should be renamed to tales_msa_width() or should be converted to a method for the width() generic.
tales_from_telltale should be renamed to tales_from_telltales because it takes results from tell_tales
repeat_to_rvd_align -> .repeat_to_rvd_align because it is internal.

# Functions that I inspected to determine what action needs to be taken

.domain_to_cluster_align -> keep
.domain_to_sim_align -> keep
.rvd_to_match_align -> keep
repeat_to_rvd_align (never called) -> keep because it builds fixture data in test_plot_tales_msa.R.
.rvds_from_annotale_file -> keep just in case. It be used on the run_annotale_predict() output, but we may as well want to write a tales_from_annotale() function, similar to tales_from_telltale(). Actually this MUST be done and advertised in website.

# Redundant functions?
 - .tidy_biostrings_msa
 Not redundant but problem with @param msa multiple sequence alignment file -> reading from file is not implemented

# Code refactoring

in plot.tales_msa the content of both 'if (!is.null(rvd_align)) {' blocks should be merged to improve legibility.



# What kind of shape conversion do we need? What do we have?

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






