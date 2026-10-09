# Passive voice in the documentation: proposed rewrites (Q224)

Prepared 2026-10-09 for the maintainer. Source: every user-facing help
page (`man/*.Rd` without `\keyword{internal}`), the articles, the Get
started page, `pkgdown/index.md` and `README.qmd`, searched for "is/are/
was/were/be/been/being + participle". That gives about 360 hits. Most
are fine and are not listed: the agent is unknown or irrelevant ("is
pinned", "are refused", "is required"), the participle is an adjective
("is truncated", "are related"), or the sentence describes a state
("is stored with the object"). Listed below are the cases where the
agent is known and an active sentence is shorter or clearer.

Mark each line **Y** (rewrite as proposed), **N** (keep) or give your own
wording. Locations are the roxygen source, `file:line` on 2026-10-09.

## Help pages

**N** 1. `R/target_predictions.R:20` and `:149` (`talvez()`, `preditale()`, `rvd_seqs`).
   Now: "Tale RVD sequences are supplied as either a fasta file ... or as
   a Biostrings XStringSet ..."
   Proposed: "The TALE RVD sequences: the path to a fasta file whose
   headers name the TALEs and whose sequences are RVDs separated by spaces
   or hyphens, or a Biostrings XStringSet in the same format."
**N** 2. `R/target_predictions.R:332` (`plot_target_preds()`).
   Now: "RVDs sequences predicted to target an EBE on the sense strand ...
   are plotted on top of the double stranded DNA sequence ... Those
   predicted to target an EBE on the opposite strand are displayed below."
   Proposed: "The plot draws the RVD sequences that target an EBE on the
   sense strand above the double-stranded DNA, each next to its EBE,
   which it highlights on that strand; those targeting the opposite strand
   go below."
3. **Y** `R/target_predictions.R:334` (same page).
   Now: "Individual RVDs are printed inside colored boxes. The color of the
   boxes indicates to which degree the RVD is predicted to have affinity
   ..."
   Proposed: "Each RVD sits in a box whose colour says how strongly the
   RVD is expected to bind the base at that position, compared with the
   other three bases."
**N** 4. `R/target_predictions.R:340` (`subj_file`).
   Now: "The fasta file of subject DNA sequences that was used to predict
   DNA binding elements."
   Proposed: "The fasta file of subject DNA sequences the predictions were
   made on." (still passive, but shorter) or "The fasta file of subject DNA
   sequences given to the predictor."
**Y** 5. `R/tales_plot.R:195` (`plot.tales_msa()`).
   Now: "Three things are decided independently, and it helps to read the
   figure that way: ..."
   Proposed: "Three settings act independently, and it helps to read the
   figure that way: ..."
**Y** 6. `R/tales_plot.R:238` (same page).
   Now: "`ref_pattern` is matched against the array names and must
   identify exactly one, otherwise the default is used with a warning"
   Proposed: "`ref_pattern` must match exactly one array name; otherwise
   the function warns and uses the default"
**Y** 7. `R/classification.R:38` and `:152` (`tales_group_hclust()`,
   `tales_group_kmedoids()`).
   Now: "The clustering is computed from `tale_distances`, but the result
   belongs on the `tales` object the distances were computed from, so that
   is what comes back."
   Proposed: "The function clusters `tale_distances` but returns the
   `tales` object those distances came from, with the groups added."
**Y** 8. `R/classification.R:358` (`talomes_heatmap()`).
   Now: "Within a group, variants are ranked by how many strains carry
   them, and the cell colour is that rank"
   Proposed: "Within a group, the cell colour gives each variant's rank by
   the number of strains that carry it"
**Y** 9. `R/functal.R:115-117` (`tales_compare_functal()`).
   Now: "Each array's repeats are turned into a position weight matrix
   (PWM) over the RVD-to-base specificity code, and PWMs are compared
   pairwise with compare_motifs()."
   Proposed: "The function turns each array's repeats into a position
   weight matrix (PWM) of the bases its RVDs prefer, and compares the PWMs
   pairwise with compare_motifs()."
**Y** 10. `R/functal.R:130-134` (same page).
    Now: "The PWMs themselves are built by tales_to_universalmotif() ...
    That conversion is exposed on its own so it can be used with any
    universalmotif function."
    Proposed: "tales_to_universalmotif() builds the PWMs (see its page for
    what drives them and how it treats an RVD missing from
    rvd_dna_specificity). It is exported so that its output can go to any
    universalmotif function."
**Y** 11. `R/telltale.R:1367-1371` (`tell_tales()`).
    Now: "Hits that are (nearly ...) adjacent are grouped in "taleArrays"
    which are considered as potential tal genes."
    Proposed: "It groups hits that are adjacent, or nearly so (see
    `min_gap`), into "taleArrays", each a potential *tal* gene."
**Y** 12. `R/telltale.R:1374-1376` (same page).
    Now: "If the correct_array parameter is turned off, the longest
    predicted open reading frame (+extend_len) for each talArray is fed to
    AnnoTALE ..."
    Proposed: "With `correct_array = FALSE`, tell_tales() gives AnnoTALE
    the longest open reading frame of each taleArray (extended by
    `extend_len`) ..." (also fixes "talArray"/"talearrays" spellings)
**Y** 13. `R/telltale.R:1381` (same page).
    Now: "If correct_array is turned on, these talearrays are passed to the
    CorrectFrameshifts function that attempts to 'correct' ..."
    Proposed: "With `correct_array = TRUE`, CorrectFrameshifts() first
    tries to correct frameshifts in the taleArrays ..."
**N** 14. `R/telltale.R:1392` (same page).
    Now: "This should be detected and reported in the tell_tales log."
    Proposed: "tales_from_telltales() warns about these arrays and leaves
    them out." (Checked: the tell_tales log does not report them;
    tales_from_telltales() does, naming the arrays. So the current
    sentence is also inaccurate.)
    Correction (Q241): this check was wrong. The log does list them,
    under "Noteworthy AnnoTale issues"; the sentence is accurate.
**N** 15. `R/telltale.R:1396` (same page).
    Now: "A tal gene that spans the junction of a circular molecule ... is
    cut in two, and is reported, if at all, ..."
    Proposed: "tell_tales() cuts a tal gene that spans the junction of a
    circular molecule ... in two, and reports it, if at all, ..."
**Y** 16. `R/tantale_data.R:143` (`tantale_genome()`).
    Now: "Four Xanthomonas oryzae pv. oryzae genome assemblies are used
    throughout the articles and examples."
    Proposed: "The articles and examples use four Xanthomonas oryzae pv.
    oryzae genome assemblies."
**N** 17. `R/tantale_data.R:152` (same page).
    Now: "It was produced by the authors of tantale, ..."
    Proposed: "The authors of tantale produced it; it is available only
    from ..."
**Y** 18. `R/tantale_setup.R:142` (`tantale_setup()`).
    Now: "Two archives, attached to releases of tantale's GitHub
    repository, are checked against a sha256 recorded in the package and
    unpacked into ..."
    Proposed: "The function checks two archives, attached to releases of
    tantale's GitHub repository, against a sha256 recorded in the package
    and unpacks them into ..."
**Y** 19. `R/tales_class.R:194` (`tales()`, `sanitize`).
    Now: "FALSE (default): none. The anomalies are warned about, so odd
    predictions can still be loaded and inspected."
    Proposed: "FALSE (default): none. tales() warns about the anomalies,
    and the odd predictions stay in the object for inspection."
**Y** 20. `R/annotale.R:21-23` (`run_annotale_predict()`).
    Now: "The whole AnnoTALE workflow can be completed by a subsequent call
    to the run_annotale_build function."
    Proposed: "run_annotale_build() completes the AnnoTALE workflow."
**Y** 21. `R/annotale.R:27` (same page).
    Now: "Path to a fasta file containing DNA sequences (e.g. a genome
    assembly) to be analyzed for TALE content."
    Proposed: "Path to a fasta file of the DNA sequences (e.g. a genome
    assembly) to search for TALEs."
**Y** 22. `R/pairwise_distances_class.R:75-79` (`pairwise_distances()`).
    Now: "A similarity is accepted and converted: a table carrying sim ...
    is folded into dissim, and those restatements are then dropped ..."
    Proposed: "The function also accepts a similarity: it converts a table
    carrying sim ... into dissim and drops the similarity columns, so that
    only one copy of the quantity remains. It keeps any further column
    (...)."

## Articles and home page

**Y** 23. `tale_mining.qmd:66`. "Each terminal region is searched with the
    profile HMM of the TALE N- or C-terminus" → "`tell_tales()` searches
    each terminal region with the profile HMM of the TALE N- or
    C-terminus".
**Y** 24. `tale_mining.qmd:97-99`. "Hits are merged, grouped into candidate
    arrays by proximity, and each array's longest ORF is handed to
    AnnoTALE ... Terminus codes are added at either end of the RVD string"
    → "`tell_tales()` merges the hits, groups them into candidate arrays by
    proximity and hands each array's longest ORF to AnnoTALE ... It then
    adds terminus codes at either end of the RVD string".
**Y** 25. `trunctale_correction.qmd:262`. "`max_comparisons = 20` is used above
    to keep this article's build time ..." → "This article uses
    `max_comparisons = 20` to keep its build time ...".
**N** 26. `tale_target_prediction.qmd:36`. "Both are wrapped as `talvez()` and
    `preditale()`, which take ..." → "tantale wraps them as `talvez()` and
    `preditale()`, which take ...".
**Y** 27. `tale_target_prediction.qmd:151-152`. "Each RVD is printed in a box
    over the base it is predicted to contact; predictions on the sense
    strand are drawn above the sequence ..." → "The plot puts each RVD in a
    box over the base it should contact, and draws the predictions on the
    sense strand above the sequence ...".
**Y** 28. `tale_classification.qmd:60`. "Discovery runs independently per genome,
    so it is written as one function ..." → "Discovery runs independently
    per genome, so the article writes it as one function ...".
**Y** 29. `tale_classification.qmd:98`. "`array_id` is prefixed by strain before
    the three objects are combined" → "The code prefixes `array_id` with
    the strain before combining the three objects".
**Y** 30. `tales_class.qmd:83`. "Two columns are required outright, because ..."
    → "A `tales` requires two columns outright, because ...".
**Y** 31. `tales_class.qmd:101-102`. "Any further column ... is preserved
    untouched and never warned about" → "`tales()` keeps any further column
    ... untouched and never warns about it".
**Y** 32. `README.qmd:42`. "Here is a snapshot of the topics that are covered:"
    → "The package covers:".
