# Column names used in NSE (dplyr verbs, ggplot aes) are invisible to
# R CMD check, which reports them as undefined globals. Declaring them here
# keeps that NOTE quiet so a *real* unresolved symbol stands out -- which
# matters: the same NOTE was burying ~40 genuinely broken calls, including
# the ones that made plot_tales_composition() fail outright.
#
# Regenerate after adding NSE columns; do not hand-edit casually.

utils::globalVariables(c(
  "# Seq-ID", ".", "aa_seq", "array_id", "domain_type", "id1", "id2", "dissim", "sim",
  "arlem_score", "max_length", "norm_arlem_score", "position_in_array", "position_in_crd", "dna_seq", "dom_code", "source_directory", "rvdSimVsRef", ".x", "rvd1", "position in uncorrected sequences", "AnnoTALELength", "Approx. p-value", "EBEstrand", "RANK", "RVD", "RVDs", "SCORE", "SEQ_ID",
  "Score", "Seq", "Sequence", "Strand", "TALBS_end",
  "TALBS_sequence", "TALBS_start", "TALE", "TAL_ID", "TAL_SEQ",   "aa_length", "ebeSeq", "group", "label",
  "matchConsensusRepeat", "matchConsensusRvd", "n", "name",
  "pident", "position",
  "position in
  uncorrected sequences",   "qcov", "query", "queryHits", "query_name",
  "repeatClusterId", "repeatID", "repeatSimVsRef", "rvd",
  "rvd2ntMatchScore", "rvdFileLength", "rvdfac", "rvds", "rvdseq", "score",
  "seqnames", "strain", "strand", "string", "subjSeqId", "subjectHits", "subtree", "taleId", "target", "target_name",
  "tcov", "value", "variable", "xPos", "yPos", "yPosOnSeq"
))
