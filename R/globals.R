# Column names used in NSE (dplyr verbs, ggplot aes) are invisible to
# R CMD check, which reports them as undefined globals. Declaring them here
# keeps that NOTE quiet so a *real* unresolved symbol stands out -- which
# matters: the same NOTE was burying ~40 genuinely broken calls, including
# the ones that made plot_tale_composition() fail outright.
#
# Regenerate after adding NSE columns; do not hand-edit casually.

utils::globalVariables(c(
  "# Seq-ID", ".", "AnnoTALELength", "Approx. p-value", "Dissim",
  "EBEstrand", "RANK", "RVD", "RVDs", "RepU1", "RepU2", "SCORE", "SEQ_ID",
  "Score", "Seq", "Sequence", "Sim", "Strand", "TAL1", "TAL2", "TALBS_end",
  "TALBS_sequence", "TALBS_start", "TALE", "TAL_ID", "TAL_SEQ", "aaSeq",
  "aaSeqLength", "arlemScore", "arrayID", "dnaSeq", "domCode",
  "domainType", "ebeSeq", "group", "isNaAaSeq", "label",
  "matchConsensusRepeat", "matchConsensusRvd", "maxLength", "n", "name",
  "normArlemScore", "pattern", "pident", "position",
  "position in
  uncorrected sequences", "positionInArray",
  "positionInCrd", "qcov", "query", "queryHits", "query_name",
  "repeatClusterId", "repeatID", "repeatSimVsRef", "rvd",
  "rvd2ntMatchScore", "rvdFileLength", "rvdfac", "rvds", "rvdseq", "score",
  "seqnames", "sourceDirectory", "strain", "strand", "string", "subj",
  "subjSeqId", "subjectHits", "subtree", "taleId", "target", "target_name",
  "tcov", "value", "variable", "xPos", "yPos", "yPosOnSeq"
))
