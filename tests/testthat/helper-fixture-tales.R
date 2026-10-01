# The shared fixture holds one genuine truncTALE (PXO86_ROI_00019, its
# C-terminus coded XXXXX), which tales() reports as an anomaly. Tests that
# merely load the fixture do not need that warning.
tales_quietly <- function(parts, ...) {
  suppressWarnings(tales(parts, ...), classes = "tantale_warning_tales_anomalous")
}
