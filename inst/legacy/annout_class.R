
# ---------------------------------------------------------------------------
# annout -- retired 2026-09.
#
# An exported S4 class: an AAStringSet with one extra data.frame slot. Its
# only role was to let sapply() inside tell_tales() return an RVD sequence and
# a domains report together as a single object per TALE. The loop is an
# lapply() returning list(rvds =, domains =) now, which does the same thing
# without an implementation detail appearing in the package's API.
#
# Nothing else ever constructed one, and no method was ever defined on it.
# ---------------------------------------------------------------------------

#' annotale output 
#' @exportClass annout
#' @import Biostrings
#' @importFrom methods setClass
annout <- setClass(
  # Set the name for the class
  Class = "annout",
  
  # Define the slots
  slots = c(
    domainsReport = "data.frame"
  ),
  
  contains = "AAStringSet",
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity = function(object) {
    val <- is.data.frame(object@domainsReport)
    return(val)
  }
)
