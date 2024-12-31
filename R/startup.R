
greet_startup_cli <- function() {
  cli_inform(c("v" = "Attaching the tantale package",
               ">" = "Email sebastien.cunnac@ird.fr for comments"),
             class = "packageStartupMessage")
}


.onAttach <- function(...){
  packageStartupMessage(greet_startup_cli(), appendLF = FALSE)
}

