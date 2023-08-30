.onAttach <- function(...){
  packageStartupMessage(logger::log_info('Loading the tantale package'))
  packageStartupMessage(logger::log_debug('Email sebastien.cunnac@ird.fr for comments'))
  logger::log_threshold('INFO', namespace = 'tantale')
}

.onLoad <- function(...){
  logger::log_layout(logger::layout_glue_colors,namespace ='tantale')
}

g <- glue::glue
m <- dplyr::mutate
s <- dplyr::select
gb <- dplyr::group_by