#' Data from the JakStat signaling
#' 
#' This is the data from the Jak2/Stat5 signaling pathway. The data is taken
#' from the works of Timmer et al. (DOI 10.1038/msb.2011.50).
#' 
#' @details  A data frame with 16 rows and 5 columns
#' \describe{
#'    \item{t}{t in minutes}
#'    \item{y1}{phosphorylated STAT5 in cytoplasma}
#'    \item{y1sd}{standard deviation of the measurement of y1 at time t}
#'    \item{y2}{total STAT5 in cytoplama}
#'    \item{y2sd}{standard deviation of the measurement of y2 at time t}
#' }
#' @source \url{http://msb.embopress.org/content/7/1/516}
"jakstatMeasurement"


#' Input of the Jak2/Stat5 signaling pathway
#' 
#' A dataframe that containing the inputs of the Jak2/Stat5 pathway that is
#' used to stimulate the system. The data is taken
#' from the works of Timmer et al. (DOI 10.1038/msb.2011.50).
#' 
#' @details  A dataframe with 16 rows and 2 columns
#' \describe{
#'    \item{t}{t in minutes}
#'    \item{u}{activation of Epo-receptor}
#' }
#' @source \url{http://msb.embopress.org/content/7/1/516}
"jakstatInput"
