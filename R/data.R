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

#' UVB signal pathway
#' 
#' A data frame containing simulated values of the UVB Signaling pathway. The error
#' of the system is synthetic and is added to the states x3 and x11. The model
#' is taken from the works of Ouyang et al. \url{https://doi.org/10.1073/pnas.1412050111}
#' 
#' @details A data frame with 8 rows and 11 columns
#' \describe{
#'     \item{t}{time in fractions of an hour}
#'     \item{y1}{total amounts of UVR8 monomers}
#'     \item{y2}{total amounts of COP1 monomers}
#'     \item{y3}{total amounts of UVR8 dimers}
#'     \item{y4}{concentration of elongated hypocotyl 5 (HY5) protein}
#'     \item{y5}{concentration measured of UVR8 monomers}
#'     \item{y1sdt}{standard deviation of the first measurement}
#'     \item{y2sdt}{standard deviation of the second measurement}
#'     \item{y3sdt}{standard deviation of the third measurement}
#'     \item{y4sdt}{standard deviation of the fourth measurement}
#'     \item{y5sdt}{standard deviation of the fifth measurement}
#' }
#' @source \url{https://doi.org/10.1073/pnas.1412050111}
"uvbData"
