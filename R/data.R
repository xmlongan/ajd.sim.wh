#' Heston SV Model Central Moment Formulas
#'
#' The first eight conditional central moment formulas for the Heston SV model.
#'
#' @format ## `fmu.hest`
#' A list with 8 elements:
#' \describe{
#'   \item{`fmu1`}{The first conditional central moment, equals 0}
#'   \item{`fmu2`}{The second conditional central moment}
#'   \item{...}{The ... conditional central moment}
#'   \item{`fmu8`}{The eighth conditional central moment}
#' }
#' Each of the formula consists of a list of 9 vectors:
#' \describe{
#'   \item{`e^{-kt}`}{Power of this term}
#'   \item{...}{Power of this term}
#'   \item{`rho`}{Power of this term}
#'   \item{`num`}{Numerator for this row coefficient}
#'   \item{`den`}{Denominator for this row coefficient}
#' }
"fmu.hest"

#' SVCJ Model Central Moment Formulas
#'
#' The first eight conditional central moment formulas for the SVCJ model.
#'
#' @format ## `fmu.svcj`
#' A list with 8 elements:
#' \describe{
#'   \item{`fmu1`}{The first conditional central moment, equals 0}
#'   \item{`fmu2`}{The second conditional central moment}
#'   \item{...}{The ... conditional central moment}
#'   \item{`fmu8`}{The eighth conditional central moment}
#' }
#' Each of the formula consists of a list of 15 vectors:
#' \describe{
#'   \item{`e^{kt}`}{Power of this term}
#'   \item{...}{Power of this term}
#'   \item{`sigma_s`}{Power of this term}
#'   \item{`num`}{Numerator for this row coefficient}
#'   \item{`den`}{Denominator for this row coefficient}
#' }
"fmu.svcj"

#' SVCJ Model Unconditional Central Moment Formulas
#'
#' The first eight unconditional central moment formulas for the SVCJ model.
#'
#' @format ## `fmu.svcj.unc`
#' A list with 8 elements:
#' \describe{
#'   \item{`fmu1`}{The first central moment, equals 0}
#'   \item{`fmu2`}{The second central moment}
#'   \item{...}{The ... central moment}
#'   \item{`fmu8`}{The eighth central moment}
#' }
#' Each of the formula consists of a list of 14 vectors:
#' \describe{
#'   \item{`e^{kt}`}{Power of this term}
#'   \item{...}{Power of this term}
#'   \item{`sigma_s`}{Power of this term}
#'   \item{`num`}{Numerator for this row coefficient}
#'   \item{`den`}{Denominator for this row coefficient}
#' }
"fmu.svcj.unc"
