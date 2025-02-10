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
#'   \item{`fmu8`}{The eight conditional central moment}
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
#'   \item{`fmu8`}{The eight conditional central moment}
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

#' #' SVVJ Model Central Moment Formulas
#' #'
#' #' The first eight conditional central moment formulas for the SVVJ model
#' #'
#' #' @format ## `fmu.svvj`
#' #' A list with 7 elements:
#' #' \describe{
#' #'   \item{`fmu2`}{A list with 2 elements}
#' #'   \item{...}{...}
#' #'   \item{`fmu8`}{A list with 5 elements}
#' #' }
#' #' Each of the formula consists of a few or several list, taking `fmu8` as
#' #' an example:
#' #' \describe{
#' #'   \item{`fmu8[[1]]`}{A list without summation part in each row}
#' #'   \item{`fmu8[[2]]`}{A list with 1 layer summation in each row}
#' #'   \item{...}{...}
#' #'   \item{`fmu8[[5]]`}{A list with 4 layer summations in each row}
#' #' }
#' #' where `fmu8[[1]]` is a list with two vectors:
#' #' \describe{
#' #'   \item{`base`}{row indexes of `svvj.base`}
#' #'   \item{`coef`}{row indexes of `svvj.coef`}
#' #' }
#' #' , `fmu8[[2]]` is a list with three vectors:
#' #' \describe{
#' #'   \item{`base`}{row indexes of `svvj.base`}
#' #'   \item{`lopq`}{row indexes of `svvj.jump$l1`}
#' #'   \item{`coef`}{row indexes of `svvj.coef`}
#' #' }
#' #' , `fmu8[[3]]` is a list with three vectors:
#' #' \describe{
#' #'   \item{`base`}{row indexes of `svvj.base`}
#' #'   \item{`lopq`}{row indexes of `svvj.jump$l2`}
#' #'   \item{`coef`}{row indexes of `svvj.coef`}
#' #' }
#' #' , `fmu8[[4]]` is a list with three vectors:
#' #' \describe{
#' #'   \item{`base`}{row indexes of `svvj.base`}
#' #'   \item{`lopq`}{row indexes of `svvj.jump$l3`}
#' #'   \item{`coef`}{row indexes of `svvj.coef`}
#' #' }
#' #' , and `fmu8[[5]]` is a list with three vectors:
#' #' \describe{
#' #'   \item{`base`}{row indexes of `svvj.base`}
#' #'   \item{`lopq`}{row indexes of `svvj.jump$l4`}
#' #'   \item{`coef`}{row indexes of `svvj.coef`}
#' #' }
#' #' Note that all (and the followed) the row indexes starting from 0,
#' #' which is different the starting index 1 adopted in R.
#' "fmu.svvj"
#'
#'
#' #' SVVJ Model Central Moment Formulas Component 1
#' #'
#' #' The formula component 1.
#' #'
#' #' @format ## `svvj.base`
#' #' A list with 9 vectors:
#' #' \describe{
#' #'   \item{`e^{-kt}`}{Power of this term}
#' #'   \item{...}{Power of this term}
#' #'   \item{`sqrt(1-rho^2)`}{Power of this term}
#' #' }
#' "svvj.base"
#'
#'
#' #' SVVJ Model Central Moment Formulas Component 2
#' #'
#' #' The formula component 2.
#' #'
#' #' @format ## `svvj.jump`
#' #' A list with 4 sub-lists:
#' #' \describe{
#' #'   \item{`l1`}{one levels of summation, explained later in Component 2-1}
#' #'   \item{`l2`}{two level of summation, explained later in Component 2-2}
#' #'   \item{`l3`}{three level of summation, explained later in Component 2-3}
#' #'   \item{`l4`}{four level of summation, explained later in Component 2-4}
#' #' }
#' #'
#' #' Please refer to
#' #' \url{http://yyschools.com/ajdmom/generated/ajdmom.ito_cond_mom.html}
#' #' for the details about the formula.
#' #'
#' #' The formula component 2-1.
#' #'
#' #' Sub-list `l1` a list with 2 vectors:
#' #' \describe{
#' #'   \item{`l1`}{power corresponding to \eqn{l_1}, see
#' #'   \url{http://yyschools.com/ajdmom/generated/ajdmom.ito_cond_mom.html}}
#' #'   \item{`o1`}{power corresponding to \eqn{o_1}}
#' #' }
#' #'
#' #' The formula component 2-2.
#' #'
#' #' Sub-list `l2` a list with 6 vectors:
#' #' \describe{
#' #'   \item{`l1`}{power corresponding to \eqn{l_1}}
#' #'   \item{`l2`}{power corresponding to \eqn{l_2}}
#' #'   \item{`o1`}{power corresponding to \eqn{o_1}}
#' #'   \item{`o2`}{power corresponding to \eqn{o_2}}
#' #'   \item{`p2`}{power corresponding to \eqn{p_2}}
#' #'   \item{`q2`}{power corresponding to \eqn{q_2}}
#' #' }
#' #'
#' #' The formula component 2-3.
#' #'
#' #' Sub-list `l3` a list with 10 vectors:
#' #' \describe{
#' #'   \item{`l1`}{power corresponding to \eqn{l_1}}
#' #'   \item{`l2`}{power corresponding to \eqn{l_2}}
#' #'   \item{`l3`}{power corresponding to \eqn{l_3}}
#' #'   \item{`o1`}{power corresponding to \eqn{o_1}}
#' #'   \item{`o2`}{power corresponding to \eqn{o_2}}
#' #'   \item{`o3`}{power corresponding to \eqn{o_3}}
#' #'   \item{`p2`}{power corresponding to \eqn{p_2}}
#' #'   \item{`p3`}{power corresponding to \eqn{p_3}}
#' #'   \item{`q2`}{power corresponding to \eqn{q_2}}
#' #'   \item{`q3`}{power corresponding to \eqn{q_3}}
#' #' }
#' #'
#' #' The formula component 2-4.
#' #'
#' #' Sub-list `l4` a list with 14 vectors:
#' #' \describe{
#' #'   \item{`l1`}{power corresponding to \eqn{l_1}}
#' #'   \item{`l2`}{power corresponding to \eqn{l_2}}
#' #'   \item{`l3`}{power corresponding to \eqn{l_3}}
#' #'   \item{`l4`}{power corresponding to \eqn{l_4}}
#' #'   \item{`o1`}{power corresponding to \eqn{o_1}}
#' #'   \item{`o2`}{power corresponding to \eqn{o_2}}
#' #'   \item{`o3`}{power corresponding to \eqn{o_3}}
#' #'   \item{`o4`}{power corresponding to \eqn{o_4}}
#' #'   \item{`p2`}{power corresponding to \eqn{p_2}}
#' #'   \item{`p3`}{power corresponding to \eqn{p_3}}
#' #'   \item{`p4`}{power corresponding to \eqn{p_4}}
#' #'   \item{`q2`}{power corresponding to \eqn{q_2}}
#' #'   \item{`q3`}{power corresponding to \eqn{q_3}}
#' #'   \item{`q4`}{power corresponding to \eqn{q_4}}
#' #' }
#' "svvj.jump"
#'
#'
#' #' SVVJ Model Central Moment Formulas Component 3
#' #'
#' #' The formula component 3.
#' #'
#' #' @format ## `svvj.coef`
#' #' A list with 2 vectors:
#' #' \describe{
#' #'   \item{`num`}{numerator}
#' #'   \item{`den`}{denominator}
#' #' }
#' "svvj.coef"
#'
#'
#' #' SVVJ Model Central Moment Formulas for the One-Jump Special Case
#' #'
#' #' The first eight conditional central moment formulas for the SVVJ model
#' #' when there are only one jump occurred. The 2|3|4 layers of summations
#' #' are shrank into only 1 layer of summation.
#' #'
#' #' @format ## `fmu.svvj.special`
#' #' Similar to `fmu.svvj`, but only contains the necessary data, i.e.,
#' #' the shrank results of the 2|3|4 layers of summations
#' "fmu.svvj.special"
#'
#'
#' #' SVVJ Model Central Moment Formulas One-Jump Special Component 2
#' #'
#' #' The formula component 2 when there is only one jump occurred.
#' #'
#' #' @format ## `svvj.jump.special`
#' #' A list with 4 sub-lists, similar to `svvj.jump`:
#' #' \describe{
#' #'   \item{`l1`}{empty (0), the same as that in `svvj.jump`}
#' #'   \item{`l2`}{shrank two levels of summation}
#' #'   \item{`l3`}{shrank three levels of summation}
#' #'   \item{`l4`}{shrank four levels of summation}
#' #' }
#' "svvj.jump.special"
#'
#'
#' #' SVVJ Model Central Moment Formulas One-Jump Special Component 3
#' #'
#' #' The formula component 3.
#' #'
#' #' @format ## `svvj.coef.special`
#' #' A list with 2 vectors:
#' #' \describe{
#' #'   \item{`num`}{numerator}
#' #'   \item{`den`}{denominator}
#' #' }
#' "svvj.coef.special"
