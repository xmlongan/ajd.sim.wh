#' Write Parallel Result to csv file
#'
#' Write the result from foreach loop to a csv file.
#'
#' @param err_dur list of two-element vectors.
#' @param fname file name fo the csv file.
#' @param append append or rewrite, default to `FALSE`.
#'
#' @export
#'
#' @examples
#' # err_dur = list(c(0.14, 10), (0.15, 11))
#' # write_to_csv(err_dur, "error-computing_time.csv")
write_to_csv <- function(err_dur, fname, append=FALSE) {
  N = length(err_dur); err = rep(0,N); dur = rep(0,N)
  for (i in 1:N) {
    err[i] = err_dur[[i]][1]
    dur[i] = err_dur[[i]][2]
  }
  df = data.frame(error=err, computing_time=dur)
  write.csv(df, file = fname, append)
}
