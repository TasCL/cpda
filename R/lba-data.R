#' Probability Density Approximation Data Sets
#'
#' An example data set to fit piecewise LBA model, converted from Holmes's 
#' (2015) MATLAB formatted Example 4 Files (Data1.mat).
#'
#' @docType data
#' @keywords dataset
#' @name d
#'
#' @source \href{http://dx.doi.org/10.1016/j.jmp.2015.08.006}{PDA paper}
#' @examples
#' \dontrun{
#' data(lba)
#' head(d)
#' }
#' ##   Response ResponseTime Block
#' ## 1        0    0.5676343     1
#' ## 2        1    0.6183177     1
#' ## 3        0    0.8806298     1
NULL

#' Probability Density Approximation Data Sets
#'
#' An example data set to fit piecewise LBA model, converted from Holmes's
#' MATLAB formatted Example 4 Files (Data1.mat).
#'
#' @docType data
#' @keywords dataset
#' @name plba
#'
#' @source \href{http://dx.doi.org/10.1016/j.jmp.2015.08.006}{PDA paper}
#' @examples
#' data(lba)
#' str(plba)
#' ## List of 4
#' ## $ DT1 : num [1:695] 0.488 0.801 0.376 0.507 0.532 ...
#' ## $ DT2 : num [1:305] 0.538 0.77 0.568 0.271 0.881 ...
#' ## $ eDT1: num [1:7020] 0.475 0.346 0.42 0.401 0.368 ...
#' ## $ eDT2: num [1:2980] 0.703 0.693 0.704 0.462 0.468 ...
NULL
