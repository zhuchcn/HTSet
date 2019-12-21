#' @rdname Extract
#' @aliases $
#' @title Extract or replace parts of an HTSet object
#' @description Operators "$" and [[ extract or replace a slot, while [ extracts
#' a subset of the \code{\link{HTSet-class}} object in a similar manner of
#' \code{\link{subset_features}} and \code{\link{subset_samples}}.
#' @param x \code{\link{HTSet-class}}.
#' @param name The name of the slot.
#' @export
#' @seealso \code{\code{HTSet-class}}
setMethod("$", "HTSet", function(x, name){
    slot(x, name)
})

#' @rdname Extract
#' @param value A replacement value for this slot, which must match the slot
#' definition.
#' @aliases $<-
#' @export
setMethod(
    "$", signature = "HTSet",
    definition = function(x, name){
        slot(x, name)
    }
)

#' @rdname Extract
#' @aliases [[<-
#' @export
setReplaceMethod(
    "$", "HTSet",
    function(x, name, value){
        slot(x, name) = value
        validObject(x)
        return(x)
    }
)

#' @export
setMethod("[[", "HTSet", function(x, i, j, ...){
    stopifnot(is.character(i))
    slot(x, i)
})

#' @rdname Extract
#' @aliases [[<-
#'
#' @export
setReplaceMethod(
    "[[", "HTSet",
    function(x, i, j, ..., value){
        stopifnot(is.character(i))
        slot(x, i) = value
        validObject(x)
        return(x)
    }
)

#' @rdname Extract
#' @aliases [
#' @param i,j indices for features (i) and samples (j) to extract.
#' @param k column index for fdata
#' @param l column index for pdata
#' @export
setMethod(
    "[", "HTSet",
    function(x, i, j, k, l, ..., drop = FALSE) {
        x@edata = x@edata[i,j]
        x@pdata = x@pdata[j,l]
        x@fdata = x@fdata[i,k]
        validObject(x)
        return(x)
    }
)
