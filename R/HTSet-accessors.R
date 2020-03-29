#' @rdname HTSet-Extract
#' @aliases $
#' @title Extract or replace parts of an HTSet object
#' @description Operators "$" and [[ extract or replace a slot, while [ extracts
#' a subset of the \code{\link{HTSet-class}} object in a similar manner of
#' \code{\link{subset_features}} and \code{\link{subset_samples}}.
#' @param x \code{\link{HTSet-class}}.
#' @param name The name of the slot.
#' @export
#' @seealso \code{\link{HTSet-class}}
setMethod("$", "HTSet", function(x, name){
    slot(x, name)
})

#' @rdname HTSet-Extract
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

.DollarNames.HTSet = function(x, pattern) {
    grep(pattern, slotNames(x), value = TRUE)
}


#' @rdname HTSet-Extract
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

#' @rdname HTSet-Extract
#' @aliases [[
#' @aliases [[,HTSet,character-method
#' @export
setMethod(
    "[[",
    signature = signature(x = "HTSet", i = "character", j = "ANY"),
    function(x, i, j, ...){
        slot(x, i)
    }
)

#' @rdname HTSet-Extract
#' @aliases [[<-
#' @aliases [[<-,HTSet,character-method
#' @export
setReplaceMethod(
    "[[", signature = signature(x = "HTSet", i = "character", j = "ANY"),
    function(x, i, j, ..., value){
        stopifnot(is.character(i))
        slot(x, i) = value
        validObject(x)
        return(x)
    }
)

#' @rdname HTSet-Extract
#' @aliases [,HTSet-method
#' @aliases [
#' @param i,j indices for features (i) and samples (j) to extract.
#' @param k column index for fdata
#' @param l column index for pdata
#' @param drop If TRUE the result is coerced to the lowest possible dimension
#' @param ... other arguments
#' @export
setMethod(
    "[", signature = signature(x = "HTSet", i = "ANY", j = "ANY"),
    function(x, i, j, k, l, ..., drop = FALSE) {
        x@edata = x@edata[i,j,drop = drop]
        x@pdata = x@pdata[j,l,drop = drop]
        x@fdata = x@fdata[i,k,drop = drop]
        validObject(x)
        return(x)
    }
)
