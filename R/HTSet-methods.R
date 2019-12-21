setGeneric("nsamples", function(x) standardGeneric("nsamples"))
setGeneric("nfeatures", function(x) standardGeneric("nfeatures"))
setGeneric("sampleNames", function(x) standardGeneric("sampleNames"))
setGeneric("sampleNames<-", function(x, value) standardGeneric("sampleNames<-"))
setGeneric("featureNames", function(x) standardGeneric("featureNames"))
setGeneric("featureNames<-", function(x, value) standardGeneric("featureNames<-"))

#' @title Get dimemsions
#' @description Get the dimensions of the edata of an HTSet object
#' @seealso \code{\link{HTSet-class}}
#' @param x A \code{\link{HTSet-class}} or derived class object.
#' @return A integer vector with length of 2. The first is the number of
#' features, and the second is number of samples.
#' @export
setMethod(
    "dim", signature = "HTSet",
    function(x){
        return(dim(x@edata))
    }
)

#' @title Get number of samples
#' @description Get the number of samples of a \code{\link{HTSet}} object
#'
#' @param x \code{\link{HTSet-class}}
#'
#' @return integer
#' @seealso \code{\link{HTSet-class}}
#' @export
setMethod(
    "nsamples", signature = "HTSet",
    function(x){
        return(ncol(x@edata))
    }
)

#' @title Get number of features
#' @description Get the number of features of a \code{\link{HTSet}} object
#'
#' @param x \code{\link{HTSet-class}}
#'
#' @return integer
#' @seealso \code{\link{HTSet-class}}
#' @export
setMethod(
    "nfeatures", signature = "HTSet",
    function(x){
        return(nrow(x@edata))
    }
)

#' @rdname sampleNames-HTSet-method
#' @title Get or set the sample names of an HTSet object
#' @description Get or set the sample names of an \code{\link{HTSet}} object.
#'
#' @param x \code{\link{HTSet-class}}
#'
#' @return A character vector.
#' @export
setMethod(
    "sampleNames", signature = "HTSet",
    function(x){
        return(colnames(x@edata))
    }
)

#' @rdname sampleNames-HTSet-method
#' @aliases sampleNames<-
#' @param value A character vector. The length must match number of samples.
#' @export
setReplaceMethod(
    "sampleNames", signature = "HTSet",
    function(x, value){
        colnames(x@edata) = value
        if(!is.null(x@pdata))
            rownames(x@pdata) = value
        validObject(x)
        return(x)
    }
)

#' @rdname featureNames-HTSet-method
#' @title Get or set the feature names of an HTSet object
#' @description Get or set the feature names of an \code{\link{HTSet}} object.
#'
#' @param x \code{\link{HTSet-class}}
#'
#' @return A character vector.
#' @export
setMethod(
    "featureNames", signature = "HTSet",
    function(x){
        return(rownames(x@edata))
    }
)

#' @rdname featureNames-HTSet-method
#' @aliases featureNames<-
#' @param value A character vector. The length must match number of features.
#' @export
setReplaceMethod(
    "featureNames", signature = "HTSet",
    function(x, value){
        rownames(x@edata) = value
        if(!is.null(x@fdata))
            rownames(x@fdata) = value
        validObject(x)
        return(x)
    }
)

#' @title Subset samples
#' @description Get a subset of samples of an HTSet object.
#' @param x \code{\link{HTSet-class}}
#' @param samples The samples to subset. Can be character, logical, or integers.
#' @return \code{\link{HTSet-class}}
#' @export
subset_samples = function(object, samples){
    stopifnot(is(object, "HTSet"))
    msg = "The \"samples\" provided not valid. Please verify."
    if(is.character(samples)){
        if(any(! samples %in% sampleNames(object)))
            stop(msg)
    }
    if(is.numeric(samples)){
        if(max(samples) > nsamples(object))
            stop(msg)
    }
    if(is.logical(samples)){
        if(length(samples) != nsamples(object))
            stop(msg)
    }
    edata = object@edata[, samples]
    pdata = object@pdata[samples, ]
    new_obj = HTSet(edata, object@fdata, pdata, object@assay)
    return(new_obj)
}

#' @title Subset features
#' @description Get a subset of features of an HTSet object.
#' @param x \code{\link{HTSet-class}}
#' @param features The features to subset. Can be character, logical, or integers.
#' @return \code{\link{HTSet-class}}
#' @export
subset_features = function(object, features){
    stopifnot(is(object, "HTSet"))
    msg = "The \"features\" provided not valid. Please verify."
    if(is.character(features)){
        if(any(! features %in% featureNames(object)))
            stop(msg)
    }
    if(is.numeric(features)){
        if(max(features) > nfeatures(object))
            stop(msg)
    }
    if(is.logical(features)){
        if(length(features) != nfeatures(object))
            stop(msg)
    }
    edata = object@edata[features, ]
    fdata = object@fdata[features, ]
    new_obj = HTSet(edata, fdata, object@pdata, object@assay)
    return(new_obj)
}

#' @title Summarize by feature variable
#' @description This function summarize the edata according to the given feature
#' variable.
#' @param object \code{\link{HTSet-class}}
#' @param feature_var character. Must come from the column names of the
#' fdata. Length equals to 1.
#' @return \code{\link{HTSet-class}}
#' @export
#' @author Chenghao Zhu
#' @importFrom magrittr %>%
#' @import dplyr
#' @import reshape2
summarize_feature = function(object, feature_var){
    stopifnot(is(object, "HTSet"))
    if(length(feature_var) != 1)
        stop("Invalid feature_var")
    if(!is.character(feature_var))
        stop("Invalid feature_var")
    if(!feature_var %in% colnames(object@fdata))
        stop(paste0("The feature_var '", feature_var, "' does not exist"))

    edata = object@edata %>%
        as.data.frame %>%
        mutate(feature_var = object@fdata[,feature_var]) %>%
        melt(id.var = "feature_var",
             variable.name = "sample_id",
             value.name = "expr") %>%
        group_by(feature_var, sample_id) %>%
        summarize(expr = sum(expr, na.rm = TRUE)) %>%
        dcast(feature_var~sample_id, value.var = "expr")
    rownames(edata) = edata$feature_var
    edata = select(edata, -feature_var) %>%
        as.matrix

    object@edata = edata
    object@fdata = NULL
    validObject(object)
    return(object)
}

#' @title Transform by sample
#' @description Transform the edata in a sample-wise manner (MARGIN = 2).
#' @param object \code{\link{HTSet-class}}
#' @param FUN function to apply.
#' @param ... additional arguments to pass to the function
#' @return \code{\link{HTSet-class}}
#' @export
transform_by_sample = function(object, FUN, ...){
    if(length(FUN(1:nfeatures(object))) != nfeatures(object))
        stop("Function invalid")

    edata = apply(object@edata, 2, FUN, ...)

    rownames(edata) = featureNames(object)
    colnames(edata) = sampleNames(object)

    object@edata = edata
    validObject(object)
    return(object)
}

#' @title Transform by feature
#' @description Transform the edata in a feature-wise manner (MARGIN = 2).
#' @param object \code{\link{HTSet-class}}
#' @param FUN function to apply.
#' @param ... additional arguments to pass to the function
#' @return \code{\link{HTSet-class}}
#' @export
transform_by_feature = function(object, FUN, ...){
    if(length(FUN(1:nsamples(object))) != nsamples(object))
        stop("Function invalid")

    edata = t(apply(object@edata, 1, FUN, ...))

    rownames(edata) = featureNames(object)
    colnames(edata) = sampleNames(object)

    object@edata = edata
    validObject(object)
    return(object)
}
