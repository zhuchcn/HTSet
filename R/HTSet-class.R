setClassUnion("listOrNULL", c("NULL", "list"))
setClassUnion("dataframeOrNULL", c("NULL", "data.frame"))

#' @title HTSet
#'
#' @description S4 class for storing high through-put experiment data. It has
#' three enssential components, the edata, fdata, and pdata. The edata is the
#' expression/abundance matrix of each feature in each sample. The features are
#' in rows and samples are in columns. The fdata is the feature data which
#' holds the information for each feature, such as gene names, protein IDs, or
#' metabolite names. The pdata is the phenotype data that describe additional
#' information for each sample, such as sex, condition, or treatment. The rows
#' in the edata must match the rows in the fdata, and the columns in the edata
#' must match the rows in the pdata.
#'
#' @slot edata matrix. Must be numeric (integer or float)
#' @slot fdata data.frame. The number of rows and row names must match those
#' of the edata.
#' @slot pdata data.frame. The number of rows and row names must match the
#' number of columns and column names of edata.
#' @slot assay list. Additional assay or study information.
#'
#' @exportClass HTSet
#' @author Chenghao Zhu
setClass(
    "HTSet",
    representation = representation(
        edata = "matrix",
        fdata = "dataframeOrNULL",
        pdata = "dataframeOrNULL",
        assay = "listOrNULL"
    ),
    validity = function(object){
        if(!is.numeric(object@edata)){
            stop("edata must be numeric.")
        }
        if(!is.null(object@fdata)){
            if (nrow(object@edata) != nrow(object@fdata)){
                stop("The row number of edata must match fdata.")
            }
            if (!all.equal(rownames(object@edata), rownames(object@fdata))){
                stop("The row names of edata must match fdata.")
            }
        }
        if(!is.null(object@pdata)){
            if (ncol(object@edata) != nrow(object@pdata)) {
                stop("The column number of edata must match the row number of pdata.")
            }
            if (!all.equal(colnames(object@edata), rownames(object@pdata))) {
                stop("The column nuames of edata must match the row names of pdata.")
            }
        }
        return(TRUE)
    }
)

#' @export
HTSet = function(edata, fdata = NULL, pdata = NULL, assay = NULL){
    new("HTSet", edata = edata, fdata = fdata, pdata = pdata, assay = assay)
}

#' @keywords internal
`%+%` = function(x, y) paste0(x,y)

#' @export
setMethod(
    "show", signature = "HTSet",
    definition = function(object){
        if(is.null(object@pdata)){
            nsamplevars = 0
        } else {
            nsamplevars = ncol(object@pdata)
        }
        if(is.null(object@fdata)){
            nfeaturevars = 0
        } else {
            nfeaturevars = ncol(object@fdata)
        }
        if(requireNamespace("crayon")){
            blue = crayon::blue
        } else {
            blue = cat
        }
        cat(
            blue("S4 class: HTSet") %+% "\n" %+%
            "  nsamples: " %+% nsamples(object) %+% "  \tnfeatures:" %+% nfeatures(object) %+% "\n" %+%
            "  sample vars: " %+% nsamplevars %+% "\tfeatures vars: " %+% nfeaturevars %+% "\n" %+%
            blue("slots:") %+% "\n" %+%
            "  $edata: " %+% nrow(object@edata) %+% "x" %+% ncol(object@edata) %+% " matrix\n" %+%
            "  $fdata: " %+% nrow(object@fdata) %+% "x" %+% ncol(object@fdata) %+% " data.frame\n" %+%
            "  $pdata: " %+% nrow(object@pdata) %+% "x" %+% ncol(object@pdata) %+% " data.frame\n" %+%
            "  $assay: " %+% length(object@assay) %+% " list" %+% "\n"
        )
    }
)
