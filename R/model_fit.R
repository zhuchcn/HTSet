#' @keywords internal
setDefault = function(x, key, args){
    if(is.null(x[[key]]))
        x[[key]] = list()
    if(!is.list(x[[key]]))
        stop("args must be a list")
    for(k in names(args)){
        if(is.null(x[[key]][[k]]))
            x[[key]][[k]] = args[[k]]
    }
    return(x)
}

#' @keywords internal
#' @title validate args in options
#' @description validate args for \code{\link{model_fit}} in given options
#' @param args list
#' @param opts character
validateArgsInOptions = function(args, opts){
    if(any(! names(args) %in% opts)){
        if (length(opts) == 1){
            msg = "any argument in \"args\" not equal to \"" %+% opts %+% "\" is ignored"
        } else {
            str_opts = "\"" %+% paste(opts[1:length(opts)-1], collapse = "\", \"") %+%
                "\" or \"" %+% opts[length(opts)] %+% "\""
            msg = "any argument in \"args\" not equal to " %+% str_opts %+% " is ignored"
        }
        warning(msg)
    }
    return(0)
}

#' @title statistic tests for HTSet
#' @description This is a wrapper of several widely used statistical method for
#' high through-put experimental data such as RNAseq. The
#' \code{\link[limma]{limma-package}} performs linear model on continous data,
#' or cooperate with the \code{\link[limma]{voom}} function to handle count
#' data. The \code{\link[edgeR]{edgeR-package}} and
#' \code{\link[DESeq2]{DESeq2-package}} performs negative binomial generalized
#' linear model on count data.
#'
#' @param object \code{\link{HTSet}}
#' @param design matrix. Number of rows must batch number of samples of object.
#' Usually the output of \code{\link{model.matrix}}.
#' @param coef character. The coefficient to perform statistical tests. Must
#' be in the column names of design.
#' @param engine character. The engine to perform statistical analysis.
#' Supported are limma, edgeR, and DESeq2.
#' @param args list. A list of argumnets to be parsed to the backend statistical
#' engine.
#' @param adjust.method character. Method used to adjust the p-values for
#' multiple testing. See \code{\link{p.adjust}} for the complete list of
#' supported methods. Default is "BH" for Benjamini-Hochberg test.
#' @param transform function. The transform to be passed to
#' \code{\link{transform_by_sample}}. No transform is performed by default.
#'
#' @return A list-like S3 class ModelFit object is returned with the elements
#' as following.
#'
#' \item{\strong{results}}{
#'   A data.frame of the statistical test results for each gene/
#'   feature.
#'
#'   \describe{
#'
#'   \item{\strong{logFC}}{estimate of the log2-fold-change corresponding to the
#'   effect.}
#'
#'   \item{\strong{mean}}{average log2-expression for the gene/feature accross
#'   all samples. Same as the AveExpr in limma's \code{\link[limma]{topTable}},
#'   and the baseMean in DESeq2's \code{\link[DESeq2]{results}}.}
#'
#'   \item{\strong{stat}}{the statistic value of the corresponding test. When
#'   using limma, this is the t-statistic value, same as the t column in the
#'   result of \code{\link[limma]{topTable}}. For edgeR, this is the f-statistic
#'   value for Quansi-likelihood test or the likelihood ratio statistic value
#'   for likelihood ratio test. Same as the column F or RT in the result of
#'   edgeR's \code{\link[edgeR]{topTags}}. As for DESeq2, this is the Wald
#'   statistic value for Wald test or the difference in deviance between the
#'   reduced model and the full model for likelihood ratio test. same as the
#'   stat column in the output of DESeq2's \code{\link[DESeq2]{results}}.}
#'
#'   \item{\strong{pval}}{raw p-value}
#'
#'   \item{\strong{padj}}{adjusted p-value}
#'   }
#' }
#'
#' \item{\strong{adjust.method}}{Method used to correct for multiple testing}
#'
#' \item{\strong{design}}{design matrix}
#'
#' \item{\strong{df}}{degree of freedoms}
#'
#' \item{\strong{distribution}}{The distribution that p values were calculated}
#'
#' \item{\strong{engine}}{package used for statistical test}
#'
#' \item{\strong{coef}}{coefficient tested}
#'
#' \item{\strong{params}}{additional parameters parsed}
#'
#' @seealso
#' \code{\link{HTSet-class}}
#' \code{\link[limma]{lmFit}}
#' \code{\link[limma]{voom}}
#' \code{\link[edgeR]{glmQLFit}}
#' \code{\link[edgeR]{glmLRT}}
#' \code{\link[DESeq2]{DESeq}}
#'
#' @importFrom stats p.adjust.methods
#' @export
#' @author Chenghao Zhu
#'
#' @examples
#' data(exrna)
#' design = model.matrix(~ Condition, data = exrna$pdata)
#' coef = "ConditionSystemic Lupus Erythematosus"
#' fit1 = model_fit(object = exrna, design = design, coef = coef, engine = "limma", args = list(voom = TRUE))
#' fit2 = model_fit(object = exrna, design = design, coef = coef, engine = "edgeR")
#' fit3 = model_fit(object = exrna, design = design, coef = coef, engine = "edgeR", args = list(model = "lrt"))
#' fit4 = model_fit(object = exrna, design = design, coef = coef, engine = "DESeq2")
model_fit = function(object, design, coef, engine, args = list(), transform,
                     adjust.method = "BH"){
    if(!is(object, "HTSet")){
        stop("model_fit only works for HTSet objects.")
    }
    if(length(coef) != 1){
        stop("Please provide only one coef.")
    }
    if(!coef %in% colnames(design)){
        stop("coef " %+% coef %+% " not found in the design matrix")
    }
    if(length(engine) != 1){
        stop("engine number > 1")
    }
    if(nrow(design) != nsamples(object)){
        stop("design matrix invalid")
    }
    if(length(adjust.method) != 1){
        stop("adjust.method not valid")
    }
    if(!adjust.method %in% p.adjust.methods){
        stop("adjust.method not valid")
    }

    if(!missing(transform)){
        stopifnot(is.function(transform))
        object = transform_by_sample(object, transform)
    }

    .fit = switch(
        engine,
        "limma"  = .fit_limma,
        "edgeR"  = .fit_edger,
        "DESeq2" = .fit_deseq2,
        stop("The engines supported are limma, edgeR, and DESeq2")
    )
    res = .fit(object, design, coef, adjust.method, args)
    return(res)
}

#' @keywords internal
.fit_limma = function(object, design, coef, adjust.method, args){
    .args = args
    if(!requireNamespace("limma"))
        stop("Package limma not installed")
    validateArgsInOptions(args, c("voom", "lmFit", "eBayes"))
    if(identical(args$voom, TRUE))
        args$voom = list()
    if(identical(args$voom, FALSE))
        args$voom = NULL

    if(is.null(args$voom)){
        y = object$edata
    } else {
        if(!is.list(args$voom))
            stop("args must be a nested list")
        if(!requireNamespace("edgeR"))
            stop("Package edgeR not installed")
        y = edgeR::DGEList(object$edata)
        y = edgeR::calcNormFactors(y)
        args$voom$counts = y
        args$design = design
        y = do.call(limma::voom, args$voom)
    }
    args = setDefault(args, "lmFit", list(object = y, design = design))
    fit = do.call(limma::lmFit, args$lmFit)

    args = setDefault(args, "eBayes", list(fit = fit))
    fit = do.call(limma::eBayes, args$eBayes)

    res = limma::topTable(fit, coef = coef, number = Inf, sort.by = "none",
                          adjust.method = adjust.method)
    structure(
        list(
            results = data.frame(
                logFC = res$logFC,
                mean  = res$AveExpr,
                stat  = res$t,
                pval  = res$P.Value,
                padj  = res$adj.P.Val,
                row.names = rownames(res)
            ),
            df = fit$df.total,
            distribution = "t",
            adjust.method = adjust.method,
            design = design,
            coef = coef,
            params = .args,
            engine = "limma"
        ),
        class = "ModelFit"
    )
}

#' @keywords internal
.fit_edger = function(object, design, coef, adjust.method, args){
    .args = args
    validateArgsInOptions(args, c("calcNormFactors", "estimateDisp", "model"))
    if(!requireNamespace("edgeR")){
        stop("Package edgeR not installed")
    }
    y = edgeR::DGEList(counts = object$edata)

    args = setDefault(args, "calcNormFactors", list(object = y))
    y = do.call(edgeR::calcNormFactors, args$calcNormFactors)

    args = setDefault(args, "estimateDisp", list(y = y, design = design))
    y = do.call(edgeR::estimateDisp, args$estimateDisp)

    if(is.null(args$model))
        args$model = "qlf"
    if(length(args$model) > 1)
        stop("args$model invalid")
    if(args$model == "qlf"){
        args = setDefault(args, "glmQLFit", list(y = y, design = design))
        fit = do.call(edgeR::glmQLFit, args$glmQLFit)

        args = setDefault(args, "glmQLFTest", list(glmfit = fit, coef = coef))
        fit = do.call(edgeR::glmQLFTest, args$glmQLFTest)
        # from edgeR::glmQLFTest
        df = list(
            df1 = fit$df.test,
            df2 = fit$df.total
        )
        distribution = "f"
    } else if (args$model == "lrt"){
        args = setDefault(args, "glmFit", list(y = y , design = design))
        fit = do.call(edgeR::glmFit, args$glmFit)

        args = setDefault(args, "glmLRT", list(glmfit = fit, coef = coef))
        fit = do.call(edgeR::glmLRT, args$glmLRT)

        df = fit$df.test
        distribution = "chisq"
    } else {
        stop("edgeR model " %+% args$model %+% " not valid")
    }
    res = edgeR::topTags(fit, n = Inf, sort.by = "none")
    structure(
        list(
            results = data.frame(
                logFC = res$table$logFC,
                mean  = res$table$logCPM,
                stat  = {if(args$model == "qlf") res$table$F else res$table$LR},
                pval  = res$table$PValue,
                padj  = res$table$FDR,
                row.names = rownames(res$table)
            ),
            df = df,
            distribution = distribution,
            adjust.method = adjust.method,
            design = design,
            coef = coef,
            params = .args,
            engine = "edgeR"
        ),
        class = "ModelFit"
    )
}

#' @keywords internal
.fit_deseq2 = function(object, design, coef, adjust.method, args){
    .args = args
    validateArgsInOptions(args, c("DESeq"))
    if(!requireNamespace("DESeq2")){
        stop("Package DESeq2 not installed")
    }
    de = DESeq2::DESeqDataSetFromMatrix(object$edata, object$pdata, design)

    args = setDefault(args, "DESeq", list(object = de, test = "Wald", useT = FALSE))
    de = do.call(DESeq2::DESeq, args$DESeq)
    resultName = DESeq2::resultsNames(de)[which(colnames(design) == coef)]
    res = DESeq2::results(de, name = resultName, pAdjustMethod = adjust.method)

    if(args$DESeq$test == "Wald"){
        df = de@rowRanges@elementMetadata$tDegreesFreedom
        if(args$DESeq$useT){
            distribution = 't'
        } else {
            distribution = "norm"
        }
    } else {
        df = rep(ncol(design - ncol(args$DESeq$reduced)), nfeatures(object))
        distribution = "chisq"
    }
    structure(
        list(
            results = data.frame(
                logFC = res$log2FoldChange,
                mean  = res$baseMean,
                stat  = res$stat,
                pval  = res$pvalue,
                padj  = res$padj,
                row.names = rownames(res)
            ),
            df = df,
            distribution = distribution,
            adjust.method = adjust.method,
            design = design,
            coef = coef,
            params = .args,
            engine = "DESeq2"
        ),
        class = "ModelFit"
    )
}

#' @export
volcanoplot = function(x, ...){
    UseMethod("volcanoplot", x)
}

#' @title volcano plot
#' @description Creates a volcano plot for the returned result of
#' \code{\link{model_fit}}
#' @param x ModelFit. Must be returned by \code{\link{model_fit}}
#' @param hline numeric. The p-value where the horizontal lines to draw. All
#' values must between 0 and 1. If not given, serious lines will be drown at
#' 0.05, 0.01, 0.001, 0.0001 ...
#' @param color character for point color
#' @param alpha numeric for point alpha
#' @param hline_type character for hline type. Default is dashed
#' @param hline_color chracter for hline color.
#' @param hline_size numeric for hline size
#' @param hline_alpha numeric for hline alpha
#' @param ... other args
#' @return ggplot
#' @export
#' @import ggplot2
#' @aliases volcanoplot
#' @author Chenghao Zhu
volcanoplot.ModelFit = function(
    x, hline, color = "grey10", alpha = 0.75, hline_type = "dashed",
    hline_color = "grey30", hline_size = 0.5, hline_alpha = 0.75, ...
){
    p = ggplot(x$results)
    if(missing(hline)){
        p = p + geom_hline(yintercept = -log(0.05), linetype = hline_type,
                           color = hline_color, size = hline_size,
                           alpha = hline_alpha)
        yint = 0.01
        while(yint > min(x$results$pval)){
            p = p + geom_hline(yintercept = -log(yint), linetype = hline_type,
                               color = hline_color, size = hline_size,
                               alpha = hline_alpha)
            yint = yint / 10
        }
    } else {
        if(!is.numeric(hline)){
            stop("invalid hline")
        }
        if(any(!between(hline, 0, 1))){
            stop("invalid hline")
        }
        for(h in hline){
            p = p + geom_hline(yintercept = -log(h), linetype = hline_type,
                               color = hline_color, size = hline_size,
                               alpha = hline_alpha)
        }
    }
    p + geom_point(aes(x = logFC, y = -log(pval)), color = color, alpha = alpha) +
        theme_classic() +
        labs(title = x$coef) +
        theme(
            plot.title = element_text(hjust = 0.5)
        )
}

#' @title p value histogram
#' @description create a histogram for p values of model fit result returned by
#' \code{\link{model_fit}}
#'
#' @param x ModelFit. Must be returned by \code{\link{model_fit}}
#' @param vline p value where vertical line to draw
#' @param bins numeric. Total number of bins
#' @param color character for bin outline colors
#' @param fill character for bin fill colors
#' @param vline_linetype character for vline linetype
#' @param vline_size numeric for vline size
#' @param vline_color character for vline color
#' @param ... other args
#'
#' @return ggplot
#' @import ggplot2
#' @export
#' @author Chenghao Zhu
hist.ModelFit = function(x, vline, bins = 40, color = "white", fill = "gray10",
                         vline_linetype = "dashed", vline_size = 0.5,
                         vline_color = "gray30", ...){
    p = ggplot(x$results) +
        geom_histogram(aes(pval), bins = bins, color = color, fill = fill)
    if(!missing(vline)){
        if(!is.numeric(vline))
            stop("vline invalid")
        if(length(vline) != 1)
            stop("vline invalid")
        if(!between(vline, 0, 1))
            stop("vline invalid")
        p = p + geom_vline(xintercept = vline, color = vline_color,
                           linetype = vline_linetype, size = vline_size)
    }
    p + theme_classic()
}
