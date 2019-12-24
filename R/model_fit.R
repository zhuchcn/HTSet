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
#' @param voom bool. Whether to call \code{\link[limma]{voom}}. Only useful when
#' using limma as the engine.
#' @param edger_model character. Choose from "qlt" for quasi-likelihood negative
#' binomial generalized log-linear model, or "lrt" for negative binomial
#' generalized linear model with likelihood ratio test. See
#' \code{\link[edgeR]{glmQLFit}} and \code{\link[edgeR]{glmLRT}}
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
#'   \item{\strong{logFC}}{estimat of the log2-fold-change corresponding to the
#'   effect.}
#'
#'   \item{\strong{mean}}{average log2-expression for the gene/feature accross
#'   all samples. Same as the AveExpr in limma's \code{\link[limma]{topTable}},
#'   and the baseMean in DESeq2's \code{\link[DESeq2]{results}}.}
#'
#'   \item{\strong{stat}}{the statistic value of the corresponding test. When
#'   using limma, this is the t-statistic value, same as the t column in the
#'   result of \code{\link[limma]{topTable}}. For edgeR, this is the f-statistic
#'   value, same as the F column in the result of edgeR's
#'   \code{\link[edgeR]{topTags}}. As for DESeq2, this is the Wald statistic
#'   value, same as the stat column in the output of DESeq2's
#'   \code{\link[DESeq2]{results}}}
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
#' \item{\strong{egine}}{package used for statistical test}
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
model_fit = function(object, design, coef, engine, transform, voom = FALSE,
                     edger_model = "qlf", adjust.method = "BH"){
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
    if(length(edger_model) != 1){
        stop("edger_model not valid")
    }

    if(!missing(transform)){
        stopifnot(is.function(transform))
        object = transform_by_sample(object, transform)
    }

    if(engine == "limma"){
        res = .fit_limma(object, design, coef, voom, adjust.method)
    } else if(engine == "edgeR"){
        res = .fit_edger(object, design, coef, edger_model, adjust.method)
    } else if(engine == "DESeq2"){
        res = .fit_deseq2(object, design, coef, adjust.method)
    } else {
        stop("The engines supported are limma, edgeR, and DESeq2")
    }
    return(res)
}

#' @keywords internal
.fit_limma = function(object, design, coef, voom, adjust.method){
    if(!requireNamespace("limma")){
        stop("Package limma not installed")
    }
    if(voom){
        if(!requireNamespace("edgeR")){
            stop("Package edgeR not installed")
        }
        y = edgeR::DGEList(object$edata)
        y = edgeR::calcNormFactors(y)
        y = limma::voom(y, design, plot = FALSE)
    } else {
        y = object$edata
    }
    fit = limma::lmFit(y, design)
    fit = limma::eBayes(fit)
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
            adjust.method = adjust.method,
            design = design,
            coef = coef,
            params = list(
                voom = voom
            ),
            engine = "limma"
        ),
        class = "ModelFit"
    )
}

#' @keywords internal
.fit_edger = function(object, design, coef, model, adjust.method){
    if(!requireNamespace("edgeR")){
        stop("Package edgeR not installed")
    }
    y = edgeR::DGEList(counts = object$edata)
    y = edgeR::calcNormFactors(y)
    y = edgeR::estimateDisp(y, design)
    if(model == "qlf"){
        fit = edgeR::glmQLFit(y, design)
        fit = edgeR::glmQLFTest(fit, coef = coef)
    } else if (model == "lrt"){
        fit = edgeR::glmFit(y, design)
        fit = edgeR::glmLRT(fit, coef = coef)
    } else {
        stop("edger_model " %+% model %+% " not valid")
    }
    res = edgeR::topTags(fit, n = Inf, sort.by = "none")
    structure(
        list(
            results = data.frame(
                logFC = res$table$logFC,
                mean  = res$table$logCPM,
                stat  = res$table$F,
                pval  = res$table$PValue,
                padj  = res$table$FDR,
                row.names = rownames(res$table)
            ),
            adjust.method = adjust.method,
            design = design,
            coef = coef,
            params = list(
                edger_model = model
            ),
            engine = "edgeR"
        ),
        class = "ModelFit"
    )
}

#' @keywords internal
.fit_deseq2 = function(object, design, coef, adjust.method){
    if(!requireNamespace("DESeq2")){
        stop("Package DESeq2 not installed")
    }
    de = DESeq2::DESeqDataSetFromMatrix(object$edata, object$pdata, design)
    de = DESeq2::DESeq(de)
    resultName = DESeq2::resultsNames(de)[which(colnames(design) == coef)]
    res = DESeq2::results(de, name = resultName, pAdjustMethod = adjust.method)
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
            adjust.method = adjust.method,
            design = design,
            coef = coef,
            params = list(),
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
