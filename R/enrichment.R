#' @title Enrichment test
#' @description Perform enrichment test of genes/features in categories with
#' Fisher's exact test or Kolmogorov-Smirnov Tests.
#'
#' @details
#' For Fisher's exact test, the question being answered is, whether the number
#' of genes/features in category A that are increased (or decreased) in
#' treatment different from the background. In other words, the Fisher's exact
#' test is used to test the independency of gene/feature category A and whether
#' genes/features are increased (or decreased).
#'
#' For Kilmogorov-Smirnov test, a similar hypothesis is tested. The p-values
#' are first tranformed into one-tail. The transformation is done using the
#' statistic values returned by \code{\link{model_fit}}. The transformed
#' p-values are then compared to an arithmetic sequence between 0 and 1 using
#' \code{\link{ks.test}}.
#'
#' @param object \code{\link{HTSet-class}}
#' @param fit ModelFit. The result of \code{\link{model_fit}}
#' @param group character. The fdata column of object that will be used as the
#' categories to perform enrichment analysis.
#' @param test character. Must be "fet" for Fisher's exact test or "kst" for
#' Kolmogorov-Smirnov test.
#' @param alternative character. One of both, greater, less, or two.sided.
#' Default is both, that both greater and less were calculated and returned.
#' @param p.cutoff numeric. The p-value cutoff to define increase or decreasing.
#' Only useful with fet. The default value of this argument is 0.05 when
#' alternative equals to two.sided, and is 1 otherwise.
#'
#' @return A list-like S3 object. \code{pval} is the pvalues for enrichment of
#' each category. For fet, the \code{odds.ratio} has the odds ratio for each
#' category, and the \code{matrix} has the N, m, n, k, and x values. For kst,
#' \code{adj.p.val} has the transformed p-values under each of the alternative
#' hypothesis of greater or less.
#'
#' @examples
#' design = model.matrix(~ Condition, data = exrna$pdata)
#' coef = "ConditionSystemic Lupus Erythematosus"
#' fit = model_fit(exrna, design, coef, engine = "limma", args = list(voom = TRUE))
#' en = enrichment_test(exrna, fit, "gene_type", "fet")
#' en = enrichment_test(exrna, fit, "gene_type", "kst")
#'
#' @seealso
#' \code{\link{HTSet-class}}
#' \code{\link{model_fit}}
#' \code{\link{fisher.test}}
#' \code{\link{phyper}}
#' \code{\link{ks.test}}
#'
#' @importFrom stats fisher.test
#' @importFrom stats ks.test
#' @export
#' @author Chenghao Zhu
enrichment_test = function(object, fit, group, test = c("fet", "kst"),
                           alternative = c("both", "greater", "less", "two.sided"),
                           p.cutoff){
    stopifnot(is(object, "HTSet"))
    stopifnot(is(fit, "ModelFit"))
    if(nfeatures(object) != nrow(fit$results))
        stop("object and fit must match")
    if(!all.equal(featureNames(object), rownames(fit$results)))
        stop("object and fit must match")
    if(length(test) != 1)
        stop("invalid test")
    test = match.arg(test, c("fet", "kst"))
    alternative = match.arg(alternative, c("both", "greater", "less", "two.sided"))
    if(missing(p.cutoff))
        p.cutoff = if(alternative == "two.sided") 0.05 else 1
    if(length(p.cutoff) != 1)
        stop("invalid p.cutoff")
    if(!between(p.cutoff, 0, 1))
        stop("invalid p.cutoff")
    if(!is.character(group))
        stop("invalid group")
    if(length(group) != 1)
        stop("invalid group")
    if(!group %in% colnames(object$fdata))
        stop("group can not been found")

    group = object$fdata[, group]
    res = {
        if(test == "fet"){
            .fet_enrichment(object, fit, group, alternative, p.cutoff)
        } else {
            .kst_enrichment(object, fit, group, alternative)
        }
    }
    return(res)
}

#' @keywords internal
.fet_enrichment = function(object, fit, group, alternative, p.cutoff){

    if(alternative == "both"){
        res = lapply(c("greater", "less"),
                     function(alt) .fet_wrapper(fit, group, alt, p.cutoff))
        names(res) = c("greater", "less")
        pval = sapply(res, function(r) r[, "p.value"])
        odds.ratio = sapply(res, function(r) r[, "odds.ratio"])
        matrix = lapply(res, function(r) r[,1:5])
    } else {
        res = .fet_wrapper(fit, group, alternative, p.cutoff)
        pval = res[, "p.value"]
        odds.ratio = res[, "odds.ratio"]
        matrix = res[, 1:5]
    }

    structure(
        list(
            pval = pval,
            odds.ratio = odds.ratio,
            matrix = matrix,
            p.cutoff = p.cutoff,
            fit = fit,
            group = group,
            alternative = alternative
        ),
        class = "EnrichmentFET"
    )
}

#' @keywords internal
.fet_wrapper = function(fit, group, alt, p.cutoff){
    if(alt != "two.sided"){
        .compare = switch(
            alt,
            "greater" = function(x){x > 0},
            "less" = function(x){x < 0}
        )
    }
    lapply(unique(group), function(ele){
        N = nrow(fit$results)
        m = sum(group == ele)
        n = N - m
        if(alt == "two.sided"){
            k = sum(fit$results$pval < p.cutoff)
            x = sum(fit$results$pval < p.cutoff & group == ele)
        } else {
            k = sum(.compare(fit$results$logFC) & fit$results$pval < p.cutoff)
            x = sum(.compare(fit$results$logFC) & fit$results$pval < p.cutoff & group == ele)
        }

        # http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
        fet = fisher.test(matrix(c(x, k-x, m-x, N-m-(k-x)), 2, 2), alternative = "greater")
        p.value = fet$p.value
        odds.ratio = fet$estimate
        names(odds.ratio) = NULL
        return(c(N = N, m = m, n = n, k = k, x = x, p.value = p.value, odds.ratio = odds.ratio))
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(unique(group))
}

#' @keywords internal
.kst_enrichment = function(object, fit, group, alternative){
    if(fit$engine == "DESeq2"){
        if(!is.null(fit$params$DESeq)){
            if(!is.null(fit$params$DESeq$test)){
                if(fit$params$DESeq$test == "LRT")
                    stop("DESeq2's LRT test is not supported for ks test yet.")
            }
        }
    }
    if(alternative == "both"){
        res = lapply(c("greater", "ess"),
                     function(alt) .kst_wrapper(fit, group, alt))
        names(res) = c("greater", "less")
        pval = sapply(res, function(r) r$pvalues)
        d = sapply(res, function(r) r$d)
        fit_pvalues = data.frame(
            raw = fit$results$pval,
            greater = res$greater$adjusted.p.values,
            less = res$less$adjusted.p.values,
            row.names = featureNames(object)
        )
    } else {
        res = .kst_wrapper(fit, group, alternative)
        pval = res$pvalues
        d = res$d
        fit_pvalues = {
            if(alternative == "two.sided")
                data.frame(
                    raw = fit$results$pval,
                    row.names = featureNames(object)
                )
            else
                data.frame(
                    raw = fit$results$pval,
                    adjusted = res$adjusted.p.values,
                    row.names = featureNames(object)
                )
        }
    }

    structure(
        list(
            pval = pval,
            d = d,
            fit = fit,
            group = group,
            fit_pvalues = fit_pvalues,
            alternative = alternative
        ),
        class = "EnrichmentKST"
    )
}

#' @keywords internal
.kst_wrapper = function(fit, group, alt){
    pvalues = {
        if(alt != "two.sided")
            .transform_pvalues(fit$results$stat, fit$distribution, fit$df, alt == "less")
        else
            fit$results$pval
    }
    res = lapply(unique(group), function(ele){
        pvals = pvalues[group == ele]
        ks = ks.test(
            pvals, seq(from = 0, to = 1, length.out = length(pvals)),
            alternative = "greater"
        )
        d = ks$statistic
        names(d) = NULL
        c("d" = d, "p" = ks$p.value)
    })
    res = do.call(rbind, res)
    rownames(res) = unique(group)
    names(pvalues) = rownames(fit)
    list(pvalues = res[, "p"], "d" = res[, "d"], adjusted.p.values = pvalues)
}

#' @keywords internal
#' @importFrom stats pt
#' @importFrom stats pnorm
#' @importFrom stats pf
#' @importFrom stats pchisq
.transform_pvalues = function(stat, dist, df, lower.tail){
    if(dist == "t"){
        pt(q = stat, df = df, lower.tail = lower.tail)
    } else if(dist == "norm"){
        pnorm(q = stat, lower.tail = lower.tail)
    } else if(dist == "f"){
        # from edgeR::glmQLFTest
        pf(q = stat, df1 = df$df1, df2 = df$df2, lower.tail = lower.tail)
    } else if(dist == "chisq"){
        pchisq(q = stat, df = df, lower.tail = lower.tail)
    }
}

#' @rdname barplot
#' @export
barplot = function(x, ...){
    UseMethod("barplot", x)
}

#' @aliases barplot
#' @rdname barplot
#' @title barplot for fet enrichment test
#' @description This function makes bar plot for top most positively and negatively
#' enriched levels.
#'
#' @param x EnrichmentFET. Returned by \code{\link{enrichment_test}}
#' with test equals to fet.
#' @param num integer. Number of categories to display. Use it when alternative
#' equals to greater, less, or two.sided.
#' @param greater integer. Number of categories which are greater. Use it when
#' alternative equals to both.
#' @param less integer. Number of categories which are less. Use it when
#' alternative equals to both.
#' @param text.size integer
#' @param fill character. Bar fill color
#' @param panel.ratio numeric. The ratio of left and right panels.
#' @param title character. The title of the plot.
#' @param ... other params not supported
#'
#' @return a ggplot object
#' @examples
#' lpd = transform_by_sample(lipidome, function(x) log(x/sum(x)))
#' design = model.matrix(~Treatment * Timepoint + Subject, data = lpd$pdata)
#' fit = model_fit(lpd, design, "TreatmentMed:TimepointPre", "limma")
#' en = enrichment_test(lpd, fit, "class", "fet")
#' barplot(en)
#' en = enrichment_test(lpd, fit, "class", "fet", alternative = "greater")
#' barplot(en)
#' @export
#' @import ggplot2
#' @import dplyr
#' @import cowplot
#' @author Chenghao Zhu
#' @seealso \code{\link{enrichment_test}}
barplot.EnrichmentFET = function(x, num = 6, greater = 3, less = 3, text.size = 10,
                              fill = "gray20", panel.ratio = 1.2, title, ...){
    df = x$matrix %>% as.data.frame
    df$var = rownames(df)
    my_theme = theme_classic() +
        theme(
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = text.size + 4)
        )
    if(x$alternative == "both"){
        if(missing(title))
            title = c("higher", "lower")
        df = df %>%
            mutate(diff.x = greater.x - less.x) %>%
            arrange(diff.x) %>%
            filter(var %in% var[1:less] | var %in% rev(var)[1:greater]) %>%
            mutate(var = factor(var, levels = var))
        ymax = max(df$greater.x, df$less.x)
        p1 = ggplot(df) +
            geom_col(aes(x = var, y = greater.x), fill = fill, width = .75) +
            geom_text(
                aes(x = var, y = greater.x + ymax/20, label = greater.x),
                size = text.size / 2.83
            ) +
            coord_flip(ylim = c(0, ymax * 1.2)) +
            labs(x = "", title = title[1]) +
            my_theme +
            theme(axis.text.y = element_text(
                color = "black", hjust = 0.5, margin = margin(0, 10, 0, 0),
                size = text.size
            ))
        p2 = ggplot(df) +
            geom_col(aes(x = var, y = less.x), fill = fill, width = .75) +
            geom_text(aes(x = var, y = less.x + ymax/20, label = less.x),
                      size = text.size / 2.83) +
            scale_y_reverse() +
            scale_x_discrete(position = "top") +
            coord_flip(ylim = c(0, ymax * 1.2)) +
            labs(x = "", title = title[2]) +
            my_theme +
            theme(axis.text.y = element_blank())
        p = plot_grid(p2, p1, ncol = 2, rel_widths = c(1, panel.ratio))
    } else {
        if(num > nrow(df))
            num = nrow(df)
        df = arrange(df, desc(x)) %>%
            head(num) %>%
            mutate(var = factor(var, levels = rev(var)))
        p = ggplot(df) +
            geom_col(aes(var, x), fill = fill, width = .75) +
            geom_text(
                aes(x = var, y = x + max(df$x)/20, label = x),
                size = text.size / 2.83
            ) +
            coord_flip() +
            my_theme +
            theme(axis.text.y = element_text(
                color = "black", hjust = 0.5, margin = margin(0, 10, 0, 0),
                size = text.size
            ))
        if(!missing(title)){
            if(length(title) > 1){
                warning("the first title is used")
                title = title[1]
            }
            p = p + ggtitle(title)
        }
    }
    return(p)
}

#' @rdname ecdf
#' @export
ecdf = function(x, ...){
    UseMethod("ecdf", x)
}

#' @export
ecdf.default = function(x, ...){
    stats::ecdf(x, ...)
}

#' @rdname ecdf
#' @title ecdf plot for kst enrichment test
#'
#' @description makes a ecdf plot for enrichment test result returned by
#' \code{\link{enrichment_test}} with test equals to kst.
#'
#' @param x EnrichmentKST. Returned by \code{\link{enrichment_test}}
#' when test equals to kst.
#' @param level character. The level to plot the ecdf plot. Must be in the
#' group.
#' @param alternative character. Either greater or less
#' @param colors character
#' @param ... other arguments not supported
#'
#' @return a ggplot2 object
#' @examples
#' lpd = transform_by_sample(lipidome, function(x) log(x/sum(x)))
#' design = model.matrix(~Treatment * Timepoint + Subject, data = lpd$pdata)
#' fit = model_fit(lpd, design, "TreatmentMed:TimepointPre", "limma")
#' en = enrichment_test(lpd, fit, "class", "kst")
#' ecdf(en, level = "PC", alternative = "greater")
#'
#' en = enrichment_test(lpd, fit, "class", "kst", alternative = "greater")
#' ecdf(en, level = "PC")
#'
#' en = enrichment_test(lpd, fit, "class", "kst", alternative = "two.sided")
#' ecdf(en, level = "PC")
#'
#' @export
#' @import ggplot2
#' @import dplyr
#' @author Chenghao Zhu
#' @seealso \code{\link{enrichment_test}}
ecdf.EnrichmentKST = function(x, level, alternative = c("greater", "less"),
                                  colors = c("blue", "red"), ...){
    stopifnot(level %in% x$group)
    stopifnot(is.character(colors))
    stopifnot(length(colors) <= 2)
    if(length(colors) == 1)
        colors = rep(colors, 2)

    alt = switch(
        x$alternative,
        "both" = match.arg(alternative, c("greater", "less")),
        "two.sided" = "raw",
        "adjusted"
    )

    df = data.frame(
        pvalue = x$fit_pvalues[,alt]
    )
    df$var = rownames(x$fit_pvalues)
    df %>%
        mutate(group = x$group) %>%
        filter(group == level) %>%
        ggplot() +
        stat_ecdf(geom = "step", aes(pvalue), color = colors[1]) +
        stat_ecdf(geom = "point", aes(pvalue), color = colors[1]) +
        stat_ecdf(geom = "step", aes(seq(0,1,length.out = length(pvalue))),
                  color = colors[2]) +
        stat_ecdf(geom = "point", aes(seq(0,1,length.out = length(pvalue))),
                  color = colors[2]) +
        labs(y = "Fn(pvalue)", title = level) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
}
