test_that("enrichment test argument validation", {
    lpd = transform_by_sample(lipidome, function(x) log(x/sum(x)))
    lpd_class = summarize_feature(lpd, "class")
    design = model.matrix(~ Treatment * Timepoint + Subject + 1, data = lpd$pdata)
    fit = model_fit(lpd, design, coef = "TreatmentMed:TimepointPre", engine = "limma")
    expect_error(enrichment_test(lpd_class, fit, test = "fet"),
                 "object and fit must match")
    expect_error(enrichment_test(lpd, fit, test = c("fet", "ks")))
    expect_error(enrichment_test(lpd, fit, test = "fet", p.cutoff = 1:5),
                 "invalid p.cutoff")
    expect_error(enrichment_test(lpd, fit, test = "fet", p.cutoff = 5),
                 "invalid p.cutoff")
    expect_error(enrichment_test(lpd, fit, test = "two.sided", p.cutoff = 1))
})

test_that("fet enrichment", {
    design = model.matrix(~ Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    fit = model_fit(exrna, design, coef, engine = "limma", args = list(voom = TRUE))
    en = enrichment_test(exrna, fit, "gene_type", "fet")
    expect_is(en, "EnrichmentFET")
    expect_equal(en$alternative, "both")

    en = enrichment_test(exrna, fit, "gene_type", "fet", "greater")
    expect_equal(en$alternative, "greater")
    expect_is(en$pval, "numeric")

    en = enrichment_test(exrna, fit, "gene_type", "fet", "less")
    expect_equal(en$alternative, "less")
    expect_is(en$pval, "numeric")

    en = enrichment_test(exrna, fit, "gene_type", "fet", "two.sided")
    expect_equal(en$p.cutoff, 0.05)
    expect_equal(en$alternative, "two.sided")

    fit = model_fit(exrna, design, coef, engine = "edgeR")
    en = enrichment_test(exrna, fit, "gene_type", "fet")
    expect_is(en, "EnrichmentFET")

    fit = model_fit(exrna, design, coef, engine = "edgeR", args = list(model = "lrt"))
    en = enrichment_test(exrna, fit, "gene_type", "fet")
    expect_is(en, "EnrichmentFET")

    fit = model_fit(exrna, design, coef, engine = "DESeq2")
    en = enrichment_test(exrna, fit, "gene_type", "fet")
    expect_is(en, "EnrichmentFET")

    fit = model_fit(exrna, design, coef, engine = "DESeq2", args = list(DESeq = list(useT = TRUE)))
    en = enrichment_test(exrna, fit, "gene_type", "fet")
    expect_is(en, "EnrichmentFET")

    reduced = model.matrix(~1, data = exrna$pdata)
    fit = model_fit(exrna, design, coef, engine = "DESeq2", args = list(DESeq = list(test = "LRT", reduced = reduced)))
    en = enrichment_test(exrna, fit, "gene_type", "fet")
    expect_is(en, "EnrichmentFET")
})

test_that("kst enrichment", {
    design = model.matrix(~ Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    fit = model_fit(exrna, design, coef, engine = "limma", args = list(voom = TRUE))
    expect_warning(en <- enrichment_test(exrna, fit, "gene_type", "kst"))
    expect_is(en, "EnrichmentKST")
    expect_equal(en$alternative, "both")

    expect_warning(en <- enrichment_test(exrna, fit, "gene_type", "kst", "greater"))
    expect_equal(en$alternative, "greater")
    expect_is(en$pval, "numeric")

    expect_warning(en <- enrichment_test(exrna, fit, "gene_type", "kst", "less"))
    expect_equal(en$alternative, "less")
    expect_is(en$pval, "numeric")

    fit = model_fit(exrna, design, coef, engine = "edgeR")
    expect_warning(en <- enrichment_test(exrna, fit, "gene_type", "kst"))
    expect_is(en, "EnrichmentKST")

    fit = model_fit(exrna, design, coef, engine = "edgeR", args = list(model = "lrt"))
    expect_warning(en <- enrichment_test(exrna, fit, "gene_type", "kst"))
    expect_is(en, "EnrichmentKST")

    fit = model_fit(exrna, design, coef, engine = "DESeq2")
    expect_warning(en <- enrichment_test(exrna, fit, "gene_type", "kst"))
    expect_is(en, "EnrichmentKST")

    fit = model_fit(exrna, design, coef, engine = "DESeq2", args = list(DESeq = list(useT = TRUE)))
    expect_warning(en <- enrichment_test(exrna, fit, "gene_type", "kst"))
    expect_is(en, "EnrichmentKST")

    reduced = model.matrix(~1, data = exrna$pdata)
    fit = model_fit(exrna, design, coef, engine = "DESeq2", args = list(DESeq = list(test = "LRT", reduced = reduced)))
    expect_error(enrichment_test(exrna, fit, "gene_type", "kst"),
                 "DESeq2's LRT test is not supported for ks test yet.")
})

test_that("enrichment fet barplot", {
    lpd = transform_by_sample(lipidome, function(x) log(x/sum(x)))
    design = model.matrix(~Treatment * Timepoint + Subject, data = lpd$pdata)
    fit = model_fit(lpd, design, "TreatmentMed:TimepointPre", "limma")
    en = enrichment_test(lpd, fit, "class", "fet")
    p = barplot(en, each = 3)
    expect_is(p, "ggplot")
})

test_that("enrichment kst ecdf", {
    lpd = transform_by_sample(lipidome, function(x) log(x/sum(x)))
    design = model.matrix(~Treatment * Timepoint + Subject, data = lpd$pdata)
    fit = model_fit(lpd, design, "TreatmentMed:TimepointPre", "limma")
    en = enrichment_test(lpd, fit, "class", "kst")
    p = ecdf(en, level = "PC", alternative = "greater")
    expect_is(p, "ggplot")
})
