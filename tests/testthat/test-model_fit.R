test_that("model_fit failure", {
    design = model.matrix(~Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    expect_error(model_fit(mtcars, design, coef, "limma"),
                 regexp = "*HTSet objects*")
    expect_error(model_fit(exrna, design, c("A", "B"), "limma"),
                 regexp = "*one coef*")
    expect_error(model_fit(exrna, design, "condition", "limma"),
                 regexp = "*not found*")
    expect_error(model_fit(exrna, design[1:3,], coef, "limma"),
                 regexp = "*design matrix invalid*")
    expect_error(model_fit(exrna, design, coef, c("limma", "DESeq2")),
                 regexp = "engine number > 1")
    expect_error(model_fit(exrna, design, coef, "limma", adjust.method = "BOOP"),
                 regexp = "adjust.method not valid")
    expect_error(model_fit(exrna, design, coef, "edgeR", edger_model = "HAHA"),
                 regexp = "edger_model.*not valid")
    expect_error(model_fit(exrna, design, coef, "super engine"),
                 regexp = "The engines supported are.*")
})

test_that("model_fit with limma", {
    # test model_fit with voom
    design = model.matrix(~Condition, data = exrna$pdata)
    res = model_fit(exrna, design, "ConditionSystemic Lupus Erythematosus",
                    "limma", voom = TRUE)
    expect_s3_class(res, "ModelFit")
    expect_equal(res$engine, "limma")
    expect_true(res$params$voom)
    expect_equal(nrow(res$results), nfeatures(exrna))
    expect_equal(nrow(res$design), nsamples(exrna))
    expect_true(all.equal(rownames(res$results), featureNames(exrna)))
    expect_true(all.equal(rownames(res$design), sampleNames(exrna)))
    expect_equal(res$coef, "ConditionSystemic Lupus Erythematosus")
    expect_equal(res$adjust.method, "BH")

    # test model_fit with voom = FALSE
    design = model.matrix(~Treatment * Timepoint + Subject,
                          data = lipidome$pdata)
    res = model_fit(lipidome, design, "TreatmentMed:TimepointPre", "limma")
    expect_s3_class(res, "ModelFit")
    expect_equal(res$engine, "limma")
    expect_false(res$params$voom)
    expect_equal(nrow(res$results), nfeatures(lipidome))
    expect_equal(nrow(res$design), nsamples(lipidome))
    expect_true(all.equal(rownames(res$results), featureNames(lipidome)))
    expect_true(all.equal(rownames(res$design), sampleNames(lipidome)))
    expect_equal(res$coef, "TreatmentMed:TimepointPre")
    expect_equal(res$adjust.method, "BH")

    res = model_fit(lipidome, design, "TreatmentMed:TimepointPre", "limma",
                    transform = log)
    expect_s3_class(res, "ModelFit")
})

test_that("model_fit with edgeR", {
    # qlf
    design = model.matrix(~Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    res = model_fit(exrna, design, coef, "edgeR")
    expect_s3_class(res, "ModelFit")
    expect_equal(res$engine, "edgeR")
    expect_equal(res$params$edger_model, "qlf")
    expect_equal(nrow(res$results), nfeatures(exrna))
    expect_equal(nrow(res$design), nsamples(exrna))
    expect_true(all.equal(rownames(res$results), featureNames(exrna)))
    expect_true(all.equal(rownames(res$design), sampleNames(exrna)))
    expect_equal(res$coef, coef)
    expect_equal(res$adjust.method, "BH")

    # lrt
    design = model.matrix(~Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    res = model_fit(exrna, design, coef, "edgeR", edger_model = "lrt")
    expect_s3_class(res, "ModelFit")
    expect_equal(res$engine, "edgeR")
    expect_equal(res$params$edger_model, "lrt")
    expect_equal(nrow(res$results), nfeatures(exrna))
    expect_equal(nrow(res$design), nsamples(exrna))
    expect_true(all.equal(rownames(res$results), featureNames(exrna)))
    expect_true(all.equal(rownames(res$design), sampleNames(exrna)))
    expect_equal(res$coef, coef)
    expect_equal(res$adjust.method, "BH")
})

test_that("model_fit with DESeq2", {
    design = model.matrix(~Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    res = model_fit(exrna, design, coef, "DESeq2")
    expect_s3_class(res, "ModelFit")
    expect_equal(res$engine, "DESeq2")
    expect_equal(length(res$params), 0)
    expect_equal(nrow(res$results), nfeatures(exrna))
    expect_equal(nrow(res$design), nsamples(exrna))
    expect_true(all.equal(rownames(res$results), featureNames(exrna)))
    expect_true(all.equal(rownames(res$design), sampleNames(exrna)))
    expect_equal(res$coef, coef)
    expect_equal(res$adjust.method, "BH")
})

test_that("volcano plot", {
    defaults = formals(volcanoplot.ModelFit)
    design = model.matrix(~Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    res = model_fit(exrna, design, coef, "limma", voom = TRUE)
    g = volcanoplot(res)
    expect_is(g, "ggplot")
    expect_is(g$layers[[length(g$layers)]]$geom, "GeomPoint")
    expect_equal(g$layers[[length(g$layers)]]$aes_params$colour, defaults$color)
    expect_equal(g$layers[[length(g$layers)]]$aes_params$alpha, defaults$alpha)
    expect_equal(g$labels$title, coef)

    expect_error(volcanoplot(mtcars))
})

test_that("p value hist", {
    defaults = formals(hist.ModelFit)
    design = model.matrix(~Condition, data = exrna$pdata)
    coef = "ConditionSystemic Lupus Erythematosus"
    res = model_fit(exrna, design, coef, "limma", voom = TRUE)
    g = hist(res)
    expect_is(g, "ggplot")
    expect_is(g$layers[[1]]$geom, "GeomBar")
    expect_equal(g$layers[[1]]$aes_params$colour, defaults$color)
    expect_equal(g$layers[[1]]$aes_params$fill, defaults$fill)
})
