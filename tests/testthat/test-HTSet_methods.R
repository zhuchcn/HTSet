load_rawdata_lipidome = function(){
    filepaths = list.files(
        system.file("testdata", "lipidome", package = "HTSet"),
        pattern = "*data.csv", full.names = TRUE
    )
    filenames = sapply(filepaths, function(x){
        name = basename(x)
        gsub(".csv$", "", name)
    })
    res = lapply(filepaths, function(x) read.csv(x, row.names = 1))
    names(res) = filenames
    res$edata = as.matrix(res$edata)
    return(res)
}
lipidome = load_rawdata_lipidome()

test_that("HTSet dim, nsamples, and nfeatures", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata, fdata, pdata)
    expect_true(all.equal(dim(lpd), dim(edata)))
    expect_equal(nsamples(lpd), nrow(pdata))
    expect_equal(nfeatures(lpd), nrow(fdata))
})

test_that("HTSet sampleNames and featureNames", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata, fdata, pdata)
    expect_true(all.equal(sampleNames(lpd), colnames(edata)))
    expect_true(all.equal(featureNames(lpd), rownames(edata)))

    new_sn = paste("test", 1:ncol(edata))
    sampleNames(lpd) = new_sn
    expect_true(all.equal(sampleNames(lpd), new_sn))

    new_fn = paste("test", 1:nrow(edata))
    featureNames(lpd) = new_fn
    expect_true(all.equal(featureNames(lpd), new_fn))
})

test_that("HTSet subset and summarize", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata, fdata, pdata)
    lpd_sum = summarize_feature(lpd, "class")
    expect_equal(nfeatures(lpd_sum), length(unique(fdata$class)))

    lpd_sub = subset_samples(lpd, 1:5)
    expect_equal(nsamples(lpd_sub), 5)
    lpd_sub = subset_samples(lpd, pdata$Treatment == "FF")
    expect_equal(nsamples(lpd_sub), 20)
    lpd_sub = subset_samples(lpd, rownames(pdata)[1:10])
    expect_equal(nsamples(lpd_sub), 10)

    expect_error(subset_samples(lpd, 1:100))
    expect_error(subset_samples(lpd, rep(TRE, 10)))
    expect_error(subset_samples(lpd, c(rownames(pdata)[1:10], "ahaha")))

    lpd_sub = subset_features(lpd, 1:5)
    expect_equal(nfeatures(lpd_sub), 5)
    lpd_sub = subset_features(lpd, fdata$class == "PC")
    expect_equal(nfeatures(lpd_sub), sum(fdata$class == "PC"))
    lpd_sub = subset_features(lpd, rownames(fdata)[1:10])
    expect_equal(nfeatures(lpd_sub), 10)

    expect_error(subset_features(lpd, 1:1000))
    expect_error(subset_features(lpd, rep(TRE, 10)))
    expect_error(subset_features(lpd, c(rownames(fdata)[1:10], "ahaha")))
})

test_that("Test transform functions", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata, fdata, pdata)

    lpd2 = transform_by_sample(lpd, function(x) x/sum(x, na.rm = TRUE))
    expect_true(sum(colSums(lpd2@edata)) == nsamples(lpd))

    lpd2 = transform_by_feature(lpd, scale, center = TRUE, scale = TRUE)
    expect_true(as.integer(sum(rowMeans(lpd2@edata))) == 0)
})
