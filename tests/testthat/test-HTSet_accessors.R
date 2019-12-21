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

test_that("HTSet accessor $", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata = edata, fdata = fdata, pdata = pdata)
    edata2 = lpd$edata
    expect_equal(class(edata2), class(edata))
    expect_equal(nrow(edata2), nrow(edata))
    expect_equal(ncol(edata2), ncol(edata))
    pdata2 = lpd$pdata
    expect_equal(class(pdata2), class(pdata))
    expect_equal(nrow(pdata2), nrow(pdata))
    expect_equal(ncol(pdata2), ncol(pdata))
    fdata2 = lpd$fdata
    expect_equal(class(fdata2), class(fdata))
    expect_equal(nrow(fdata2), nrow(fdata))
    expect_equal(ncol(fdata2), ncol(fdata))

    expect_error(lpd$hahaha)

    edata2 = apply(edata2, 2, function(x) x/sum(x))
    lpd$edata = edata2
    expect_true(all.equal(rowSums(lpd$edata), rowSums(edata2)))
    pdata2$fake_var = "fake"
    lpd$pdata = pdata2
    expect_true(all(lpd$pdata$fake_var == "fake"))
    lpd$pdata$fake_var = "fake2"
    expect_true(all(lpd$pdata$fake_var == "fake2"))
    lpd$pdata$fake_var2 = "fake3"
    expect_true(all(lpd$pdata$fake_var2 == "fake3"))
    fdata2$fake_var = "fake"
    lpd$fdata = fdata2
    expect_true(all(lpd$fdata$fake_var == "fake"))
    lpd$fdata$fake_var = "fake2"
    expect_true(all(lpd$fdata$fake_var == "fake2"))
    lpd$fdata$fake_var2 = "fake3"
    expect_true(all(lpd$fdata$fake_var2 == "fake3"))
})

test_that("HTSet accessor [[", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata = edata, fdata = fdata, pdata = pdata)
    edata2 = lpd[["edata"]]
    expect_equal(class(edata2), class(edata))
    expect_equal(nrow(edata2), nrow(edata))
    expect_equal(ncol(edata2), ncol(edata))
    pdata2 = lpd[["pdata"]]
    expect_equal(class(pdata2), class(pdata))
    expect_equal(nrow(pdata2), nrow(pdata))
    expect_equal(ncol(pdata2), ncol(pdata))
    fdata2 = lpd[["fdata"]]
    expect_equal(class(fdata2), class(fdata))
    expect_equal(nrow(fdata2), nrow(fdata))
    expect_equal(ncol(fdata2), ncol(fdata))

    expect_error(lpd[['hahaha']])

    edata2 = apply(edata2, 2, function(x) x/sum(x))
    lpd[["edata"]] = edata2
    expect_true(all.equal(rowSums(lpd@edata), rowSums(edata2)))
    pdata2$fake_var = "fake"
    lpd[["pdata"]] = pdata2
    expect_true(all(lpd@pdata$fake_var == "fake"))
    lpd[["pdata"]]$fake_var = "fake2"
    expect_true(all(lpd[["pdata"]]$fake_var == "fake2"))
    lpd[["pdata"]]$fake_var2 = "fake3"
    expect_true(all(lpd[["pdata"]]$fake_var2 == "fake3"))
    fdata2$fake_var = "fake"
    lpd[["fdata"]] = fdata2
    expect_true(all(lpd[["fdata"]]$fake_var == "fake"))
    lpd[["fdata"]]$fake_var = "fake2"
    expect_true(all(lpd[["fdata"]]$fake_var == "fake2"))
    lpd[["fdata"]]$fake_var2 = "fake3"
    expect_true(all(lpd[["fdata"]]$fake_var2 == "fake3"))
})

test_that("HTSet accessor [", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata = edata, fdata = fdata, pdata = pdata)
    lpd2 = lpd[1:5,]
    expect_equal(nfeatures(lpd2), 5)
    lpd2 = lpd[,1:6]
    expect_equal(nsamples(lpd2), 6)
    lpd2 = lpd[1:7,1:7]
    expect_equal(nsamples(lpd2), 7)
    expect_equal(nfeatures(lpd2), 7)
    lpd2 = lpd[,,1:3]
    expect_equal(ncol(lpd2@fdata), 3)
    lpd2 = lpd[,,,1:4]
    expect_equal(ncol(lpd2@pdata), 4)
})
