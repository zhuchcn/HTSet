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

test_that("HTSet class is initiated successfully", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    lpd = HTSet(edata = edata, fdata = fdata, pdata = pdata)
    expect_s4_class(lpd, "HTSet")
    expect_true(is.matrix(lpd@edata))
    expect_s3_class(lpd@pdata, "data.frame")
    expect_s3_class(lpd@fdata, "data.frame")
    expect_null(lpd@assay)
    expect_equal(nrow(lpd@edata), nrow(lpd@fdata))
    expect_equal(ncol(lpd@edata), nrow(lpd@pdata))

    # when pdata is not given
    lpd = HTSet(edata = edata, fdata = fdata)
    expect_s4_class(lpd, "HTSet")
    expect_true(is.null(lpd@pdata))

    # when fdata is not given
    lpd = HTSet(edata = edata, pdata = pdata)
    expect_s4_class(lpd, "HTSet")
    expect_true(is.null(lpd@fdata))
})

test_that("HTSet will fail", {
    edata = lipidome$edata
    fdata = lipidome$fdata
    pdata = lipidome$pdata
    pdata2 = pdata[1:20,]
    expect_error(HTSet(edata, fdata, pdata2))
    pdata2 = pdata
    rownames(pdata2) = NULL
    expect_error(HTSet(edata, fdata, pdata2))

    fdata2 = fdata[1:20,]
    expect_error(HTSet(edata, fdata2, pdata))
    fdata2 = fdata
    rownames(fdata2) = NULL
    expect_error(HTSet(edata, fdata2, pdata))

    edata2 = as.character(edata)
    expect_error(HTSet(edata2, fdata, pdata))
})

