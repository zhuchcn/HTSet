#' @title format lipid names
#' @description Format lipid names into a unified format, such as PC 36:1. This
#' is useful for downstream processing such as calculation $EOD_18$ or ACL.
#' @param x character.
#' @return character
#' @import stringr
#' @export
#' @author Chenghao Zhu
format_lipid_name = function(x){
    x = gsub("^(.+)or.+$","\\1",x)
    x = gsub("^(.+); .+", "\\1", x)
    x = gsub("\\s+$", "", x)
    x = gsub(";$", "", x)
    # if fatty acids are listed separated
    pattern = "\\d{1,2}:\\d{1}\\/{1}\\d{1,2}:\\d{1}"
    for(i in 1:length(x)){
        if(grepl(pattern, x[i])){
            x[i] = gsub(
                pattern,
                str_extract(x[i],pattern) %>%
                    str_split("\\/", simplify = T) %>%
                    str_split("\\:", simplify = T) %>%
                    apply(2, function(xx) sum(as.integer(xx))) %>%
                    paste(collapse =":"),
                x[i])
        }
    }

    # remove Z, H, and E
    x = gsub("\\(\\d{1,2}[ZHE]\\)","",x)
    # Cholesterol
    x = gsub("^cholesterol$", "Cholesterol", x)
    # CE
    x = gsub("^.*?(\\d{1,2}\\:\\d{1}).*?(Cholesteryl ester|Cholesteryl Ester|CHOLESTERYL ESTER|cholesteryl ester).*$", "CE \\1", x)
    x = gsub("^.*?(\\d{1,2}\\:\\d{1}).*?(CE|ce)\\W*$", "CE \\1", x)
    x = gsub("^.*?(Cholesteryl ester|cholesteryl ester|Cholesteryl Ester|CHOLESTERYL ESTER).*?(\\d{1,2}\\:\\d{1}).*$", "CE \\2", x)
    x = gsub("^(CE|ce)\\W*?(\\d{2}\\:\\d{1}).*$", "CE \\2", x)
    # DG
    x = gsub('DAG', 'DG', x)
    x = gsub('Diacylglycerol', 'DG', x)
    x = gsub('diacylglycerol', 'DG', x)
    x = gsub('DIACYLGLYCEROL', 'DG', x)
    x = gsub("^\\W*DG\\W*?(\\d{1,2}\\:\\d{1}).*$", "DG \\1", x)
    # MG
    x = gsub('MAG', "MG", x)
    x = gsub("Monoacylglycerol", "MG", x)
    x = gsub("monoacylglycerol", "MG", x)
    x = gsub("MONOACYLGLYCEROL", "MG", x)
    x = sub("^\\W*MG\\W*?(\\d{1,2}\\:\\d{1}).*$", "MG \\1", x)
    # FA
    x = gsub("^FA\\W*?(\\d{1,2}\\:\\d{1}).*$", "FA \\1", x)
    # DCER
    x = gsub('^DCER\\W*?([0-9]{1,2}:[0-9]{1}).*$', 'Dihydroceramide \\1', x)
    # HCER
    x = gsub('^HCER\\W*?([0-9]{1,2}:[0-9]{1}).*$', 'Hexosylceramide \\1', x)
    # Gal-Gal-Cer
    x = gsub("LCER", 'Lactosylceramide', x)
    x = gsub("^Gal-Gal-Cer\\W*?d(\\d{1,2}:\\d{1}).*$", "Gal-Gal-Cer \\1 d",x)
    x = gsub("^Lactosylceramide\\W*?d*(\\d{1,2}:\\d{1}).*$", "Gal-Gal-Cer \\1 d",x)
    # GlcCer
    x = gsub("^GlcCer\\W*?d(\\d{1,2}:\\d{1}).*$", "GlcCer \\1 d", x)
    # Ceramide
    x = gsub("CERAMIDE", 'Ceramide', x)
    x = gsub('CER', 'Ceramide', x)
    x = gsub("^Ceramide\\W*?d*(\\d{1,2}:\\d{1}).*$", "Ceramide \\1 d",x)
    x = gsub("^Cer\\W*?d(\\d{2}\\:\\d{1}).*$", "Ceramide \\1 d",x)
    # phosphoglyceropipids
    x = gsub("[Ll]yso(PCEASI)", "L\\1", x)
    x = gsub("^(L{0,1}P[CEASI])\\W*?(\\d{1,2}:\\d{1}).*$", "\\1 \\2", x)
    x = gsub("^(L{0,1}P[CEASI])\\W*?([opOP])-(\\d{1,2}:\\d{1}).*$", "\\1 \\3 \\2", x)
    x = gsub("^Plasmenyl-(L{0,1}P[CEASI])\\W*?(\\d{1,2}:\\d{1}).*$", "\\1 \\2 p", x)
    x = gsub("Ox(L{0,1}P[CEASI])\\W*?(\\d{2}:\\d{1}).*$", "\\1 \\2 ox", x)
    # SM
    x = gsub("^SM\\W*?d*(\\d{1,2}:\\d{1}).*$","SM \\1 d", x)
    # TG
    x = gsub("TAG", "TG", x)
    x = gsub("^TG\\W*?(\\d{1,2}:\\d{1}).*$","TG \\1", x)
    return(x)
}

#' @title Extract the carbon chain length
#' @description Extract the total number of carbons from a lipid annotation
#' name. For example, the carbon chain length of lipid PC 34:2 is 34. The lipid
#' names must be first cleaned by \code{\link{format_lipid_name}}.
#' @param x character. The lipid annotation name. Must be first cleaned by
#' \code{\link{format_lipid_name}}.
#' @return integer
#' @export
#' @author Chenghao Zhu
nCarbons = function(x){
    ifelse(
        grepl("Cholesterol", x), 0,
        as.integer(str_split_fixed(str_split_fixed(x, " ", n = 2)[,2], ":", 2)[,1])
    )
}

#' @title number of fatty acyls
#' @description This function returns the number of fatty acyls of a given
#' lipid annotation name. For example, the number of fatty acyl of lipid PC 34:2
#' is 2. The lipid annotation names must be first cleaned by
#' \code{\link{format_lipid_name}}.
#' @param x character. The lipid annotation name. Must be first cleaned by
#' \code{\link{format_lipid_name}}.
#' @return integer
#' @export
#' @author Chenghao Zhu
nFattyAcyls = function(x){
    class = str_split_fixed(x, " ", 2)[,1]
    return(ifelse(
        class == 'Cholesterol', 0,
        ifelse(
            class == 'TG', 3,
            ifelse(
                class %in% c('LPC', "LPE", "LPG", "LPI", "LPA", "MG", "FA", "CE"),
                1, 2
            )
        )
    ))
}

#' @title number of double bonds
#' @description This function returns the number of double bonds of a given
#' lipid annotation name. For example, the number of double bonds of lipid PC 34:2
#' is 2. The lipid annotation names must be first cleaned by
#' \code{\link{format_lipid_name}}.
#' @param x character. The lipid annotation name. Must be first cleaned by
#' \code{\link{format_lipid_name}}.
#' @return integer
#' @export
#' @author Chenghao Zhu
nDoubleBonds = function(x){
    ifelse(
        grepl("Cholesterol", x), 0,
        as.integer(str_split_fixed(str_split_fixed(x, " ", 3)[,2], ":", 2)[,2])
    )
}

#' @title summarize average carbon chain length
#' @description Calculate the ACL (average carbon chain length) of each lipid class.
#' The ACL of lipid class k in sample i is calculated using the equation below:
#'
#' \ifelse{html}{\out{
#'   <center><i>
#'   ACL<sub>i,k</sub> = Σ (conc<sub>i,j,k</sub>*nc<sub>j,k</sub>)
#'   / Σ (conc<sub>i,j,k</sub>*nfa<sub>j,k</sub>)
#'   </i></center>
#' }}{\deqn{
#'   ACL_{i,k} = (\sum conc_{i,j,k} \times nc_{j,k}) / (\sum conc_{j,k} \times nfa_{j,k})
#' }}
#'
#' Where \ifelse{html}{\out{<i>conc<sub>i,j,k</sub></i>}}{\eqn{conc_{i,j,k}}}
#' is the mol concentration of lipid \emph{j} that belongs to class \emph{k} in
#' sample \emph{i}; \ifelse{html}{\out{<i>nc<sub>j,k</sub></i>}}{\eqn{nc_{j,k}}}
#' is the number of carbon of lipid \emph{j} that belongs to class \emph{k}; and
#' \ifelse{html}{\out{<i>nfa<sub>j,k</sub></i>}}{\eqn{nfa_{j,k}}} is the number
#' of fatty acyls of lipid \emph{j} that belongs to class \emph{k}.
#'
#' @param object \code{\link{HTSet-class}}
#' @param name character. The column name of the feature variable contains the
#' annotation name in fdata. The annotation name must be first cleaned by the
#' \code{\link{format_lipid_name}}.
#' @param class character. The name of the feature variable that contains the
#' lipid class.
#' @return \code{\link{HTSet-class}}
#' @references
#' Zhu, C. et al., The HDL lipidome is widely remodeled by fast food versus
#' Mediterranean diet in 4 days. \emph{Metabolomics}, 2019 Aug 17;15(8):114.
#' doi: 10.1007/s11306-019-1579-1.
#' @export
#' @author Chenghao Zhu
summarize_ACL = function(object, name, class){
    if(!is(object, "HTSet"))
        stop("object must be a HTSet object")
    if(!name %in% colnames(object@fdata))
        stop("name " %+% name %+% " not found")
    if(!class %in% colnames(object@fdata))
        stop("class " %+% class %+% " not found")

    nCB = nCarbons(object@fdata[, name])
    nFA = nFattyAcyls(object@fdata[, name])
    class = object@fdata[, class]

    mat1 = apply(object@edata * nCB, 2, function(col){
        tapply(col, class, function(x) sum(x, na.rm = TRUE))
    })
    mat2 = apply(object@edata * nFA, 2, function(col){
        tapply(col, class, function(x) sum(x, na.rm = TRUE))
    })
    mat = mat1 / mat2
    mat = mat[rownames(mat) != "Cholesterol",]
    Overall = apply(object@edata * nCB, 2, function(x) sum(x, na.rm = TRUE)) /
        apply(object@edata * nFA, 2, function(x) sum(x, na.rm = TRUE))
    mat = rbind(Overall, mat)
    HTSet(mat, pdata = object@pdata)
}

#' @title Equivalent Double Bond per 18 carbons
#' @description Calculate the \ifelse{html}{\out{EOD<sub>18</sub>}}{\eqn{EOD_{18}}}
#' (equivalen double bonds per 18 carbons) of each lipid class.
#'
#' \ifelse{html}{EOD<sub>18</sub>}{\eqn{EOD_{18}}} of class \emph{k} in sample
#' \emph{i} is calculated using the equation below:
#'
#' \ifelse{html}{\out{
#'   <center><i>
#'   EOD<sub>18 j,k</sub> = Σ (conc<sub>i,j,k</sub>*ndb<sub>j,k</sub>)
#'   / Σ (conc<sub>i,j,k</sub>*nc<sub>j,k</sub>) * 18
#'   </i></center>
#' }}{\deqn{
#'   EOD_{18 i,k} = (\sum conc_{i,j,k} \times ndb_{j,k}) / (\sum conc_{i,j,k} x nc_{j,8}) \times 18
#' }}
#'
#' Where \ifelse{html}{\out{<i>conc<sub>i,j,k</sub></i>}}{\eqn{conc_{i,j,k}}}
#' is the mol concentration of lipid \emph{j} that belongs to class \emph{k} in
#' sample \emph{i}; \ifelse{html}{\out{<i>ndb<sub>j,k</sub></i>}}{\eqn{nfa_{j,k}}}
#' is the number of double bonds of lipid \emph{j} that belongs to class
#' \emph{k}; and \ifelse{html}{\out{<i>nc<sub>j,k</sub></i>}}{\eqn{nc_{j,k}}}
#' is the number of carbon of lipid \emph{j} that belongs to class \emph{k}.
#'
#' @param object \code{\link{HTSet-class}}
#' @param name character. The column name of the feature variable contains the
#' annotation name in fdata. The annotation name must be first cleaned by the
#' \code{\link{format_lipid_name}}.
#' @param class character. The name of the feature variable that contains the
#' lipid class.
#' @return \code{\link{HTSet-class}}
#' @references
#' Zhu, C. et al., The HDL lipidome is widely remodeled by fast food versus
#' Mediterranean diet in 4 days. \emph{Metabolomics}, 2019 Aug 17;15(8):114.
#' doi: 10.1007/s11306-019-1579-1.
#' @export
#' @author Chenghao Zhu
summarize_EOD = function(object, name, class){
    if(!is(object, "HTSet"))
        stop("object must be a HTSet object")
    if(!name %in% colnames(object@fdata))
        stop("name " %+% name %+% " not found")
    if(!class %in% colnames(object@fdata))
        stop("class " %+% class %+% " not found")

    nDB = nDoubleBonds(object@fdata[, name])
    nCB = nCarbons(object@fdata[, name])
    class =object@fdata[, class]

    mat1 = apply(object@edata * nDB, 2, function(col){
        tapply(col, class, function(x) sum(x, na.rm = TRUE))
    })
    mat2 = apply(object@edata * nCB, 2, function(col){
        tapply(col, class, function(x) sum(x, na.rm = TRUE))
    })
    mat = mat1 / mat2 * 18
    mat = mat[rownames(mat) != "Cholesterol",]
    Overall = apply(object@edata * nDB, 2, function(x) sum(x, na.rm = TRUE))/
        apply(object@edata * nCB, 2, function(x) sum(x, na.rm = TRUE)) * 18
    mat = rbind(Overall, mat)
    HTSet(mat, pdata = object$pdata)
}

#' @title summarize odd chain lipids
#' @description Summarize odd chain lipids abundance of each lipid class.
#' @param object \code{\link{HTSet-class}}
#' @param name character. The column name of the feature variable contains the
#' annotation name in fdata. The annotation name must be first cleaned by the
#' \code{\link{format_lipid_name}}.
#' @param class character. The name of the feature variable that contains the
#' lipid class.
#' @return \code{\link{HTSet-class}}
#' @export
#' @author Chenghao Zhu
summarize_odd_chain = function(object, name, class){
    if(!is(object, "HTSet"))
        stop("object must be a HTSet object")
    if(!name %in% colnames(object@fdata))
        stop("name " %+% name %+% " not found")
    if(!class %in% colnames(object@fdata))
        stop("class " %+% class %+% " not found")

    nCB = nCarbons(object@fdata[, name])

    object = subset_features(object, nCB %% 2 == 1)
    Overall = apply(object@edata, 2, function(x) sum(x, na.rm = TRUE))
    object = summarize_feature(object, class)

    edata = rbind(Overall, object@edata)

    HTSet(edata, pdata = object@pdata)
}

#' @title Calculate lipid classes ratios
#' @description Calculate ratios between lipid classes. Those calculated lipid
#' ratios are often indicators for enzyme activities or have diagnostic
#' significance. The ratios calculated are PC/LPC, PE/LPE, CE/free cholesterol,
#' TG/DG, SM / Cer, CE/TG, surface lipid / core lipid, PC / surface lipid, and
#' SM / surface lipid.
#' @param object \code{\link{HTSet-class}}
#' @param name character. The column name of the feature variable contains the
#' annotation name in fdata. The annotation name must be first cleaned by the
#' \code{\link{format_lipid_name}}.
#' @param class character. The name of the feature variable that contains the
#' lipid class.
#' @return \code{\link{HTSet-class}}
#' @export
#' @author Chenghao Zhu
summarize_lipid_ratios = function(object, name, class){
    if(!is(object, "HTSet"))
        stop("object must be a HTSet object")
    if(!name %in% colnames(object@fdata))
        stop("name " %+% name %+% " not found")
    if(!class %in% colnames(object@fdata))
        stop("class " %+% class %+% " not found")

    object = summarize_feature(object, class)
    featureNames(object) = tolower(featureNames(object))
    edata = object@edata

    keys = list(
        c('ce',  "cholesterol"),
        c("pc",  "lpc"),
        c("pe",  "lpe"),
        c("pg",  'lpg'),
        c("pi",  "lpi"),
        c("pa",  "lpa"),
        c("tg",  "dg"),
        c("tag", "dag"),
        c("sm",  "cer"),
        c("sm",  "ceremide"),
        c("tg",  "ce"),
        c("tag", "ce")
    )

    ratios = NULL
    for(key in keys){
        if( key[1] %in% rownames(edata) & key[2] %in% rownames(edata)){
            ratios = rbind(ratios, edata[key[1],]/edata[key[2],])
            rownames(ratios)[nrow(ratios)] = key[1] %+% "/" %+% key[2]
        }
    }

    surface_lipids = c(
        "pc", "pe", "ps", "pi", "pg", "pa", "lpc", "lpe", "lps", "lpi", "lpa",
        "dg", "dag", "sm", "cer", "ceramide", "mg", "mag", "cholesterol",
        "dcer", 'ffa', 'hcer', 'lcer'
    )
    core_lipids = c("ce", "tg", 'tag')
    pl = c("pc", "pe", "ps", "pi", "pg", "pa", "lpc", "lpe", "lps", "lpi",
           "lpa", "sm")

    surface = colSums(edata[rownames(edata) %in% surface_lipids,])
    core = colSums(edata[rownames(edata) %in% core_lipids,])
    pl = colSums(edata[rownames(edata) %in% pl,])
    ratios = rbind(edata["pc",]/pl, ratios)
    ratios = rbind(edata["sm",]/pl, ratios)
    ratios = rbind(edata["pc",]/surface, ratios)
    ratios = rbind(edata["sm",]/surface, ratios)
    ratios = rbind(surface / core, ratios)
    rownames(ratios)[1:5] = c("surface/core", "sm/surface", "pc/surface", "sm/pl", "pc/pl")

    HTSet(ratios, pdata = object@pdata)
}

#' @title m/z to molecular weight
#' @description Calculate the molecular weight from m/z and adduct ion species.
#' The calculation is based on the adduct table extracted from the West Coast
#' Metabolomics Center's websit. See \code{\link{wcmc_adduct}}.
#' @param species character value of the adduct ion species
#' @param mz numeric value of the m/z
#' @return numeric
#' @author Chenghao Zhu
#' @seealso \code{\link{wcmc_adduct}}
#' @export
mz2molwt = function(species, mz){
    wcmc_adduct = wcmc_adduct()
    multiply = wcmc_adduct[species,]$Mult
    plus = wcmc_adduct[species,]$Mass
    (mz - plus) / multiply
}

#' @title show wcmc's adduct table
#' @description Show West Coast Metabolomics Center's ion adduct table.
#' @export
#' @import utils
#' @author Chenghao Zhu
#' @seealso \code{\link{mz2molwt}}
wcmc_adduct = function(){
    file = system.file("extdata", "wcmc-adduct.tsv", package = "HTSet")
    adduct = read.table(file, row.names = 1, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
    adduct$Mult = as.numeric(adduct$Mult)
    adduct$Mass = as.numeric(adduct$Mass)
    return(adduct)
}
