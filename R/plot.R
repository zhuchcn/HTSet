#' @title make a boxplot
#' @description Plotting a boxplot of a specified feature
#' @param object \code{\link{HTSet-class}}
#' @param feature character. The feature to use as y axis.
#' @param x character. The sample variable to be used as x axis.
#' @param cols character. The sample variable to split vertically.
#' @param rows character. The sample variable to split horizontally.
#' @param line_by characte. The sample variable to draw lines.
#' @param fill_by character. The sample variable to fill colors to boxes.
#' @param color_by character. The sample variable to color points and/or points.
#' @param points boolean. Whether to draw points on top of boxes.
#' @param point_color character. The points color. Ignored if color_by is
#' specified.
#' @param line_color character. The lines color. Ignored if color_by is
#' specified.
#' @param box_width numeric. The box width. Default: 0.75
#' @param whiskers boolean. Whether to draw wiskers.
#' @param whisker_width numeric. The width of wiskers.
#'
#' @return ggplot
#'
#' @export
#' @author Chenghao Zhu
#' @import ggplot2 dplyr
plot_boxplot = function(
    object, feature, x, rows = NULL, cols = NULL, line_by = NULL,
    fill_by = NULL, color_by = NULL, points = FALSE, point_color = "black",
    line_color = "black", box_width = 0.75, whiskers = TRUE,
    whisker_width = 0.2
) {
    stopifnot(is(object, "HTSet"))
    stopifnot(feature %in% featureNames(object))

    sample_vars = c(x, rows, cols, line_by, fill_by, color_by)
    sample_vars = unique(sample_vars)
    sample_vars = sample_vars[!is.na(sample_vars)]

    df = data.frame(
        feature = object@edata[feature,]
    ) %>%
        cbind(object@pdata[, sample_vars])

    # main plot
    p = ggplot(df, aes_string(x = x, y = "feature"))

    # boxplot fill color
    if(is.null(fill_by)) {
        p = p + geom_boxplot(width = box_width)
    } else {
        p = p + geom_boxplot(aes_string(fill = fill_by), width = box_width)
    }

    # whiskers
    if(whiskers) {
        p = p + stat_boxplot(geom = "errorbar", width = whisker_width)
    }

    # points
    if(points){
        if(is.null(color_by)) {
            p = p + geom_point(color = point_color)
        } else {
            p = p + geom_point(aes(color = color_by))
        }
    }

    # lines
    if(!is.null(line_by)) {
        if(is.null(color_by)) {
            p = p + geom_line(aes_string(group = line_by), color = line_color)
        } else {
            p = p + geom_line(aes_string(group = line_by, color = color_by))
        }
    }

    # facets
    if(!is.null(cols) & !is.null(rows)) {
        p = p + facet_grid(rows = vars(!!sym(rows)), cols = vars(!!sym(cols)))
    } else if(!is.null(cols)) {
        p = p + facet_grid(cols = vars(!!sym(cols)))
    } else if(!is.null(rows)) {
        p = p + facet_grid(rows = vars(!!sym(rows)))
    }

    p = p + theme_bw()
    return(p)
}
