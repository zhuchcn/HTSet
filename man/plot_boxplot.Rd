% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_boxplot}
\alias{plot_boxplot}
\title{make a boxplot}
\usage{
plot_boxplot(
  object,
  feature,
  x,
  rows = NULL,
  cols = NULL,
  line_by = NULL,
  fill_by = NULL,
  color_by = NULL,
  points = FALSE,
  point_color = "black",
  line_color = "black",
  box_width = 0.75,
  whiskers = TRUE,
  whisker_width = 0.2
)
}
\arguments{
\item{object}{\code{\link{HTSet-class}}}

\item{feature}{character. The feature to use as y axis.}

\item{x}{character. The sample variable to be used as x axis.}

\item{rows}{character. The sample variable to split horizontally.}

\item{cols}{character. The sample variable to split vertically.}

\item{line_by}{characte. The sample variable to draw lines.}

\item{fill_by}{character. The sample variable to fill colors to boxes.}

\item{color_by}{character. The sample variable to color points and/or points.}

\item{points}{boolean. Whether to draw points on top of boxes.}

\item{point_color}{character. The points color. Ignored if color_by is
specified.}

\item{line_color}{character. The lines color. Ignored if color_by is
specified.}

\item{box_width}{numeric. The box width. Default: 0.75}

\item{whiskers}{boolean. Whether to draw wiskers.}

\item{whisker_width}{numeric. The width of wiskers.}
}
\value{
ggplot
}
\description{
Plotting a boxplot of a specified feature
}
\author{
Chenghao Zhu
}
