coord_radar2 <- function(theta = "x", start = 0, direction = 1, clip = "off") {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") {
    "y"
  } else {
    "x"
  }
  ggproto("CoordRadar", ggplot2::CoordPolar,
          theta = theta,
          r = r, start = start, clip = clip,
          direction = sign(direction), is_linear = function(coord) TRUE
  )
}

ggRadar2 <- function(data, mapping = NULL, rescale = TRUE, legend.position = "top",
                     colour = "red", alpha = 0.3, size = 3, ylim = NULL, scales = "fixed",
                     use.label = FALSE, interactive = FALSE, clip = "off", ...) {
  data <- as.data.frame(data)
  (groupname <- setdiff(names(mapping), c("x", "y")))
  groupname
  mapping
  length(groupname)
  if (length(groupname) == 0) {
    groupvar <- NULL
  }
  else {
    groupvar <- ggiraphExtra:::getMapping(mapping, groupname)
  }
  groupvar
  facetname <- colorname <- NULL
  if ("facet" %in% names(mapping)) {
    facetname <- ggiraphExtra:::getMapping(mapping, "facet")
  }
  (colorname <- setdiff(groupvar, facetname))
  if ((length(colorname) == 0) & !is.null(facetname)) {
    colorname <- facetname
  }
  data <- ggiraphExtra:::num2factorDf(data, groupvar)
  (select <- sapply(data, is.numeric))
  if ("x" %in% names(mapping)) {
    xvars <- ggiraphExtra:::getMapping(mapping, "x")
    xvars
    if (length(xvars) < 3) {
      warning("At least three variables are required")
    }
  }
  else {
    xvars <- colnames(data)[select]
  }
  (xvars <- setdiff(xvars, groupvar))
  if (rescale) {
    data <- ggiraphExtra:::rescale_df(data, groupvar)
  }
  temp <- sjlabelled::get_label(data)
  cols <- ifelse(temp == "", colnames(data), temp)
  if (is.null(groupvar)) {
    id <- ggiraphExtra:::newColName(data)
    data[[id]] <- 1
    longdf <- reshape2::melt(data, id.vars = id, measure.vars = xvars)
  }
  else {
    cols <- setdiff(cols, groupvar)
    longdf <- reshape2::melt(data, id.vars = groupvar, measure.vars = xvars)
  }
  temp <- paste0("plyr::ddply(longdf,c(groupvar,'variable'), dplyr::summarize,mean=mean(value,na.rm=TRUE))")
  df <- eval(parse(text = temp))
  colnames(df)[length(df)] <- "value"
  df
  groupvar
  if (is.null(groupvar)) {
    id2 <- ggiraphExtra:::newColName(df)
    df[[id2]] <- "all"
    id3 <- ggiraphExtra:::newColName(df)
    df[[id3]] <- 1:nrow(df)
    df$tooltip <- paste0(df$variable, "=", round(
      df$value,
      1
    ))
    df$tooltip2 <- paste0("all")
    p <- ggplot(data = df, aes_string(
      x = "variable", y = "value",
      group = 1
    )) +
      ggiraph::geom_polygon_interactive(aes_string(tooltip = "tooltip2"),
                                        colour = colour, fill = colour, alpha = alpha
      ) +
      ggiraph::geom_point_interactive(aes_string(
        data_id = id3,
        tooltip = "tooltip"
      ), colour = colour, size = size)
  }
  else {
    if (!is.null(colorname)) {
      id2 <- ggiraphExtra:::newColName(df)
      df[[id2]] <- df[[colorname]]
    }
    id3 <- ggiraphExtra:::newColName(df)
    df[[id3]] <- 1:nrow(df)
    df$tooltip <- paste0(
      groupvar, "=", df[[colorname]], "<br>",
      df$variable, "=", round(df$value, 1)
    )
    df$tooltip2 <- paste0(groupvar, "=", df[[colorname]])
    p <- ggplot(data = df, aes_string(
      x = "variable", y = "value",
      colour = colorname, fill = colorname, group = colorname
    )) +
      ggiraph::geom_polygon_interactive(aes_string(tooltip = "tooltip2"),
                                        alpha = alpha
      ) +
      ggiraph::geom_point_interactive(aes_string(
        data_id = id3,
        tooltip = "tooltip"
      ), size = size)
  }
  p
  if (!is.null(facetname)) {
    formula1 <- as.formula(paste0("~", facetname))
    p <- p + facet_wrap(formula1, scales = scales)
  }
  p <- p + xlab("") + ylab("") + theme(legend.position = legend.position)
  
  p <- p + coord_radar2(clip = clip)
  if (!is.null(ylim)) {
    p <- p + expand_limits(y = ylim)
  }
  p
  
  p
}
