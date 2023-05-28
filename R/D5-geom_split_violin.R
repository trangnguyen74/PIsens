
# The code for ggproto class .GeomSplitViolin and function .geom_split_violin here are adapted from https://stackoverflow.com/a/45614547/2121635
# The adaptation is to add a "nudge" parameter to insert some space between the two halves of the violin. This param is borrowed from the gghalves package at https://github.com/erocoar/gghalves.


########################################
#### .geom_split_violin

#' (internal) geom function for split violin plot
#'
#' @inheritParams ggplot2::geom_violin
#' @param nudge Add space between the two halves of the violin
#' @keywords internal

.geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              nudge = 0,
                              ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {

    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = stat,
                   geom = .GeomSplitViolin,
                   position = position,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(trim = trim,
                                 scale = scale,
                                 nudge = nudge,
                                 draw_quantiles = draw_quantiles,
                                 na.rm = na.rm,
                                 ...))
}




########################################
#### .GeomSplitViolin

#' (internal) ggproto class for split violin plot
#'
#' @keywords internal
#'
.GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin",
    ggplot2::GeomViolin,
    draw_group = function(self,
                          data,
                          ...,
                          nudge = 0,
                          draw_quantiles = NULL) {

        data <- transform(data,
                          xminv = x - violinwidth * (x - xmin),
                          xmaxv = x + violinwidth * (xmax - x))
        grp <- data[1, "group"]
        newdata <- plyr::arrange(transform(data,
                                           x = if (grp %% 2 == 1) xminv else xmaxv),
                                 if (grp %% 2 == 1) y else -y)
        newdata <- rbind(newdata[1, ],
                         newdata,
                         newdata[nrow(newdata), ],
                         newdata[1, ])
        newdata[c(1,
                  nrow(newdata) - 1,
                  nrow(newdata)),
                "x"] <- round(newdata[1, "x"])

        newdata$x <- ifelse(newdata$group %% 2 == 1,
                            newdata$x - nudge,
                            newdata$x + nudge)

        if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {

            stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))

            quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                                 draw_quantiles)
            aesthetics <- data[rep(1, nrow(quantiles)),
                               setdiff(names(data), c("x", "y")),
                               drop = FALSE]
            aesthetics$alpha <- rep(1, nrow(quantiles))
            both <- cbind(quantiles, aesthetics)
            quantile_grob <- ggplot2::GeomPath$draw_panel(both,
                                                          ...)
            ggplot2:::ggname("geom_split_violin",
                             grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata,
                                                                            ...),
                                            quantile_grob))
        }
        else {
            ggplot2:::ggname("geom_split_violin",
                             ggplot2::GeomPolygon$draw_panel(newdata,
                                                             ...))
        }
    }
)



