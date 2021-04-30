#' ggplot theme developed for gte21 project
#' @return An object of class \code{\link[ggplot2]{theme}()}.
#'
#' @importFrom ggplot2 theme_void theme element_text element_rect element_line element_blank margin rel unit
#' @export
#' @family themes gte21
#' @rdname theme_gte21
theme_gte21 <- function(base_size=10) {
  theme_void() +
    theme(line = element_line(colour = "black"),
          rect = element_rect(fill = NA, colour = "black",
                              linetype = 1),
          text = element_text(colour = "black",family = "Arial"),
          
          ## Axis
          axis.line = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(margin = margin(t = base_size * 0.5,
                                                     unit = "pt")),
          axis.text.x.top = element_text(vjust = 0, margin = margin(b = base_size * 0.5, unit = "pt")),
          axis.text.y = element_text(hjust = 1, margin = margin(r = base_size * 0.25,
                                                                unit = "pt")),
          axis.ticks = element_line(lineend="round"),
          axis.title = element_text(size = rel(2)),
          axis.title.x = element_text(margin = margin(t=base_size*0.25)),
          axis.title.y = element_text(angle = 90),
          axis.ticks.length = unit( base_size * 0.3, "points"),
          
          #legend
          legend.background = element_rect(linetype = 1),
          legend.spacing = unit(base_size * 0.5, "points"),
          legend.key = element_rect(linetype = 0),
          legend.key.size = unit(0.5, "lines"),
          legend.key.height = NULL,
          legend.key.width = NULL,
          legend.text = element_text(size = rel(0.75)),
          legend.text.align = NULL,
          legend.title = element_text(size = rel(1),  hjust = 0.5),
          legend.title.align = NULL,
          legend.position = "right",
          legend.direction = "vertical",
          legend.justification = "center",
          
          # panel
          panel.background = element_rect(fill = NA),
          panel.border = element_rect(linetype = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1, "lines"),
          
          # strip
          strip.background = element_rect(fill = NA,
                                          colour = NA, linetype = 0),
          strip.text = element_text(size = rel(1.25)),
          strip.text.x = element_text(),
          strip.text.y = element_text(angle = -90),
          
          # overall/panel
          plot.background = element_rect(fill = NA,
                                         colour = NA),
          plot.title = element_text(size = rel(2.5),
                                    hjust = 0.5,margin = margin(b=base_size),
                                    face = "bold"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          aspect.ratio = NULL,
          complete = TRUE)
}


.gte21_palettes <- list(
  binary = c("red2","gray"),
  binary2 = c("lightgray","darkgray"),
  diverging = colorRamps::matlab.like2(10),
  diverging.wide = grDevices::rainbow(20),
  categories1 = RColorBrewer::brewer.pal(12, name = "Set3"),
  categories2 = head(colorRamps::primary.colors(no.white = T),15)[-1],
  increasing = grDevices::topo.colors(10)
)

#' color palettes used for gte21 project
#'
#' binary: red, gray
#'
#' diverging: colorRamps::matlab.like2
#'
#' diverging.wide: grDevices::rainbow
#'
#' categories1: colorRamps::primary.colors
#'
#' increasing: grDevices::topo.colors
#'
#' @export
#' @family color gte21
gte21_pal <-  function(palette = "categories1", reverse = FALSE, ...) {
  pal <- .gte21_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  grDevices::colorRampPalette(pal, ...)
}

#' color scales used for gte21 project
#'
#' @importFrom ggplot2 discrete_scale scale_color_gradientn
#' @export
#' @family color gte21
#' @rdname scale_color_gte21
scale_color_gte21 <- function(palette = "categories1", discrete = TRUE, reverse = FALSE, ...) {
  pal <- gte21_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("gte21_color_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' @importFrom ggplot2 discrete_scale scale_color_gradientn
#' @export
#' @family color gte21
#' @rdname scale_color_gte21
scale_fill_gte21 <- function(palette = "categories1", discrete = TRUE, reverse = FALSE, ...) {
  pal <- gte21_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("gte21_fill_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}
