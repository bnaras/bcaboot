## Generate hex sticker for bcaboot
library(ggplot2)
library(hexSticker)
library(sysfonts)
library(showtext)

font_add_google("Raleway", "raleway")
showtext_auto()

## BCa transformation curves: Phi(z0 + (z0+z)/(1-a*(z0+z)))
## Wide z range so sigmoids flatten well before the hex edges.
z <- seq(-3.5, 3.5, length.out = 400)

bca_curve <- function(z, z0, a) {
    numer <- z0 + z
    denom <- 1 - a * (z0 + z)
    pnorm(z0 + numer / denom)
}

curves <- do.call(rbind, lapply(seq(1, 5), function(i) {
    z0_val <- seq(-0.35, 0.35, length.out = 5)[i]
    a_val  <- seq(-0.14, 0.14, length.out = 5)[i]
    data.frame(
        x = z,
        y = bca_curve(z, z0_val, a_val),
        group = i
    )
}))

palette <- c("#e76f51", "#e63946", "#a8307e", "#6a0572", "#4361ee")

## Symmetric viewport: wide enough that all curves are flat at edges
p <- ggplot(curves, aes(x = x, y = y, group = factor(group),
                         colour = factor(group))) +
    geom_line(linewidth = 1.3, alpha = 0.8) +
    scale_colour_manual(values = palette) +
    coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-0.05, 1.05)) +
    theme_void() +
    theme(legend.position = "none",
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA))

sticker(
    p,
    package = "bcaboot",
    s_x = 1, s_y = 0.97,
    s_width = 1.65, s_height = 1.0,
    p_size = 20,
    p_x = 1, p_y = 1.52,
    p_color = "#3d3240",
    p_family = "raleway",
    p_fontface = "bold",
    h_fill = "#f5f0eb",
    h_color = "#a8307e",
    h_size = 1.4,
    filename = "man/figures/logo.png",
    dpi = 300
)

cat("Logo saved to man/figures/logo.png\n")
