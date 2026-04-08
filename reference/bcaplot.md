# Plots of bca confidence limits

`bcaplot` uses the output of `bcajack`, `bcajack2`, `bcanon`, or
`bcapar` to plot bca and standard confidence limits for the parameter of
interest.

## Usage

``` r
bcaplot(
  vl,
  main = " ",
  xlab = "coverage",
  ylab = "limits",
  alpha = c(0.025, 0.05, 0.1, 0.16),
  ylim,
  xlim,
  add = 0,
  sub = "black=bca, green=standard",
  sw = 1,
  ...
)
```

## Arguments

- vl:

  output of `bcajack`, `bcajack2`, `bcanon`, or `bcapar`

- main:

  The main caption (can be empty)

- xlab:

  The x axis label (supplied if not specified)

- ylab:

  The y axis labels (supplied if not specified)

- alpha:

  Coverages are \\1-2\alpha\\, e.g. `alpha=c(.025,.05)` plots intervals
  `[.025,.975]` and `[.05,.95]`. Default is `alpha=c(.025,.05,.1,.16)`
  giving coverages .95,.90,.80,.68

- ylim:

  y axis plot limits set automatically if not provided

- xlim:

  x axis plot limits set automatically if not provided

- add:

  `add=1` adds a new plot of bca limits (in red) to an existing plot

- sub:

  subtitle (can be empty)

- sw:

  `sw=1` draws light vertical dashed lines showing the bca intervals

- ...:

  further args for plot

## Details

confidence interval endpoints are plotted vertically versus two-sided
coverages \\1-2\alpha\\. Bca limits in black, Standard limits in green
(dashed.). If `vl$lims` includes the column "jacksd" of jackknife
internal standard deviations then these are indicated by vertical red
bars centered at the bca limit points.
