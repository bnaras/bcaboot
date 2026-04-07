# Automatic Construction of Bootstrap Confidence Intervals

Bootstrap confidence intervals depend on three elements: (a) the
cumulative distribution of the bootstrap replications, (b) the
bias-correction, and (c) the acceleration number that measures the rate
of change in the standard deviation of the estimate as the data changes.
The first two of these depend only on the bootstrap distribution, and
not how it is generated: parametrically or non-parametrically.
Therefore, the only difference in a parametric bca analysis would lie in
the nonparametric estimation of the acceleration, often a negligible
error.

## See also

Useful links:

- <https://bnaras.github.io/bcaboot/>

- <https://github.com/bnaras/bcaboot>

- Report bugs at <https://github.com/bnaras/bcaboot/issues>

## Author

**Maintainer**: Balasubramanian Narasimhan <naras@stat.stanford.edu>

Authors:

- Bradley Efron <brad@stat.stanford.edu>
