# Blood and other measurements in diabetics

The `diabetes` data frame has 442 rows and 3 columns. These are the data
used in the Efron et al "Least Angle Regression" paper.

## Format

This data frame contains the following columns:

- **x** a matrix with 10 columns

- **y** a numeric vector

- **x2** a matrix with 64 columns

## Source

<https://web.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>

## Details

The x matrix has been standardized to have unit L2 norm in each column
and zero mean. The matrix x2 consists of x plus certain interactions.

## References

Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) *Annals of Statistics*
