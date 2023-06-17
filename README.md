# MINMI - sample MINimum Monte carlo Inversion to estimate extinction time

This package estimates extinction time using a dataset of fossil ages,
together with measurement error standard deviations. The endpoint of the interval, the oldest date at which a fossil could potentially have been
included in the dataset, is also required. Recovery potential is assumed to be constant over the interval, and uncertainty estimating
fossil ages is assumed to be independent across specimens and normally distributed. Inversion of the sample minimum is used to find
and exact confidence interval for extinction time, under the above assumptions.

## How it works

## Installation instructions

Install the **development** version of `rminmi` from [GitHub](https://github.com/):
```{r install}
# install.packages("remotes")
remotes::install_github("victorwctsang/rminmi")

library(rminmi)
```

### Spot a bug?
Thanks for finding the bug! We would appreciate it if you can pop over to our [Issues page](https://github.com/victorwctsang/rminmi/issues) and describe how to reproduce the bug!

## References

## Credits



