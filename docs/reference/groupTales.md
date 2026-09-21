# Grouping TALEs

Classifying Tal groups by hierchical clustering or by k-medoids
clustering based on their similarity.

## Usage

``` r
groupTales(
  taleSim,
  plotTree = FALSE,
  k = NULL,
  k_test = NULL,
  method = "k-medoids"
)
```

## Arguments

- taleSim:

  a *three columns Tals similarity table* as obtained with
  [`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md)
  in the 'tal.similarity' slot of the returned object.

- plotTree:

  logical indicating whether to plot hclust tree or not. If the method
  is "k-medoids", no tree will be plotted (but instead, a plot of
  silhoutte value).

- k:

  integer indicating number of groups you want Tals to be classified. Or
  only in case that method is "k-medoids", k = "auto" to automatically
  pick the optimum k or k = NULL to interactively pick it. Do not always
  trust the automatic picking, it is better to choose k interactively or
  test with different values.

- k_test:

  integer vector of 2 indicating the range of k to test, only available
  when method = "k-medoids". Note that the minimum value for k is 2.

- method:

  one of two methods: "hclust" (see
  [`cutree`](https://rdrr.io/r/stats/cutree.html)) and "k-medoids" (see
  [`pam`](https://rdrr.io/pkg/cluster/man/pam.html)).

## Value

a data frame containing name of tals from taleSim and their classified
groups.
