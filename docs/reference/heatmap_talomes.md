# Heatmap plotting of rvd sequence variants

The function creates a graphical presentation from a tale annotation
table. The output is like a heatmap that presents rvd sequence variants
in Tal groups as column and respective strains as rows (or vice versa).
It is different from a typical heatmap that it can display more than one
value in a cell; for example, if one strain has 2 rvd sequence variants
belong to 1 group, it will be displayed by 2 colors in 1 cell.

## Usage

``` r
heatmap_talomes(
  tale_annotation,
  col,
  row,
  value,
  truncTaleLab = NULL,
  extraCol = NULL,
  x.lab = "TALE Group",
  y.lab = "Strain",
  title = "RVD sequences variants",
  plot.type = "all",
  mapcol = viridis::viridis(10),
  mar.side = c(5, 5, 3, 3),
  sepwid = 5,
  sepcol = "white",
  inner_sepcol = "white",
  save.path = NULL
)
```

## Arguments

- tale_annotation:

  a data frame containing at least 3 columns for Tal groups, strain
  names, and rvd seqs, and 1 row is 1 Tal.

- col:

  "character", column name of `tale_annotation` to be displayed as rows
  in the heatmap (e.g. tal groups).

- row:

  "character", column name of `tale_annotation` to be displayed as
  columns in the heatmap (e.g. strain names).

- value:

  "character", column name for rvdseqs in the `tale_annotation`

- truncTaleLab:

  (optional, default = NULL) "character", column name of
  `tale_annotation` labeling the truncTales by TRUE/FALSE value. The
  truncTales are labeled by "T" in the heatmap cells, but if this
  argument is called.

- extraCol:

  (optional, default = NULL) "character", column name of
  `tale_annotation` containing other information (e.g. origin). It will
  be presented in a side bar on the right of the heatmap.

- x.lab, y.lab, title:

  character for x axix, y axis names and title

- mapcol:

  character vector of colors for the cells.

- mar.side:

  margin of the heatmap for row dendrogram, col dendrogram, rownames,
  colnames, respectively. (by default, c(5, 5, 3, 3)).

- sepwid:

  numeric value for the width of separator between adjacent cells.

- sepcol:

  character of color for the separator between adjacent cells

- inner_sepcol:

  character of color for the separator between colors within 1 cell if
  there are more than 1.

- save.path:

  (optional) file path to save the plot, format of the image depends on
  the file extension. If save.path is NULL, the heatmap will be printed.
  If save.path is specified, the image file will be created.
