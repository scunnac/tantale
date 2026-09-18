# Predict TALE target boxes

Runs a TALE target (EBE) prediction tool over a set of TALE RVD
sequences and a set of DNA sequences, and returns the predictions as a
tibble.

## Usage

``` r
tales_predict_targets(x, subj_file, method = c("talvez", "preditale"), ...)
```

## Arguments

- x:

  TALE RVD sequences: a
  [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object, the path to a fasta file, or a `BStringSet`. A `tales` is
  rendered with
  [`tales_rvd_strings`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md),
  which drops the termini.

- subj_file:

  Path to a fasta file of subject DNA sequence(s).

- method:

  Which backend to use, `"talvez"` (default) or `"preditale"`.

- ...:

  Passed to the chosen backend. See
  [`talvez`](https://scunnac.github.io/tantale/dev/reference/talvez.md)
  (`opt_param`, `talvez_dir`, `conda_bin`) and
  [`preditale`](https://scunnac.github.io/tantale/dev/reference/preditale.md)
  (`opt_param`, `predictor_path`); both accept `output_dir`. Note the
  two take different `opt_param` defaults, since the options are the
  tools' own.

## Value

A tibble of EBE predictions, with column names homogenised across
backends, plus a `method` column recording which tool produced them.

## Details

`talvez` and `preditale` are independent programs, but the package
already normalises their outputs to a shared set of column names, so
they are interchangeable backends of one operation rather than two
separate functions. This is that operation;
[`talvez`](https://scunnac.github.io/tantale/dev/reference/talvez.md)
and
[`preditale`](https://scunnac.github.io/tantale/dev/reference/preditale.md)
remain available and documented individually, and are where each tool's
own options and citation live.

## See also

[`talvez`](https://scunnac.github.io/tantale/dev/reference/talvez.md),
[`preditale`](https://scunnac.github.io/tantale/dev/reference/preditale.md)

Other target prediction:
[`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md),
[`preditale()`](https://scunnac.github.io/tantale/dev/reference/preditale.md),
[`talvez()`](https://scunnac.github.io/tantale/dev/reference/talvez.md)
