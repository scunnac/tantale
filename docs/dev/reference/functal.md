# Run functal from QueTAL to build a phylogenetic tree of TALE RVD sequences.

A R wrapper around the [QueTAL](https://doi.org/10.3389/fpls.2015.00545)
'functal' perl script.

## Usage

``` r
functal(
  tal_file,
  tree_format = "fan",
  output_prefix = "FuncTALE",
  output_dir = getwd(),
  functal_path = system.file("tools", "QueTAL_v1.1", "FuncTAL", "FuncTAL_v.1.1.pl",
    package = "tantale", mustWork = T),
  conda_bin = "auto"
)
```

## Arguments

- tal_file:

  Path to a QueTAL-formatted file of TALE RVD sequences.

- tree_format:

  Tree layout passed to functal's `-n` option (default `"fan"`).

- output_prefix:

  Prefix used for functal's output file names.

- output_dir:

  Directory where output will be copied (default: current working
  directory).

- functal_path:

  Path to the functal perl script if you want to use another version
  than the one provided with tantale.

- conda_bin:

  Path to your Conda binary file if you need to specify a non-standard
  location, otherwise leave to "auto".

## Value

Returns invisibly the exit code of the shell call to functal (ie '0' if
successful).

## See also

Other external TALE tools:
[`run_annotale_build()`](https://scunnac.github.io/tantale/dev/reference/run_annotale_build.md),
[`run_annotale_predict()`](https://scunnac.github.io/tantale/dev/reference/run_annotale_predict.md)
