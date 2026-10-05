# Download AnnoTALE's catalogue of TALE classes

A TALE class groups the TALEs of different strains that are similar
enough to be considered the same effector, and carries the systematic
name the literature uses for it (`TalAA`, `TalAB`, ...). AnnoTALE
curates that catalogue and publishes it; this is an R wrapper around its
'AnnoTALEcli.jar loadAndView' call, which downloads the current
definition.
[`run_annotale_assign`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md)
consumes what it writes.

## Usage

``` r
run_annotale_load_classes(
  output_dir = tempfile("annotale_classes_"),
  class_builder = NULL,
  opt_param = "",
  java_args = "-Xms512M -Xmx8G",
  annotale_jar = .tantale_tool("annotale")
)
```

## Arguments

- output_dir:

  Directory where the catalogue will be written, created if it does not
  exist. The default is a new directory under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), which R deletes
  when the session ends: given how long the download takes, give a path
  you keep.

- class_builder:

  Path to a class builder XML to read instead of downloading, as written
  by an earlier call. `NULL`, the default, downloads the current
  definition.

- opt_param:

  A single string of further options for the tool, as `key=value` pairs
  separated by spaces. The keys this function sets itself (`c`, `cb`,
  `outdir`) are refused: use `class_builder` and `output_dir`.

- java_args:

  A single string of options for the Java virtual machine, placed before
  `-jar`. The default raises the heap to 8 GB, which the catalogue
  needs; rebuilding it was seen to hold 1.8 GB resident.

- annotale_jar:

  Path to the AnnoTALE jar file. The default is the copy
  [`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
  downloads; give a path to use another version.

## Value

The path of the class builder XML, invisibly, so the call can be passed
straight to
[`run_annotale_assign`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md).
Beside it, `output_dir` also holds
`Lists_of_classes,_strains_and_TALEs/`, whose `List_of_classes.txt`
gives every class as plain text, one line per member with its aligned
repeat-variable diresidues, strain and systematic name. That file
answers "which class is this TALE in" without going through
[`run_annotale_assign`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md)
at all, when the TALE is already catalogued.

## Details

The download is slow and bulky: about a quarter of an hour, and some 440
MB of output, most of it the class builder itself, which carries a
cached alignment between every pair of catalogued TALEs. The file is
therefore kept once and reused rather than fetched per call, which is
what `output_dir` is for: point it at a directory you keep, and pass the
returned path to
[`run_annotale_assign`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md)
as often as you like.

Record the date along with any class you take from it. A class name is
stable, but the catalogue grows, and the index AnnoTALE appends to the
class when it names a TALE (`TalAH30`, the thirtieth member of `TalAH`)
is only true of the catalogue that produced it.

## See also

[`run_annotale_assign`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md),
which places TALEs into these classes.

Other external TALE tools:
[`run_annotale_assign()`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md),
[`run_annotale_build()`](https://scunnac.github.io/tantale/reference/run_annotale_build.md),
[`run_annotale_predict()`](https://scunnac.github.io/tantale/reference/run_annotale_predict.md)

## Examples

``` r
if (FALSE) { # interactive()
# Slow: downloads and rebuilds the catalogue, about 15 minutes.
classes <- run_annotale_load_classes(output_dir = tempfile("annotale_classes_"))
}
```
