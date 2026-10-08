# Assign TALEs to AnnoTALE's published classes

An R wrapper around the 'AnnoTALEcli.jar assign' call. AnnoTALE places
each TALE in the class of the catalogue it belongs to, and opens a new
class for any that fit none. Where
[`run_annotale_build`](https://scunnac.github.io/tantale/reference/run_annotale_build.md)
groups the TALEs you hand it among themselves, knowing nothing of what
the rest of the world calls them, this one answers "which published
class is this TALE in", so its names can be compared with the
literature's.

## Usage

``` r
run_annotale_assign(
  fasta_file,
  class_builder,
  output_dir = tempfile("annotale_assign_"),
  strain = NULL,
  accession = NULL,
  opt_param = "",
  java_args = "-Xms512M -Xmx8G",
  annotale_jar = .tantale_tool("annotale")
)
```

## Arguments

- fasta_file:

  Path to a FASTA file of the TALEs to assign, one record per TALE (or
  per TALE part), in one of these forms:

  - full-length TALE coding sequences, such as the
    `Predict/TALE_DNA_sequences_*.fasta` file
    [`run_annotale_predict`](https://scunnac.github.io/tantale/reference/run_annotale_predict.md)
    writes;

  - the corresponding protein sequences, such as its
    `Predict/TALE_protein_sequences_*.fasta`;

  - the parts of each TALE, one record per N-terminus, repeat and
    C-terminus, as in its `Analyze/TALE_DNA_parts.fasta` or
    `Analyze/TALE_Protein_parts.fasta`;

  - RVD sequences, the RVDs separated by hyphens (`NI-HD-NG-NN-...`), as
    in its `Analyze/TALE_RVDs.fasta` or the shipped
    `Sample_TALEs_RVDSeqs_AnnoTALE.fasta`
    (`system.file("extdata", package = "tantale")`).

- class_builder:

  Path to the class builder XML holding the classes to assign against,
  as
  [`run_annotale_load_classes`](https://scunnac.github.io/tantale/reference/run_annotale_load_classes.md)
  returns.

- output_dir:

  Directory where output will be written, created if it does not exist.
  The default is a new directory under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), which R deletes
  when the session ends: give a path to keep the results.

- strain:

  The strain the TALEs come from. AnnoTALE uses it to build the
  systematic names it proposes, so without it they are unnamed.

- accession:

  The genome's accession number, recorded in the report.

- opt_param:

  A single string of further options, as `key=value` pairs separated by
  spaces. The keys this function sets itself (`t`, `c`, `s`, `a`,
  `outdir`) are refused: use `fasta_file`, `class_builder`, `strain`,
  `accession` and `output_dir`.

- java_args:

  A single string of options for the Java virtual machine, placed before
  `-jar`. The default raises the heap to 8 GB, which reading the
  catalogue needs; raise `-Xmx` if it runs out of memory.

- annotale_jar:

  Path to the AnnoTALE jar file. The default is the copy
  [`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
  downloads; give a path to use another version.

## Value

`output_dir`, invisibly. It holds the per-TALE assignment report, a
table of original against proposed systematic names, a report and figure
for each class that changed or was created, and a class builder extended
with the TALEs given.

## See also

[`run_annotale_load_classes`](https://scunnac.github.io/tantale/reference/run_annotale_load_classes.md),
which fetches the classes;
[`run_annotale_build`](https://scunnac.github.io/tantale/reference/run_annotale_build.md),
which builds classes from scratch instead.

Other external TALE tools:
[`run_annotale_build()`](https://scunnac.github.io/tantale/reference/run_annotale_build.md),
[`run_annotale_load_classes()`](https://scunnac.github.io/tantale/reference/run_annotale_load_classes.md),
[`run_annotale_predict()`](https://scunnac.github.io/tantale/reference/run_annotale_predict.md)

## Examples

``` r
if (FALSE) { # interactive()
# Needs a Java runtime and the catalogue, which takes a while to fetch.
classes <- run_annotale_load_classes(output_dir = tempfile("annotale_classes_"))
predicted <- run_annotale_predict(tantale_genome("MAI1"))
tales <- list.files(file.path(predicted, "Predict"),
                    pattern = "^TALE_DNA_sequences_", full.names = TRUE)
run_annotale_assign(tales, class_builder = classes, strain = "MAI1")
}
```
