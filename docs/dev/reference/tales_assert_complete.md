# Assert that a tales object holds complete arrays

Checks that every array carries all of its parts: `position_in_array`
running contiguously from 1, and — when `domain_type` and
`position_in_crd` are present — `position_in_crd` equal to
`position_in_array` minus the number of non-repeat parts before it.

## Usage

``` r
tales_assert_complete(x, arg = "x")
```

## Arguments

- x:

  A `tales` object.

- arg:

  Name of the argument being checked, for the error message.

## Value

`x`, invisibly.

## Details

This is a **precondition**, not an invariant:
`filter(x, domain_type == "repeat")` legitimately produces a valid
`tales` that is no longer complete. Alignment requires completeness,
because the mapping back from a MAFFT alignment is positional — the k-th
aligned residue is the k-th part fed in.

## See also

Other tales objects:
[`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/dev/reference/format.tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/format.tales_msa.md),
[`is_tales()`](https://scunnac.github.io/tantale/dev/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/dev/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/dev/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/print.tales_msa.md),
[`summary.tales()`](https://scunnac.github.io/tantale/dev/reference/summary.tales.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/summary.tales_msa.md),
[`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)

## Examples

``` r
rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
                         package = "tantale")
x <- as_tales(rvd_fasta, sep = "-")
tales_assert_complete(x) # every array is 1..n already -- no error

if (FALSE) { # \dontrun{
# Keeping only the repeats breaks completeness, which tales_align() needs.
repeats_only <- x[x$position_in_array > 1, ]
tales_assert_complete(repeats_only)
} # }
```
