# Compare TALEs by predicted DNA-binding specificity (FuncTAL)

Quantifies how TALE arrays relate by the DNA sequence their repeats are
predicted to bind, rather than by domain sequence identity
([`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)).
Each array's repeats are turned into a position weight matrix (PWM) over
the RVD-to-base specificity code, and PWMs are compared pairwise with
[`compare_motifs`](https://rdrr.io/pkg/universalmotif/man/compare_motifs.html).

## Usage

``` r
tales_compare_functal(
  x,
  method = "PCC",
  tryRC = FALSE,
  min.overlap = 1,
  normalise.scores = TRUE,
  score.strat = "a.mean",
  nthreads = 1
)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object carrying an `rvd` column.

- method:

  One of `compare_motifs()`'s comparison metrics (`"PCC"` default).
  `"PCC"`, `"WPCC"`, `"SW"`, `"ALLR"`, `"ALLR_LL"` and `"BHAT"` are
  similarities (higher is more similar); the rest are already distances.
  Either kind is turned into a proper `dissim` in the result.

- tryRC:

  Also score each pair's reverse-complement alignment and keep the
  better one. `FALSE` by default: an RVD array's specificity code has a
  fixed reading direction (5' to 3' target, N- to C-terminal repeat
  order), so comparing against a reverse complement asks a different,
  narrower biological question – do these two TALEs target opposite
  strands of overlapping sites – worth asking on purpose, not folded
  silently into every comparison.

- min.overlap:

  Minimum aligned width to accept, as `compare_motifs()` defines it.
  Defaults to `1` (any overlap at all), not `compare_motifs()`'s own
  default of `6` – generic TF motifs are usually longer than 6
  positions, but a TALE array legitimately has as few as half a dozen
  repeats, so the upstream default would silently refuse to compare some
  real, short arrays.

- normalise.scores:

  Penalise alignments that leave much of either motif unaligned. `TRUE`
  by default: TALE array lengths vary widely (see `min.overlap`), so an
  unnormalised score would favour matches that are mostly one array
  hanging off the end of a much longer one.

- score.strat:

  How `compare_motifs()` combines per-column scores into one alignment
  score. `"a.mean"` (its own default): a sum would scale with array
  length, favouring long-vs-long comparisons for no biological reason.

- nthreads:

  Passed to `compare_motifs()`.

## Value

A
[`tale_distances`](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
object – interchangeable with
[`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)'s,
so it can be handed directly to
[`tales_group_hclust`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md)/[`tales_group_kmedoids`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md).
No tree is built here, matching
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)'s
own division of labour: compare here, cluster/tree there.

## Details

This is a reimplementation, not a port, of QueTAL's FuncTAL: the
original Perl tool could not be kept working (it needs `Bio::Perl`,
dropped by BioPerl's 1.7 reorganisation; see
`dev/restructuring-notes.md` §12b), so the comparison itself was rebuilt
on `universalmotif` rather than patched. **Results diverge from the
original FuncTAL tool, and this is by design, not an approximation to be
improved away.** FuncTAL scored two RVD arrays by flattening their
entire padded, overlapping alignment (positions and bases together) into
one vector and taking a single Pearson correlation. `compare_motifs()`
instead correlates matched columns individually and combines the column
scores (`score.strat`). Verified empirically to disagree on real data
before writing this function, not assumed to differ only in magnitude.

The PWMs themselves are built by
[`tales_to_universalmotif`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)
– see its docs for exactly what drives them (only `rvd`, in repeat
order, termini dropped) and how an RVD outside
[`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md)
is handled. That conversion is exposed on its own precisely so it is not
locked inside this one comparison.

Only a handful of
[`compare_motifs`](https://rdrr.io/pkg/universalmotif/man/compare_motifs.html)'s
many options are exposed here, chosen for what actually varies across
TALE arrays (their differing repeat counts) rather than mirrored
wholesale; see `?compare_motifs` for the rest, several of which are
worth revisiting – ledger §12b lists them as follow-ups.

## See also

[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
comparing by domain sequence instead;
[`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md),
the conversion step this composes – called directly, any other
`universalmotif` function ( `motif_tree()`, `view_motifs()`,
`scan_sequences()`, `merge_motifs()`, ...) can be run on the same PWMs.

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md),
[`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md),
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/reference/tales_domain_distances.md),
[`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md),
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
fn <- tales_compare_functal(x)
fn
#> # A tibble: 16 × 3
#>    id1       id2       dissim
#>    <chr>     <chr>      <dbl>
#>  1 ROI_00001 ROI_00001  0    
#>  2 ROI_00002 ROI_00001  0.763
#>  3 ROI_00003 ROI_00001  0.785
#>  4 ROI_00004 ROI_00001  0.714
#>  5 ROI_00001 ROI_00002  0.763
#>  6 ROI_00002 ROI_00002  0    
#>  7 ROI_00003 ROI_00002  0.696
#>  8 ROI_00004 ROI_00002  0.728
#>  9 ROI_00001 ROI_00003  0.785
#> 10 ROI_00002 ROI_00003  0.696
#> 11 ROI_00003 ROI_00003  0    
#> 12 ROI_00004 ROI_00003  0.712
#> 13 ROI_00001 ROI_00004  0.714
#> 14 ROI_00002 ROI_00004  0.728
#> 15 ROI_00003 ROI_00004  0.712
#> 16 ROI_00004 ROI_00004  0    
```
