# Multiple alignment of TALE arrays

The [previous
article](https://scunnac.github.io/tantale/articles/tale_classification.md)
grouped TALE arrays by overall similarity. Grouping says arrays are
related. Saying *how* (which repeats correspond, where one array has an
insertion or a deletion relative to another) needs an alignment.

Code

``` r
library(tantale)
library(dplyr)
```

> **Where this fits**
>
> This article is a walkthrough: pick one real classification group and
> align it. The mechanics (what a `tales_msa` actually is, coercion
> between `tales` and `tales_msa`, and the plotting options in full) are
> a deep dive in [the `tales_msa` class
> article](https://scunnac.github.io/tantale/articles/tales_msa_class.md),
> which continues directly from the alignment built here.

The AnnoTALE tool this package builds on can assign TALEs to classes,
but cannot insert gaps while doing so. TALE arrays evolve substantially
by whole-repeat duplication and deletion, so a method that cannot
represent a gap is blind to exactly the events that matter most (see the
[QueTAL paper](https://doi.org/10.3389/fpls.2015.00545)).
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md)
drives MAFFT in text mode instead, treating each distinct domain as one
alignable symbol, which lets it open and close gaps freely.

## 1 The same group, from the previous article

Same three genomes (MAI1, BAI3, BAI3-1-1), same correction settings, and
the same comparison and grouping as [the classification
article](https://scunnac.github.io/tantale/articles/tale_classification.md),
reused here from its cached result; see [the getting-started
article](https://scunnac.github.io/tantale/articles/getting_started.md)
for how these four articles are linked.

This group holds one member from each of MAI1, BAI3 and BAI3-1-1 (the
same locus, present in all three related African strains) and, usefully
for this article, its three copies are not identical:

Code

``` r
picked_group <- grouped |> filter(group == 6)
picked_group |> distinct(array_id)
#> # A tibble: 3 × 1
#>   array_id          
#>   <chr>             
#> 1 MAI1_ROI_00007    
#> 2 BAI3_ROI_00007    
#> 3 BAI3-1-1_ROI_00006
```

## 2 Aligning the group

Code

``` r
msa <- tales_align(picked_group, residue_col = "dom_code")
msa
#> <tales_msa> 3 arrays, 18 alignment positions
#>   layers: rvd, dom_code   |   namespace: 5d8762602564ac93   |   8 other columns
#>                       dom_code
#>   MAI1_ROI_00007       84  35  43   5  41  31  24  49  35  54  18  40  31  ...
#>   BAI3_ROI_00007       84  35  43   5  41  31  24  49  35  54  18  40   -  ...
#>   BAI3-1-1_ROI_00006   85  35  43   5  41  31  24  49  35  54  18  40   -  ...
```

Code

``` r
plot(msa, tale_distances = cmp$tale_distances, domain_distances = cmp$domain_distances)
```

[![](tale_msa_files/figure-html/fig-msa-default-1.png)](https://scunnac.github.io/tantale/articles/tale_msa_files/figure-html/fig-msa-default-1.png "Figure 1: Alignment of one TALE locus across three related strains, coloured by domain cluster. MAI1’s array (marked #) is the reference.")

Figure 1: Alignment of one TALE locus across three related strains,
coloured by domain cluster. MAI1’s array (marked \#) is the reference.

[Figure 1](#fig-msa-default) already shows something
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md)’s
single distance number could not: **BAI3 and BAI3-1-1 both lack four
repeats that MAI1 has**, at alignment positions 13-16. The gap is
internal: all three arrays end with the same 20-residue half-repeat
(column 17), which lines up across the gap. BAI3-1-1 derives from BAI3,
so the two share this difference.

## 3 Does a scoring matrix change the alignment?

[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md)
accepts a `domain_distances` argument: a substitution cost matrix MAFFT
uses when scoring which repeats to match against each other. Without it,
every mismatch counts as equally bad. Passing one is optional, and it is
worth checking what it actually changes.

### 3.1 Aligning on `dom_code`

Without one, every distinct domain is as different from every other as
any other pair, and MAFFT has no notion that two repeats might be more
or less alike:

Code

``` r
msa_plain <- tales_align(picked_group, residue_col = "dom_code")
#> Now running MAFFT (Copyright 2002-2007 Kazutaka Katoh) on TALE array sequences.
```

With one, `domain_distances`, the domain-level protein similarity
already computed by
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
tells MAFFT how alike two repeats actually are:

Code

``` r
msa_scored <- tales_align(picked_group, residue_col = "dom_code",
                          domain_distances = cmp$domain_distances)
```

Code

``` r
identical(as.matrix(msa_plain), as.matrix(msa_scored))
#> [1] TRUE
```

On this group the two alignments are identical. The domain distances
count the residues a half-repeat lacks against a full repeat, so the
final half-repeat that BAI3 and BAI3-1-1 share with MAI1 stays opposite
its identical copy (column 17), and the four-repeat gap stays internal.
On more divergent arrays, the matrix lets MAFFT prefer matching similar
repeats over dissimilar ones; whether that changes the alignment depends
on the data, so it is worth checking as done here.

### 3.2 Aligning on `rvd`

The built-in RVD similarity matrix works the same way, opted into with
`domain_distances = "rvd"` instead of a distance table. It scores
DNA-binding specificity instead of protein sequence, so it is worth
checking separately:

Code

``` r
msa_rvd_plain  <- tales_align(picked_group, residue_col = "rvd")
msa_rvd_scored <- tales_align(picked_group, residue_col = "rvd", domain_distances = "rvd")
identical(as.matrix(msa_rvd_plain), as.matrix(msa_rvd_scored))
#> [1] TRUE
```

On this group, aligning on RVDs also gives the *same* alignment whether
or not the matrix is supplied. The RVD alphabet is smaller than the
domain-code one (6 distinct RVDs against 12 distinct repeat codes here),
and these three arrays are closely related, so there was little room for
a scoring matrix to change anything. Whether either matrix matters
depends on the data.

## 4 Other views of the same alignment

Colouring by similarity to the reference array instead of by cluster
membership makes graded divergence visible where cluster identity would
only say “different”:

Code

``` r
plot(msa, fill_type = "domain_sim",
    tale_distances = cmp$tale_distances, domain_distances = cmp$domain_distances)
```

[![](tale_msa_files/figure-html/fig-msa-sim-1.png)](https://scunnac.github.io/tantale/articles/tale_msa_files/figure-html/fig-msa-sim-1.png "Figure 2: The same alignment, coloured by protein-sequence similarity to the reference array.")

Figure 2: The same alignment, coloured by protein-sequence similarity to
the reference array.

In [Figure 2](#fig-msa-sim), every domain is identical to MAI1’s except
BAI3-1-1’s N-terminus, which has its own `dom_code` at about 99.7%
similarity. The cluster colouring of [Figure 1](#fig-msa-default) groups
it with the other two N-termini, so the difference only shows here.

With a consensus panel attached, the column-by-column question “does
this array agree with the group” is answered directly above the
alignment itself:

Code

``` r
plot(msa, consensus = TRUE,
    tale_distances = cmp$tale_distances, domain_distances = cmp$domain_distances)
```

[![](tale_msa_files/figure-html/fig-msa-consensus-1.png)](https://scunnac.github.io/tantale/articles/tale_msa_files/figure-html/fig-msa-consensus-1.png "Figure 3: The same alignment again, with a consensus row attached above it.")

Figure 3: The same alignment again, with a consensus row attached above
it.

## 5 Next

Once related arrays are aligned, the biologically motivated next
question is what DNA sequence their RVDs are predicted to bind, covered
in [the final
article](https://scunnac.github.io/tantale/articles/tale_target_prediction.md)
of this series. To go deeper on the alignment mechanics used here first,
see [the `tales_msa` class
article](https://scunnac.github.io/tantale/articles/tales_msa_class.md).
