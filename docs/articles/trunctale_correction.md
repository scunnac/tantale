# Genuine truncTALEs and frameshift correction

Some TALEs are naturally short. [Ji et
al. (2016)](https://doi.org/10.1038/ncomms13435) and [Read et
al. (2016)](https://doi.org/10.3389/fpls.2016.01516) independently
characterised truncated TAL effectors in *Xanthomonas oryzae*:
“interfering TALEs” (iTALEs) in the first paper, “truncTALEs” in the
second, for what both describe as the same underlying phenomenon. Both
report loss of the C-terminal transcription activation domain while the
central repeat region is retained; Read et al. additionally describe
loss of part of the N-terminal region and the second nuclear
localisation signal in most examples, and a novel 28-amino-acid repeat
length not otherwise seen in TALEs. Functionally, the characterised
examples suppress the rice resistance gene *Xa1*. Read et al.’s Tal2h
did not bind any tested candidate DNA target, consistent with a decoy or
dominant-negative role for the plant immune receptor rather than
DNA-binding transactivation. Both papers document truncTALEs directly in
PXO86 (named `Tal3`/`Tal6` in Ji et al.), the genome used throughout
this article.

That biology creates a hazard for automatic frameshift correction.
[Mining TALE
sequences](https://scunnac.github.io/tantale/articles/tale_mining.md)
already shows that neither of tantale’s two correction paths is a silver
bullet against real assembly errors. This article asks a different
question: what do they do to a TALE that is short for real biological
reasons, with nothing wrong to fix?

> **Where this fits**
>
> This article extends the correction methods introduced in [Mining TALE
> sequences](https://scunnac.github.io/tantale/articles/tale_mining.html#correcting-frameshifts-two-ways):
> read that section first. This article is a side branch of that
> walkthrough, for when you need to reason about a TALE that is
> genuinely short.

Two concrete questions follow from the hazard above. Does
`tell_tales(correct_array = TRUE)` treat a genuine truncTALE the way it
is already known to treat a genuine assembly error, extending it
regardless of whether extension reflects real biology? And does
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md),
the genome-wide correction route, behave any differently? PXO86,
documented in the literature as carrying truncTALEs and shipped with the
package, answers both directly on real discovered arrays.

Code

``` r
library(tantale)
library(dplyr)
```

Code

``` r
out <- fs::dir_create(file.path(tempdir(), "trunctale_correction"))
```

## 1 Two genuine truncTALEs in one genome

PXO86 is one of the good-quality assemblies shipped with the package;
BAI3-1-1, used in the mining article, is the deliberately error-prone
one. PXO86 also carries two documented truncTALEs, which makes it the
right genome to answer the questions above. **On a good assembly, a
short array deserves to be taken at face value.**

Code

``` r
pxo86_fa <- tantale_genome("PXO86")
pxo86_raw_dir <- file.path(out, "PXO86_raw")
```

Code

``` r
invisible(tell_tales(subject_file = pxo86_fa, output_dir = pxo86_raw_dir))
pxo86 <- suppressWarnings(tales_from_telltales(pxo86_raw_dir))
```

Comparing the amino-acid length of every array’s N- and C-terminus finds
them immediately: two arrays with a shortened C-terminus and a reduced
N-terminus:

Code

``` r
pxo86 %>%
  filter(domain_type %in% c("N-terminus", "C-terminus")) %>%
  mutate(aa_length = nchar(aa_seq)) %>%
  select(array_id, domain_type, aa_length) %>%
  tidyr::pivot_wider(names_from = domain_type, values_from = aa_length) %>%
  arrange(`C-terminus`)
#> # A tibble: 18 × 3
#>    array_id  `N-terminus` `C-terminus`
#>    <chr>            <int>        <int>
#>  1 ROI_00019          230           42
#>  2 ROI_00001          230          183
#>  3 ROI_00002          283          286
#>  4 ROI_00004          288          286
#>  5 ROI_00006          288          286
#>  6 ROI_00007          288          286
#>  7 ROI_00008          288          286
#>  8 ROI_00009          288          286
#>  9 ROI_00010          288          286
#> 10 ROI_00011          288          286
#> 11 ROI_00012          288          286
#> 12 ROI_00013          288          286
#> 13 ROI_00014          288          286
#> 14 ROI_00015          288          286
#> 15 ROI_00016          288          286
#> 16 ROI_00018          288          286
#> 17 ROI_00003          288          297
#> 18 ROI_00017          288          297
```

`ROI_00019`’s C-terminus is barely a seventh of the 286-297 aa every
other array carries. `ROI_00001`’s is shorter than normal, at 183 aa,
but far less drastically so. Both arrays share the same reduced
N-terminus, 230 aa against 283-288 aa elsewhere, in line with the
partial N-terminal loss Read et al. describe.
[`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md)
reports all four termini:

Code

``` r
tales_anomalies(pxo86)
#> # A tibble: 4 × 4
#>   array_id  check                 kind         detail                           
#>   <chr>     <chr>                 <chr>        <chr>                            
#> 1 ROI_00001 terminus_noncanonical noncanonical N-terminus is not a canonical TA…
#> 2 ROI_00001 terminus_noncanonical noncanonical C-terminus is not a canonical TA…
#> 3 ROI_00019 terminus_noncanonical noncanonical N-terminus is not a canonical TA…
#> 4 ROI_00019 terminus_noncanonical noncanonical C-terminus is not a canonical TA…
```

`NTERM` and `CTERM` mark a canonical TALE terminal domain, one that can
be expected to do its usual job; any other terminus is coded `XXXXX`
(see
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md)).
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)
decides it from a search of each terminus with the protein profile of
its TALE domain, and `array_report.tsv` holds the evidence:

Code

``` r
readr::read_tsv(file.path(pxo86_raw_dir, "array_report.tsv"),
                show_col_types = FALSE) %>%
  filter(array_id %in% c("ROI_00001", "ROI_00019")) %>%
  select(array_id, nterm_aa_evalue, nterm_aa_cover, nterm_aa_domains,
         cterm_aa_evalue, cterm_aa_cover, cterm_aa_far_gap) %>%
  mutate(across(ends_with("_evalue"), ~ format(.x, digits = 2))) %>%
  knitr::kable()
```

| array_id  | nterm_aa_evalue | nterm_aa_cover | nterm_aa_domains | cterm_aa_evalue | cterm_aa_cover | cterm_aa_far_gap |
|:----------|:----------------|---------------:|-----------------:|:----------------|---------------:|-----------------:|
| ROI_00019 | 2.4e-129        |          0.844 |                2 | 5.9e-18         |          0.133 |              242 |
| ROI_00001 | 2.4e-129        |          0.844 |                2 | 7.4e-123        |          0.656 |               96 |

Each short terminus matches its profile strongly: these segments are
TALE termini. But a canonical terminus must cover at least 90% of the
profile (`terminus_min_cover`), and none of these does. The two 230-aa
N-termini align with the N-terminal profile everywhere except profile
positions 106-150, an internal deletion of about 45 residues, so the
match comes in two pieces (`nterm_aa_domains`) covering 84% of the
profile. The 183-aa C-terminus of `ROI_00001` aligns over its full
length with profile positions 1-183, and the 42-aa C-terminus of
`ROI_00019` with positions 1-37 over its first 37 residues. These
C-termini lack the far end of the profile, 96 and 242 positions
respectively (`cterm_aa_far_gap`). The 90% threshold comes from the
curated TALEs of `tale_annotations`: their canonical termini cover 93%
of the profile or more, the truncTALE termini 84% or less.

Function is what the threshold stands for. The C-terminal region of a
TALE carries the nuclear localisation signals and, at its far end, the
acidic activation domain that recruits the host’s transcription
machinery: both C-termini here stop well before it, consistent with the
loss of the activation domain that Ji et al. and Read et al. report for
truncTALEs. The N-terminal region carries the type III secretion signal
and, next to the repeats, the degenerate repeats that bind the thymine
preceding the target; which function, if any, the 45-residue deletion
affects is not something sequence alone can settle.

Both arrays are real TALEs, and truncTALEs have a function of their own:
they suppress the *Xa1* resistance. So `tales(x, sanitize = TRUE)`,
which removes arrays whose data are inconsistent or incomplete, keeps
them; `sanitize = "canonical"` removes them when an analysis needs
canonical TALEs only:

Code

``` r
pxo86_canonical <- suppressWarnings(tales(pxo86, sanitize = "canonical"))
setdiff(unique(pxo86$array_id), unique(pxo86_canonical$array_id))
#> [1] "ROI_00001" "ROI_00019"
```

## 2 Not the same kind of short

`array_report.tsv`’s `cterm_dna_hit` column already distinguishes the
two arrays, and the reason is worth tracing back to discovery:

Code

``` r
pxo86_report <- readr::read_tsv(
  file.path(pxo86_raw_dir, "array_report.tsv"), show_col_types = FALSE
)
pxo86_report %>%
  filter(array_id %in% c("ROI_00001", "ROI_00019")) %>%
  select(array_id, nterm_dna_hit, cterm_dna_hit, cterm_aa_length, orf_coverage)
#> # A tibble: 2 × 5
#>   array_id  nterm_dna_hit cterm_dna_hit cterm_aa_length orf_coverage
#>   <chr>     <lgl>         <lgl>                   <dbl>        <dbl>
#> 1 ROI_00019 TRUE          FALSE                      42           93
#> 2 ROI_00001 TRUE          TRUE                      183           83
```

[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)
finds TALE loci with three separate profile HMMs (one each for the
N-terminus, the repeat unit, and the C-terminus), and `nterm_dna_hit`
and `cterm_dna_hit` record whether the N- and C-terminus profiles
matched somewhere in an array’s merged hits, independently of how the
final protein turns out. `hits_report.tsv`, one of the three reports
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)
always writes (the other two are `domains_report.tsv` and
`array_report.tsv` itself), carries every individual hit with its
profile name and `frameshift_count`, and shows exactly what happened
downstream of each array’s last complete repeat:

Code

``` r
hits <- readr::read_tsv(file.path(pxo86_raw_dir, "hits_report.tsv"), show_col_types = FALSE)

hits %>%
  filter(array_id %in% c("ROI_00001", "ROI_00019")) %>%
  group_by(array_id) %>%
  summarise(
    cterm_hit_found = any(grepl("C-terminus", query_name)),
    cterm_hit_frameshift_count = frameshift_count[grepl("C-terminus", query_name)][1],
    .groups = "drop"
  )
#> # A tibble: 2 × 3
#>   array_id  cterm_hit_found cterm_hit_frameshift_count
#>   <chr>     <lgl>                                <dbl>
#> 1 ROI_00001 TRUE                                     2
#> 2 ROI_00019 FALSE                                   NA
```

`ROI_00001` has a real C-terminus hit, full length at ~286 codons, but
the nhmmer alignment needed two internal reframings to call it at all.
The repeat immediately before it is short and frameshifted the same way.
A C-terminus-shaped signal is there, just out of frame.

`ROI_00019` has no C-terminus hit. The hits that make up an array must
score above a threshold (`cterm_min_score`, 200 bits for the C-terminus
profile), set for finding arrays. Weaker matches next to the repeats are
not part of the array, but `array_report.tsv` still describes them:

Code

``` r
pxo86_report %>%
  filter(array_id %in% c("ROI_00001", "ROI_00019")) %>%
  select(array_id, cterm_dna_hit, cterm_dna_score, cterm_dna_cover, cterm_dna_pieces)
#> # A tibble: 2 × 5
#>   array_id  cterm_dna_hit cterm_dna_score cterm_dna_cover cterm_dna_pieces
#>   <chr>     <lgl>                   <dbl>           <dbl>            <dbl>
#> 1 ROI_00019 FALSE                    94.3           0.204                2
#> 2 ROI_00001 TRUE                   1074.            1                    1
```

Downstream of `ROI_00019`’s last repeat, the C-terminus profile matches
in two pieces, the best at 94 bits, covering together a fifth of the
profile; a complete C-terminus scores around 1100 bits. The raw nhmmer
output (`hmmer_search_out.txt`) places them at the two ends of the
profile, positions 1-97 and 783-861. The DNA thus keeps the start and
the end of a C-terminus, with about 685 nt between them that match
nothing, and the protein stops 42 residues into the first piece.

That difference carries forward into correction. Frameshift correction,
by either path, acts on the DNA span
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)
merges from these hits. `ROI_00001`’s span already contains a
full-length, C-terminus-shaped template, out of frame; `ROI_00019`’s
span holds no C-terminus hit, only these fragments beside it. Only one
array gives a correction tool a template to reframe into.

## 3 Does correction respect that difference?

### 3.1 `correct_array = TRUE`

Code

``` r
pxo86_decipher_dir <- file.path(out, "PXO86_decipher")
```

Code

``` r
invisible(tell_tales(
  subject_file = pxo86_fa, output_dir = pxo86_decipher_dir,
  correct_array = TRUE, max_comparisons = 20
))
#> Finding the closest reference amino acid sequences:
#> ================================================================================
#> 
#> Time difference of 6.23 secs
#> ================================================================================
#> 
#> Time difference of 38.26 secs
```

Code

``` r
decipher_report <- readr::read_tsv(
  file.path(pxo86_decipher_dir, "array_report.tsv"), show_col_types = FALSE
)
decipher_report %>%
  filter(array_id %in% c("ROI_00001", "ROI_00019")) %>%
  select(array_id, nterm_dna_hit, cterm_dna_hit, cterm_aa_length, orf_coverage)
#> # A tibble: 2 × 5
#>   array_id  nterm_dna_hit cterm_dna_hit cterm_aa_length orf_coverage
#>   <chr>     <lgl>         <lgl>                   <dbl>        <dbl>
#> 1 ROI_00019 TRUE          FALSE                      42           93
#> 2 ROI_00001 TRUE          TRUE                      216           86
```

`ROI_00001` is extended: its C-terminus grows from 183 to 216 aa, and
`orf_coverage` improves from 83% to 86%. `ROI_00019` is untouched. This
is [Section 2](#sec-mechanism)’s distinction playing out directly:
[`DECIPHER::CorrectFrameshifts()`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
reframes `ROI_00001` into the C-terminus-shaped template nhmmer already
found, and has no comparable template to work with for `ROI_00019`.

> **`max_comparisons` does not change this**
>
> `max_comparisons = 20` is used above to keep this article’s build time
> reasonable. The outcome is identical at every value tried, including
> the default (50) and the full reference set of 494 sequences (measured
> once, about 22 minutes, and not re-run here):
>
> | `max_comparisons` | `ROI_00001` C-terminus | `ROI_00019` C-terminus |
> |-------------------|------------------------|------------------------|
> | 20                | 216 aa                 | 42 aa                  |
> | 50 (default)      | 216 aa                 | 42 aa                  |
> | all 494           | 216 aa                 | 42 aa                  |
>
> Whatever makes
> [`DECIPHER::CorrectFrameshifts()`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
> extend `ROI_00001` is insensitive to how many reference sequences it
> is allowed to consider.

### 3.2 `correct_tales()`

Code

``` r
pxo86_java_fa <- file.path(out, "PXO86_java_corrected.fa")
```

Code

``` r
corrections <- correct_tales(
  uncorrected_path = pxo86_fa, corrected_path = pxo86_java_fa,
  return_corrections = TRUE
)
```

Genome-wide,
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md)
needed only 1 correction:

Code

``` r
corrections
#> # A tibble: 1 × 4
#>   seqName     posInOriginSeq type      substitution
#>   <chr>                <dbl> <chr>     <chr>       
#> 1 NZ_CP007166        2024199 insertion - -> g
```

Code

``` r
pxo86_java_dir <- file.path(out, "PXO86_java_tell")
```

Code

``` r
invisible(tell_tales(
  subject_file = pxo86_java_fa, output_dir = pxo86_java_dir, correct_array = FALSE
))
java_report <- readr::read_tsv(
  file.path(pxo86_java_dir, "array_report.tsv"), show_col_types = FALSE
)
java_report %>%
  filter(array_id %in% c("ROI_00001", "ROI_00019")) %>%
  select(array_id, nterm_dna_hit, cterm_dna_hit, cterm_aa_length, orf_coverage)
#> # A tibble: 2 × 5
#>   array_id  nterm_dna_hit cterm_dna_hit cterm_aa_length orf_coverage
#>   <chr>     <lgl>         <lgl>                   <dbl>        <dbl>
#> 1 ROI_00019 TRUE          FALSE                      42           93
#> 2 ROI_00001 TRUE          TRUE                      183           83
```

A single-base insertion, and it lands inside `ROI_00001`’s own span,
consistent with [Section 2](#sec-mechanism), since that is the array
with frameshifted hits for the correction to act on. It has no visible
effect: `ROI_00001`’s C-terminus stays at 183 aa and its `orf_coverage`
at 83%, against the 33 amino acids `correct_array = TRUE` adds.
`ROI_00019` is unaffected, as before.

## 4 Summary

| Array     | What it actually is                                       | correct_array = TRUE    | correct_tales() |
|:----------|:----------------------------------------------------------|:------------------------|:----------------|
| ROI_00019 | early stop, only fragments of a C-terminus downstream     | unchanged               | unchanged       |
| ROI_00001 | genuine frameshift, C-terminus hit present (out of frame) | over-corrected (+33 aa) | unchanged       |

Two genuine truncTALEs in the same genome respond differently to the
same correction tool, because they arose from different underlying
DNA-level events. In `ROI_00019`, the protein stops early and only
fragments of a C-terminus remain downstream, too weak to count as a hit,
so a frameshift corrector has no template to reframe into. A genuine,
evolved frameshift leaves exactly the kind of disrupted-but-recognisable
signal these tools exist to fix. `correct_array = TRUE` cannot tell that
signal from a sequencing error, and “fixes” it.

**On this evidence,
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md)
is the more transparent choice for a genome suspected of carrying real
truncTALEs.** It left both arrays’ sequence alone, where
`tell_tales(correct_array = TRUE)` rewrote the genuine frameshift’s
C-terminus as though it were an assembly error. One genome and two
arrays is a narrow base for a general rule. In practice: check the
terminus columns and `orf_coverage` in `array_report.tsv` before
correcting, and verify any change a correction makes to a documented or
suspected truncTALE’s sequence.

- A short C-terminus (`cterm_aa_cover` well below 1) with
  `cterm_dna_hit = FALSE`, only weak DNA fragments (`cterm_dna_cover`
  low) and `orf_coverage` in the normal range is very likely a genuine
  short protein. Any sequence a correction adds to it was never there.
- A short C-terminus with `cterm_dna_hit = TRUE`, a full-length DNA
  match (`cterm_dna_cover` near 1) but reduced `orf_coverage` may carry
  a genuine frameshift that is part of the strain’s biology.
  `correct_array = TRUE` will extend it regardless of whether that is
  what your analysis needs.

Neither correction path can make the underlying biological judgement
call for you.

## 5 Next

This article stands on its own. To continue the walkthrough it extends,
[return to mining TALE
sequences](https://scunnac.github.io/tantale/articles/tale_mining.md).
