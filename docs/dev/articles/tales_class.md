# The tales class

Every function in tantale that describes, compares or plots a set of
TALEs takes or returns a `tales` object. This article is a deep dive
into what that object is, how to build one, and what its most important
column – `dom_code` – actually means.

> **Where this fits**
>
> This article extends [Mining TALE
> sequences](https://scunnac.github.io/tantale/dev/articles/tale_mining.md):
> read it once that walkthrough has given you a `tales` object to look
> at, and return there afterwards to continue with classification. It is
> not itself a numbered step – its sibling deep dive, on the `tales_msa`
> class, comes later.

Code

``` r
library(tantale)
library(dplyr)
```

## 1 One row per part, not per TALE

A `tales` object is a long table with **one row per part** of a TALE
array: one row for the N-terminus, one for each central repeat, one for
the C-terminus. It is not one row per TALE.

That shape is forced by the biology. A TALE array is built from a
variable number of repeats – as few as a handful, more than thirty in
some natural alleles – so “one row per TALE” would need a ragged number
of columns, one per possible repeat position. A long table sidesteps
that entirely: an array with 12 repeats and one with 27 are both just
more rows, distinguished by `array_id`.

Building one from a real
[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
run looks like this. The example below uses the small sequence set
shipped with the package, so it runs in a few seconds without any
correction step:

Code

``` r
out <- file.path(tempdir(), "tales_class_article")
invisible(tell_tales(
  subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                             package = "tantale", mustWork = TRUE),
  output_dir = out))
```

Code

``` r
x <- suppressWarnings(tales_from_telltale(out))
x
#> <tales> 4 arrays, 96 parts
#>   layers: rvd   |   6 other columns
#>              rvd
#>   ROI_00001  NTERM NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN N ...
#>   ROI_00002  NTERM NN HD NI NN HD NG HD HD NG NG NI NG NI NG CTERM
#>   ROI_00003  NTERM NN ND NN NI NK NN HD NN NG NG N* HD N* HD NI NN HD NG H ...
#>   ROI_00004  NTERM NI HD NN NS NN NG HD NG HD NG NN NG HD NS HD NI NG HD H ...
```

Two columns locate a part within its array, and they answer different
questions:

- `position_in_array` numbers **every** part, termini included, starting
  at 1 (the N-terminus, when present).
- `position_in_crd` numbers only the **repeats**, starting at 1 for the
  first one, and is `NA` on the two termini.

They agree wherever there is no N-terminus to offset them, and differ by
exactly the number of non-repeat parts that come before a given position
otherwise. This is also why keeping only the repeats of an array
(`filter(x, domain_type == "repeat")`) is a perfectly valid `tales` and
not a broken one: `position_in_crd` stays contiguous from 1 even though
`position_in_array` no longer starts there. Completeness – every part
present, `position_in_array` running `1:n` – is a precondition some
functions need (alignment, chiefly, checked by
[`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md)),
not an invariant of the class itself.

## 2 The column contract

Two columns are required outright, because without them a row cannot be
located at all: `array_id` and `position_in_array`. A third requirement
is looser – at least one of `rvd` or `dom_code` must be present, since a
`tales` with neither has no residues to describe.

Everything else is optional, but validated for shape when it is there:

| column              | holds                                                                                                       |
|---------------------|-------------------------------------------------------------------------------------------------------------|
| `domain_type`       | `"N-terminus"`, `"repeat"`, or `"C-terminus"`                                                               |
| `position_in_crd`   | position among repeats only, `NA` on termini                                                                |
| `aa_seq`, `dna_seq` | the part’s amino acid / nucleotide sequence                                                                 |
| `seqnames`          | the source contig, for parts read off a genome                                                              |
| `dom_code`          | see below – minted by comparison, not by discovery                                                          |
| `source_directory`  | which [`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md) run produced this row |

A `tales` built straight from a fasta of RVD sequences (below) has none
of these beyond `rvd` itself, and that is a legitimate, if minimal,
object. Any further column – a strain name, a clade, a host range – is
preserved untouched and never warned about: metadata is expected to
accumulate on this table via ordinary
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) calls.

## 3 Three ways to build one

**From a
[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
run**, as above, with
[`tales_from_telltale()`](https://scunnac.github.io/tantale/dev/reference/tales_from_telltale.md).
This is the richest route: it carries sequences, source coordinates and
both position columns.

**From bare sequences**, with
[`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md)
– a fasta of RVD or repeat-code strings, sequence names becoming
`array_id`:

Code

``` r
rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
                         package = "tantale")
cat(readLines(rvd_fasta, n = 2), sep = "\n")
#> >TalA_BAI3
#> NN-NG-NN-HD-HD-NI-N*-NG-HD-NI-NG-NN-HD-NI-NG-NI-NG-NN-NG-HD-NI-NI-NG-HD-NN-NG

as_tales(rvd_fasta, sep = "-")
#> <tales> 11 arrays, 258 parts
#>   layers: rvd
#>                  rvd
#>   TalA_BAI3      NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN NG  ...
#>   TalA_CFBP1947  NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN NG  ...
#>   ...            ...
#>   TalA_MAI95     NN N* NN HD HD NI NG NN HD NS NG NI N* NN NG HD NI NI NG  ...
#>   TalA_MAI99     NN N* NN HD HD NI NG NN HD NS NG NI N* NN NG HD NI NI NG  ...
```

This result is deliberately column-poor – `array_id`,
`position_in_array` and `rvd`, nothing else – because that is all a bare
sequence file can supply. It is still a complete, valid `tales`.

**From a comparison**, with
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md)
(see below): the same object, given back with a `dom_code` column added.

All three produce the same class, checked by the same
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md),
usable with the same functions afterwards.

## 4 It is a tibble

A `tales` is a tibble with extra structure layered on top, not a
different kind of thing. Ordinary dplyr verbs work on it directly:

Code

``` r
x %>%
  filter(domain_type == "repeat") %>%
  count(array_id)
#> # A tibble: 4 × 2
#>   array_id      n
#>   <chr>     <int>
#> 1 ROI_00001    26
#> 2 ROI_00002    14
#> 3 ROI_00003    26
#> 4 ROI_00004    22
```

Subsetting that would break the class – dropping `array_id`, say –
quietly degrades the result to a plain tibble rather than erroring, so a
`tales` never lies about what it is. The mechanics of that degradation
are worth a closer look once `dom_code` is in the picture – see below.

[`print()`](https://rdrr.io/r/base/print.html) and
[`format()`](https://rdrr.io/r/base/format.html) follow the same idea
Biostrings uses for sequence sets: a header with the essentials (arrays,
parts, which residue columns are present), then the first and last
arrays rendered as sequences, eliding the middle.
[`summary()`](https://rdrr.io/r/base/summary.html) adds the slower
numbers – how redundant the domains are, how many repeats per array,
whether every array has both termini – and folds in
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md):

Code

``` r
summary(x)
#> <tales> summary
#>   arrays / parts            4 / 96
#>   distinct RVDs             8
#>   repeats per array         min 14   median 24   max 26
#>   arrays with both termini  4 of 4
#>   source sequences          2
#>   anomalies                 none
```

## 5 Anomalies are reported, not refused

Real TALE predictions are messy: a terminus can be duplicated,
misplaced, or simply not found.
[`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md)
loads such arrays anyway rather than refusing them – a class that
insists on clean input forces cleaning outside the package and throws
away exactly the signal a user would want to inspect.

“Anomalous” is deliberately a narrower category than “wrong”. An array
missing a terminus, for instance, is not flagged here at all: it may
genuinely sit at the edge of a contig. What *is* flagged is a part
arrangement that cannot be biologically real – two N-termini in one
array, say, or one that is not at the start:

Code

``` r
odd <- tibble::tibble(
  array_id = c("a1", "a1", "a1", "a2", "a2", "a2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L, 3L),
  domain_type = c("N-terminus", "N-terminus", "repeat",
                  "N-terminus", "repeat", "C-terminus"),
  rvd = c("NTERM", "NTERM", "HD", "NTERM", "HD", "CTERM")
)
x_odd <- suppressWarnings(tales(odd))
tales_anomalies(x_odd)
#> # A tibble: 2 × 3
#>   array_id check               detail                      
#>   <chr>    <chr>               <chr>                       
#> 1 a1       terminus_duplicated more than one N-terminus    
#> 2 a1       terminus_misplaced  N-terminus not at position 1
```

`tales(sanitize = TRUE)` drops the flagged arrays instead of merely
warning about them:

Code

``` r
tales(odd, sanitize = TRUE)
#> Warning: Dropped 1 array with biological anomalies.
#> ✖ Array: "a1"
#> ℹ Reasons: terminus_duplicated and terminus_misplaced
#> <tales> 1 array, 3 parts
#>   layers: rvd   |   1 other column
#>       rvd
#>   a2  NTERM HD CTERM
```

Structural problems – a duplicated key, an `NA` where one is not allowed
– are a different matter and remain hard errors either way: the table
cannot be interpreted at all, rather than merely describing something
odd.

## 6 Projections: turning parts back into sequences

Some consumers want a `tales` as a table; others want each array back as
a single string. Three functions bridge the gap:

- [`tales_rvd_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md)
  – one RVD string per array, the form target-prediction tools consume.
  Drops the termini by default (`rvd_only = TRUE`), since prediction
  concerns the repeat domain only.
- [`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md)
  – the `dom_code` sibling, for ARLEM and MAFFT’s text mode. Keeps the
  termini by default, since alignment wants them as anchors, and
  separates codes with a space rather than a dash, because that is what
  those tools split on.
- [`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md)
  – not a per-array string but the lookup table itself: one row per
  distinct `dom_code`, with the amino acid sequence and RVD it stands
  for.

Code

``` r
xa <- tales_assign_domain_codes(x)
tales_coded_strings(xa)["ROI_00002"]
#> BStringSet object of length 1:
#>     width seq                                               names               
#> [1]    44 40 19 9 16 36 24 1 7 24 12 26 28 26 16 11 46      ROI_00002
tales_domain_codes(xa) %>% head(4)
#> # A tibble: 4 × 3
#>   dom_code aa_seq                             rvd  
#>   <chr>    <chr>                              <chr>
#> 1 1        LPPDQVVAIASNGGGKQALETVQRLLPVLCQAHG NG   
#> 2 10       LTPAQVVAIASNDGGKQALETVQRLLPVLCQAHG ND   
#> 3 11       LTPAQVVAIASNGGGKQALE               NG   
#> 4 12       LTPAQVVAIASNGGGKQALETVQRLLPVLCQAHG NG
```

## 7 What a `dom_code` is

This is the column the rest of the comparison machinery is built on, and
it is worth being precise about.

### 7.1 It names a domain, not a repeat

A TALE part is an N-terminus, a repeat, or a C-terminus, and **all three
get `dom_code`s on the same footing**. “Domain” is the word used here
specifically because it has to cover a repeat and a terminus
indiscriminately – “repeat code” would be the wrong name for a column
that also numbers C-termini. Two parts with identical amino acid
sequences get the same code, whatever kind of part they are.

The distinction is not pedantic. In the small example used throughout
this article:

Code

``` r
n_distinct(xa$dom_code)
#> [1] 47
table(xa$domain_type[match(unique(xa$dom_code), xa$dom_code)])
#> 
#> C-terminus N-terminus     repeat 
#>          4          4         39
```

47 distinct codes, and only 39 of them are repeats – the rest are
distinct N- and C-termini.

> **Not a repeat count**
>
> Describing the total as a repeat count would overstate the repeat
> diversity by a substantial margin – this exact error was made, and
> caught, while writing
> [`summary.tales()`](https://scunnac.github.io/tantale/dev/reference/summary.tales.md).

### 7.2 It is what makes a TALE alignable

Aligning TALE arrays residue by residue is meaningless: the repeats are
near-identical (~34 amino acids, differing mainly at two positions), so
almost anything matches almost anything else. Giving each distinct
repeat a symbol turns an array into a *sequence of repeat units*
instead, which can be aligned the way a protein sequence is – with
insertions and deletions of whole repeats rather than of individual
residues. This is why
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md)
works on `dom_code` (or `rvd`), never on `aa_seq` directly.

### 7.3 It is finer than an RVD

The RVD is residues 12-13 of a repeat, and it is what determines which
DNA base the repeat binds. Two repeats can share an RVD – the same
binding specificity – while differing elsewhere in the protein, and they
get **different** `dom_code`s. So `rvd` is the functional layer and
`dom_code` the identity layer, and the class carries both because they
answer different questions.

In this example, the RVD `HD` (the code for cytosine) is shared by
several distinct repeat proteins:

Code

``` r
reps <- xa %>% filter(domain_type == "repeat")
reps %>% filter(rvd == "HD") %>% distinct(dom_code, aa_seq) %>% head(4)
#> # A tibble: 4 × 2
#>   dom_code aa_seq                            
#>   <chr>    <chr>                             
#> 1 24       LTPDQVVAIASHDGGKQALETVQRLLPVLCQAHG
#> 2 32       LTPEQVVAIASHDGGKQALETVQRLLPVLCQAHG
#> 3 7        LTPAQVVAIASHDGGKQALETVQRLLPVLCQAHG
#> 4 8        LTPAQVVAIASHDGGKQALETVQRLLPVLCQVHG
```

A terminus, meanwhile, has no RVD at all: its `rvd` column holds a
placeholder – `"NTERM"`, `"CTERM"`, or `"XXXXX"` for a terminus whose
CDS was detected but not identified (see
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md))
– while its `dom_code` is a real identifier of a real sequence. Another
reason the two layers are not interchangeable.

### 7.4 How they are computed

Plainly: group the parts by `aa_seq`, number the groups.
[`dplyr::cur_group_id()`](https://dplyr.tidyverse.org/reference/context.html),
nothing cleverer, done by
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md).

### 7.5 The numbers are only meaningful within one call

Because the numbering depends on which distinct sequences are present,
code `42` from one comparison is not code `42` from another – add an
array, remove one, or run it on a different day’s data, and the same
protein can land on a different number. This is not a wart to work
around; it is a fact to know before using `dom_code` for anything, and
the package enforces it rather than merely documenting it.

Every object
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md)
produces carries a namespace stamp, readable with
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
derived from the exact set of sequences that produced it:

Code

``` r
tales_namespace(xa)
#> [1] "7377e3f80aa36898"
```

Mixing codes from two different calls is refused rather than silently
computing a wrong answer. The rest of this section demonstrates that,
one step at a time.

Domain distances over the whole example set are the baseline to mix in
something incompatible with:

Code

``` r
dd_all <- tales_domain_distances(xa)
```

Splitting the example array set in two and re-assigning codes to one
part gives it a *different* namespace, even though every sequence it
contains was already present in the full set:

Code

``` r
sub <- xa[xa$array_id != "ROI_00001", ][setdiff(names(xa), "dom_code")]
xa_sub <- tales_assign_domain_codes(sub)

tales_namespace(xa_sub) == tales_namespace(xa)
#> [1] FALSE
```

> **Mixing namespaces is a caught error, not a silent one**
>
> Using `xa_sub`’s codes against `dd_all` – distances keyed by the other
> namespace – does not compute a plausible-looking wrong answer. It
> refuses outright:
>
> Code
>
> ``` r
> tales_tale_distances(xa_sub, dd_all)
> #> Error in `.assert_same_namespace()`:
> #> ! `x` and `domain_distances` come from different runs.
> #> ✖ Namespaces "39eacafcbbcbfa90" and "7377e3f80aa36898".
> #> ℹ Domain codes are only meaningful within the call that minted them.
> #> ℹ Recompute the distances from this `x` with `tales_domain_distances()`.
> ```

## 8 Where `dom_code` leads: comparing TALEs

Minting codes is the first of three steps that together compare a set of
TALEs –
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md)
first, then
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md)
(how different the distinct domains are from each other, pairwise), then
[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md)
(how different the *arrays* are, aligning their domain-code sequences
using the domain distances as the substitution cost).
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
runs all three and returns the corresponding `tales`, `domain_distances`
and `tale_distances` objects together:

Code

``` r
cmp <- tales_compare(x)
names(cmp)
#> [1] "tales"            "domain_distances" "tale_distances"
```

Each step is independently useful – the domain-level distances alone
answer questions about repeat diversity, without ever running the
array-level alignment – and each is documented in its own right; see
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
for the full chain.

## 9 Subsetting keeps the class only while the contract holds

A `tales` is a tibble, so every subsetting route a tibble supports works
on it too: `[`,
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html),
[`select()`](https://dplyr.tidyverse.org/reference/select.html),
[`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html),
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html). What
differs from a plain tibble is what happens to the *class* along the
way.

**Row subsetting never breaks anything.** Every invariant a `tales`
checks is closed under keeping a subset of rows – which is exactly why
completeness ([Section 1](#sec-one-row-per-part), above) is a
*precondition* of specific functions rather than part of the class
itself. Filtering to one array, or to its repeats only, is still a valid
`tales`:

Code

``` r
one_array <- xa[xa$array_id == "ROI_00002", ]
is_tales(one_array)
#> [1] TRUE
nrow(one_array)
#> [1] 16
```

**Column subsetting is where it gets interesting.** Dropping a column
the contract needs – `array_id`, or both residue columns at once –
leaves something that can no longer be described as a `tales`, and the
class says so by quietly stepping out of the way rather than by
erroring:

Code

``` r
no_id <- xa %>% select(-array_id)
is_tales(no_id)
#> [1] FALSE
class(no_id)
#> [1] "tbl_df"     "tbl"        "data.frame"
```

The same happens whichever route removes the load-bearing column –
[`select()`](https://dplyr.tidyverse.org/reference/select.html),
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html), or
plain `[, j]`. This is deliberate: a `tales` object that kept claiming
to be one after losing its key would let every downstream function trust
a promise the object could no longer keep. Dropping an optional column
(`seqnames`, `aa_seq`, …) has no such consequence, and the class
survives untouched.

[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) earns
one extra check the others do not need: it can overwrite `array_id` or
`position_in_array` *in place*, which plain subsetting cannot do, so the
result could keep every contract column while no longer being a valid
key. That case is caught separately, immediately after the mutate:

Code

``` r
xa %>% mutate(position_in_array = 1L)
#> Error in `.tales_check_key()`:
#> ! array_id and position_in_array must together be unique.
#> ✖ 92 duplicated rows in 4 arrays: "ROI_00001", "ROI_00002", "ROI_00003", and
#>   "ROI_00004"
```

What is **not** re-checked on every verb is the row-level biology –
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md)
does not run again on each
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html) call,
because anomalies are a property to inspect deliberately (see above),
not a cost to pay on every pipe.

**Attributes travel with a valid subset**, which matters most for
`dom_code_namespace`. Cutting `xa` down to one array with `[` keeps its
namespace, because the namespace describes *the run the codes came
from*, not which rows happen to be kept right now:

Code

``` r
tales_namespace(one_array) == tales_namespace(xa)
#> [1] TRUE
```

That is only true because `one_array` was cut with `[`, not re-coded –
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md)
on a subset mints a *different* namespace, as the mismatch example above
already showed.

A `tales_msa` carries one further attribute this way,
`alignment_position`, and losing *that* while keeping everything else
degrades it by exactly one step – to a plain `tales`, not all the way to
a tibble. [The `tales_msa` class
article](https://scunnac.github.io/tantale/dev/articles/tales_msa_class.md)
covers that graded degradation in full, together with the reverse
direction: how a `tales` is promoted into an alignment in the first
place.

## 10 Next

This article stands on its own; where to go next depends on what you
were after. To continue the walkthrough this extends, [return to
classifying TALE
sequences](https://scunnac.github.io/tantale/dev/articles/tale_classification.md).
To go deeper on alignment, [the `tales_msa`
class](https://scunnac.github.io/tantale/dev/articles/tales_msa_class.md)
picks up directly from [Section 7](#sec-dom-code) and
[Section 9](#sec-subsetting), above.
