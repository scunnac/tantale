# Derive aa_seq from dna_seq by translation

[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
needs protein sequences, but an object may carry only the DNA. TALE part
coding sequences are in frame, so translating them recovers `aa_seq`
exactly.

## Usage

``` r
.translate_parts(dna)
```

## Arguments

- dna:

  A character vector of in-frame coding sequences.

## Value

A character vector of amino-acid sequences.

## Details

Two details are easy to get silently wrong. `no.init.codon = TRUE` is
required: TALE repeats begin on `CTG`/`TTG`, which are alternative start
codons, so the default forces the first residue to `M` – that alone
accounted for 865 of 955 mismatches on the reference fixture. And the
C-terminal parts carry a trailing stop codon, which is stripped.

With both handled, translation reproduces the stored `aa_seq` for all
955 parts of the reference fixture.
