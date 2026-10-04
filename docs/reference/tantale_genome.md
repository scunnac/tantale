# Path to one of the example genomes

Four *Xanthomonas oryzae* pv. *oryzae* genome assemblies are used
throughout the articles and examples. They are not part of the package
itself:
[`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
downloads them, with the Java tools, when called with `install = TRUE`.
This function returns where one of them is.

## Usage

``` r
tantale_genome(strain = c("MAI1", "BAI3", "BAI3-1-1", "PXO86"))
```

## Arguments

- strain:

  Which genome.

## Value

The path of a FASTA file.

## Details

MAI1 (GenBank CP025609.1), BAI3 (CP025610.1) and PXO86 (RefSeq
NZ_CP007166.1) are complete genomes of African and Asian strains; the
sequences are those of the records, with shortened FASTA headers.
BAI3-1-1 is an unpublished assembly of a BAI3 derivative that carries
sequencing and assembly errors in its TALE loci. It was produced by the
authors of tantale, is available only from the package's `genomes-1`
release, and serves to illustrate frameshift correction.

## See also

[`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md),
which downloads them.

## Examples

``` r
# \donttest{
# Needs the genomes downloaded by tantale_setup(install = TRUE).
tantale_genome("MAI1")
#> [1] "/home/cunnac/snap/codium/495/.local/share/R/tantale/genomes-1/MAI1.fa"
# }
```
