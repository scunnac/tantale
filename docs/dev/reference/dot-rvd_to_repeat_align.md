# View an RVD alignment in terms of repeat codes

The inverse direction of `repeat_to_rvd_align()`: given an alignment
computed on RVDs, substitute each non-gap cell with the corresponding
repeat code. The two alphabets differ greatly – a few dozen RVDs against
hundreds of mostly-singleton repeat codes – so aligning on one and
viewing as the other is a genuinely different result from aligning on
the other directly.

## Usage

``` r
.rvd_to_repeat_align(rvd_msa_by_group, repeat_vecs)
```

## Arguments

- rvd_msa_by_group:

  A character matrix of aligned RVDs, rows named by array.

- repeat_vecs:

  A named list of repeat-code vectors, one per row of
  `rvd_msa_by_group`.

## Value

A character matrix with the dimensions and dimnames of
`rvd_msa_by_group`.

## Details

The mapping is **positional**: the k-th non-gap cell of a row is taken
to be the k-th element of that row's repeat vector. That is only well
defined when the two agree in length, which is now checked.
