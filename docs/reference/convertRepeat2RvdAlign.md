# Substitute Distal repeat IDs for RVDs in a TALE alignment matrix.

Substitute Distal repeat IDs for RVDs in a TALE alignment matrix.

## Usage

``` r
convertRepeat2RvdAlign(repeatAlign, repeat2RvdMapping)
```

## Arguments

- repeatAlign:

  A multiple TALE repeat sequences alignment in the form of a matrix as
  returned by
  [`buildRepeatMsa`](https://scunnac.github.io/tantale/reference/buildRepeatMsa.md)
  or the `SeqOfRepsAlignments` element in the return object of the
  [`buildDisTalGroups`](https://scunnac.github.io/tantale/reference/buildDisTalGroups.md)
  function.

- repeat2RvdMapping:

  The return value of the
  [`getRepeat2RvdMapping`](https://scunnac.github.io/tantale/reference/getRepeat2RvdMapping.md)
  function or the
  [`getRepeat2RvdMappingFromDistalr`](https://scunnac.github.io/tantale/reference/getRepeat2RvdMappingFromDistalr.md)
  function if you used the
  [`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)
  function.

## Value

A TALE alignment matrix made up of RVD sequences.
