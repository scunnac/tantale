
# Suggested modifications

tales_anomalies() returned df should be sorted on `array_id` and `check`

plots.tales() the legend for the 'Part and length' section is not what I had written initially. Why is the N-teminus, repeat, ... prefix are repeated when the outline provide the domain type.




# Findings

Now that we have debugged the way termini codes were added to the RVD sequences, it turns out that even a `max_comparisons=5` is enough to correct all BAI3-1-1 arrays. However, it seems that one array just disappear... What happens to the the 8th one? we should look into the array_report.tsv.
We need to revise the documentation, and our preconceived ideas about this parameter. You could rerun the test you once ran to try and define adapted values but this time the read out is the 'new' rvd_string column and the number of arrays both in the sanitize and non sanitized output.
I tried to re-frame the tale mining article in this direction but I wonder if I nailed it correctly. An important consequence of this is that now tell_tales in correction mode runs pretty fast, may be even faster that the java implementation for similar of even best results. Can you try to get a sense of that using the BAI3-1-1 genome?

This:
```
#' @param max_comparisons How many references each array may be aligned
#'   against. \code{NULL} means all of them.
```
should be changed to the values found 'optimal' best results at the fastest speed.

# For discussion

  - The more I think about it the more I believe the cut off alignment length with the termini protein sequences with the hmm should be increased, potentially to the length of the profile + o - a small margin. Biologically a termini that do align with the expected profile indicates that it is related to it. But biologicall, it is pretty unlikely that a severely truncated terminus is going to perform its biological function. Hence, we should assign a N or CTERM tag only if this promize can reasonably be realized and give a XXXXX otherwise. What is you opinion on that?


  - At some point we adressed a question regarding the sorting of the arrays in an object or a plot. I cannot remember. Now, I wonder if it should be done alphabetically on array_id when ploting `tales`. What do you think?


# Note

I have added 'comments' for you in the qmd file. The comments can be found because the line starts with @CLAUDE and are single liners.
