
# Minor thinkgs that need a fix

- tales_predict_targets(): depending on the method argument, the order of the column in the returned object is not the same. This comes from inconsistency at the level of the individual prediction functions (see their exemple). Those should be updated to have a matching column order. Tool specific columns being last. A bonnus of the function could be to enable running both methods in one call and return a combinded results table.

- In tales_consensus() documentation indicate how one can obtain the aling object from a tales_msa (ie as.matrix).

- the exemple plot in plot_target_preds is awful. `filter_range` should be a lot shorter.

- talomes_heatmap()'s default value for `group_col` should be "group"

- review the documentation for tantale-package. Some facts are obsolete or clumsily stated.

- This kind of wording 'Each TALE given is placed' (the passive form) is found too often in the website and package documentation. When found, a decision must be made as to whether a direct wording is not more desirable.

- unless I am wrong, the various obects (DNA, or protein sequnces, rvd sequence, parts) that the `fasta_file` argument can take are not garanteed to work by a specific test. This may need to be fixed. Also the exact nature of what can be passed to this argument needs more details in the doc (what type of objects concretely?).


# Questions/ideas

- Rather that a msa, can one extract from maftt distance between sequence? That could be an alternative to ARLEM.

- Shouldn't a tales_msa hold a reference to the actual residues layer (rvd or dom code and distance object reference) used to compute it?

- talvez() is currently to take custom rdv <-> nt association matrices. A mecanisms enabling users to suply could be used introduced in a subsequent release of the package as a new feature.

- if nothing is supplied to do something usefull with the extended class builder written by run_annotale_assign. This function is of limited value if this object carries more information than what is aleady covered by the written tables. Even though a bit cumbersome, a parser function could be created to load this info into R. The desirability and the details of the implementaton need to be discussed.

- shouldn't the default values of the `java_args` argument be homogenized across the java functions of tantale?

