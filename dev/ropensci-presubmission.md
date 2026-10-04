Submitting Author Name: Sébastien Cunnac
Submitting Author Github Handle: <!--author1-->@scunnac<!--end-author1-->
Other Package Authors Github handles: (comma separated, delete if none) <!--author-others-->@vibaotram<!--end-author-others-->
Repository:  <!--repourl-->https://github.com/scunnac/tantale<!--end-repourl-->
Submission type: <!--submission-type-->Pre-submission<!--end-submission-type-->
Language: <!--language-->en<!--end-language-->

---

-   Paste the full DESCRIPTION file inside a code block below:

```
PASTE THE DESCRIPTION FILE HERE AT POSTING TIME
```


## Scope 

- Please indicate which category or categories from our [package fit policies](https://ropensci.github.io/dev_guide/policies.html#package-categories) or [statistical package categories](https://stats-devguide.ropensci.org/overview.html#overview-categories) this package falls under. (Please check one or more appropriate boxes below):

    **Data Lifecycle Packages**

	- [ ] data retrieval
	- [ ] data extraction
	- [x] data munging
	- [ ] data deposition
    - [ ] data validation and testing
	- [ ] workflow automation
	- [ ] version control
	- [ ] citation management and bibliometrics
	- [x] scientific software wrappers
	- [ ] field and lab reproducibility tools
	- [ ] database software bindings
	- [ ] geospatial data
	- [ ] translation
    
     **Statistical Packages**

	- [ ] Bayesian and Monte Carlo Routines
	- [ ] Dimensionality Reduction, Clustering, and Unsupervised Learning
	- [ ] Machine Learning
	- [ ] Regression and Supervised Learning
	- [ ] Exploratory Data Analysis (EDA) and Summary Statistics
	- [ ] Spatial Analyses
	- [ ] Time Series Analyses
	- [ ] Probability Distributions


- Explain how and why the package falls under these categories (briefly, 1-2 sentences).  Please note any areas you are unsure of:

  tantale is a toolkit for the analysis of TAL effectors (TALEs) of
  *Xanthomonas*. Part of it is its own analysis code: `tell_tales()` finds
  TALE genes in genome assemblies and corrects frameshifts in error-prone
  long-read ones; `tales_align()` aligns TALEs repeat by repeat;
  `tales_compare_distal()` and `tales_compare_functal()` reimplement two
  published TALE distances; `tales_group_hclust()` and
  `tales_group_kmedoids()` group TALEs into families. These rest on S3
  classes (`tales`, `tales_msa`, `pairwise_distances`) with methods to
  validate, subset, combine and plot them. The other part wraps the
  programs of the field, AnnoTALE, PrediTALE and TALEcorrection (Java),
  TALVEZ (Perl), and MAFFT, HMMER and MMseqs2 (in a conda environment
  that `tantale_setup()` builds), and reads their results into the same
  objects.

  We ticked the two categories that fit the wrappers and the classes. We
  are unsure how the analysis code is judged: it is specific to TALEs and
  serves the same workflow, but it is more than a wrapper, and the
  package also has plotting methods (`plot.tales()`, `plot.tales_msa()`,
  `talomes_heatmap()`). Does this fit the scope?

- If submitting a statistical package, have you already [incorporated documentation of standards into your code via the **srr** package](https://stats-devguide.ropensci.org/pkgdev.html#pkgdev-srr)?

  Not applicable.

-   Who is the target audience and what are scientific applications of this package?  

  Plant pathologists and microbiologists who work on *Xanthomonas*, the
  bacteria behind bacterial blight and leaf streak of rice, citrus canker
  and many other crop diseases. TALEs are injected into plant cells, where
  they bind promoters through a DNA-binding domain of near-identical repeats
  and switch on host genes; the repeats make TALE genes hard to assemble
  and to annotate. Applications: inventorying the TALE repertoire of newly
  sequenced strains (including Oxford Nanopore assemblies), classifying
  TALEs into families across strains, and predicting their target genes in
  the host genome.

-   Are there other R packages that accomplish the same thing? If so, how does yours differ or meet [our criteria for best-in-category](https://ropensci.github.io/dev_guide/policies.html#overlap)?

  We know of no R package for TALE analysis. The field relies on
  stand-alone programs and web servers (AnnoTALE, the QueTAL suite,
  PrediTALE, TALVEZ), each with its own input and output formats. tantale
  wraps them or reimplements them in R, adds TALE discovery in error-prone
  assemblies and repeat-level alignment, and keeps the results in one set
  of objects, so that a whole analysis runs in one R session.

-   (If applicable) Does your package comply with our [guidance around _Ethics, Data Privacy and Human Subjects Research_](https://devguide.ropensci.org/policies.html#ethics-data-privacy-and-human-subjects-research)?

  Not applicable: the package handles bacterial and plant sequence data.

-  Any other questions or issues we should be aware of?:

  1. **External software.** The main workflow needs a Java runtime and a
     conda environment with pinned versions (MAFFT 7.453, HMMER 3.3.2;
     later MAFFT versions change the alignment of repeat strings). The
     Java tools and four example genomes are downloaded from this
     repository's GitHub releases by `tantale_setup(install = TRUE)` into
     `tools::R_user_dir()`; the source package itself is about 1 MB.
  2. **Tests.** Tests that need these programs fail when they are missing,
     because we prefer a failing test to one that silently does not run.
     Our GitHub Actions workflow installs them first. How should this fit
     with the pkgcheck run of the review bot?
  3. **Platforms.** Linux is tested. The pinned bioconda builds also exist
     for macOS x86_64, which we have not yet tested; there are none for
     Windows or Apple Silicon, so the package declares `OS_type: unix`.
  4. **Redistribution.** The wrapped programs keep their own licences
     (GPL-3 for the Jstacs tools), listed in `inst/COPYRIGHTS`. TALVEZ and
     the QueTAL code are redistributed with their authors' permission.
  5. **Dependencies.** 40 packages in Imports, about 15 from Bioconductor.
     Each has a call site in the package.

## Use of Generative AI

- [x] Generative AI tools were used to produce some of the material in this submission.

If so, please describe usage, and include links to any relevant aspects of your repository. See [our blog post](https://ropensci.org/blog/2026/02/26/ropensci-ai-policy/) for background. (Explicit advice is not yet included in our _Dev Guide_; we are hoping to update very soon, and ask your cooperation and transparency in the meantime.)

  The 2026 restructuring
  of the package (refactoring, tests, documentation and website articles)
  was done with Claude (Anthropic) through Claude Code, working under the
  maintainer's direction; every plan was agreed and every change reviewed
  by the maintainer. The decisions and their reasons are recorded in
  [`dev/restructuring-notes.md`](https://github.com/scunnac/tantale/blob/main/dev/restructuring-notes.md),
  and the working instructions given to the assistant in
  [`dev/CLAUDE.md`](https://github.com/scunnac/tantale/blob/main/dev/CLAUDE.md).
