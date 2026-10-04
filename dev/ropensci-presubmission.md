Submitting Author Name: Sébastien Cunnac
Submitting Author Github Handle: <!--author1-->@scunnac<!--end-author1-->
Other Package Authors Github handles: (comma separated, delete if none) <!--author-others-->@vibaotram<!--end-author-others-->
Repository:  <!--repourl-->https://github.com/scunnac/tantale<!--end-repourl-->
Submission type: <!--submission-type-->Pre-submission<!--end-submission-type-->
Language: <!--language-->en<!--end-language-->

---

-   Paste the full DESCRIPTION file inside a code block below:

```
Package: tantale
Title: Mining and Analysis of Transcription Activator-Like Effectors
Version: 0.99.0
Authors@R: c(
    person(given = c("Bao", "Tram"), family ="Vi",
    email = "vbt576@gmail.com", role = c("aut"),
    comment = c(ORCID = "0000-0002-4319-5544")),
    person(given = "Sebastien", family = "Cunnac",
    email = "sebastien.cunnac@ird.fr", role = c("aut", "cre"),
    comment = c(ORCID = "0000-0002-3695-491X")),
    person(given = c("Alvaro", "L."), family = "Pérez-Quintero",
    role = "cph",
    comment = "Author of the bundled TALVEZ 3.2 and of the FuncTAL table behind rvd_dna_specificity"),
    person(given = "Molly", family = "Megraw", role = "cph",
    comment = "Author of the PlantTFBS Java classes bundled with TALVEZ (simplescancode/)"),
    person(given = c("Artemis", "G."), family = "Hatzigeorgiou", role = "cph",
    comment = "Author of the PlantTFBS Java classes bundled with TALVEZ (simplescancode/)")
          )
Maintainer: Sebastien Cunnac <sebastien.cunnac@ird.fr>
Description: A toolkit for the study of transcription activator-like
    effectors (TALEs) of Xanthomonas. tell_tales() finds TALE genes in
    genome assemblies and can correct the frameshifts that make them hard
    to annotate in error-prone long-read assemblies. TALEs are held in S3
    classes that record each one as a series of parts (N-terminus, repeats
    with their repeat-variable diresidues (RVDs), C-terminus), with methods
    to validate, subset, combine and plot them. TALEs are aligned repeat by
    repeat, compared with R implementations of the 'DisTAL' and 'FuncTAL'
    distances <doi:10.3389/fpls.2015.00545>, and grouped into families by
    hierarchical clustering or k-medoids. Wrappers for 'AnnoTALE'
    <doi:10.1038/srep21077>, 'TALEcorrection'
    <doi:10.1186/s12864-023-09228-1>, 'TALVEZ'
    <doi:10.1371/journal.pone.0068464> and 'PrediTALE'
    <doi:10.1371/journal.pcbi.1007206> read their results into the same
    objects; the last two predict target sites in host genomes.
License: MIT + file LICENSE
Copyright: tantale authors, except the bundled third-party programs and
    data listed with their licences in inst/COPYRIGHTS
URL: https://scunnac.github.io/tantale/, https://github.com/scunnac/tantale
BugReports: https://github.com/scunnac/tantale/issues
Depends:
    R (>= 4.0.0)
SystemRequirements: Java (>= 8), Perl, and conda, mamba or micromamba
    (which tantale_setup() uses to install MAFFT, HMMER and MMseqs2)
Imports:
    magrittr, fs, dplyr, ggplot2 (>= 3.5.0), tibble, stringr, readr, glue,
    IRanges, Biostrings, pwalign, GenomicRanges, BSgenome, plyranges,
    systemPipeR, DECIPHER,
    ape, ggtree, tidytree, aplot, universalmotif,
    gplots, cluster, ggnewscale, matrixStats,
    reticulate, cli, digest, methods, rlang, tidyr, stats, utils, grDevices, graphics,
    BiocGenerics, BiocParallel, GenomeInfoDb, S4Vectors, rtracklayer
Suggests:
    knitr, DT, sessioninfo,
    biomartr,
    rmarkdown,
    quarto,
    covr,
    testthat (>= 3.0.0), withr, codetools
VignetteBuilder: quarto
Config/Needs/website:
    pkgdown,
    tidyverse/tidytemplate
Config/testthat/edition: 3
Config/testthat/parallel: false
Encoding: UTF-8
OS_type: unix
biocViews: Software
Config/roxygen2/version: 8.0.0
Roxygen: list(markdown = TRUE)
LazyData: true
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

  tantale wraps the programs the TALE field relies on (AnnoTALE,
  PrediTALE, TALEcorrection, TALVEZ, MAFFT, HMMER, MMseqs2), installing
  them and reading their output into validated S3 classes.

  Unsure: it also holds TALE-specific analysis code of its own (finding
  TALE genes in error-prone assemblies, repeat-level alignment, two
  published distance measures, clustering into families) and plotting
  methods for its classes. Does that fit the scope?

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

  1. **External software.** The workflow needs a Java runtime and a conda
     environment with pinned versions (MAFFT 7.453, HMMER 3.3.2; later
     MAFFT versions align TALE repeat strings differently). The Java
     tools and the example genomes are downloaded from this repository's
     releases by `tantale_setup(install = TRUE)`; the source package is
     about 1 MB. Only Unix is supported, since the pinned bioconda builds
     exist for Linux and Intel macOS alone.
  2. **Tests.** Tests needing those programs fail rather than skip when
     they are missing, since a check that silently does not run seems
     worse than none. Our GitHub Actions workflow installs them first.
     How should that work with the review bot's pkgcheck run?

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
