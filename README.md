<p align="right">
  <img src="./man/figures/tantale_logo_small.gif">

## An integrated collection of functions for [TALE](https://en.wikipedia.org/wiki/Transcription_activator-like_effector) minning and analysis with the R language


Because there are so many concurrent and complementary tools running on different platforms (java, perl), analyzing TALEs in (mostly Xanthomonas) genomes can turn into a nightmarish experience. Furthermore, making sense of the output is difficult because there is no easy way to graphically represent the various objects of the analysis.

With `tantale`, we compiled and extended our previous code wrapping TALE analysis tools into an integrated R interface that further provides an extensive list of utilities for easy plotting. This enables a moderately proficient R programmer to perform entire analysis pipelines directly in R and access result objects for custom manipulations.

Here is a snapshot of the topics that are or will (hopefully) be covered in the future:


- A TALE-oriented OOP framework:
    - A TALE class and associated methods (to be done)


- TALE mining in bacterial sequences:
    - Wrapper around [AnnoTALE](https://doi.org/10.1038/srep21077) and [correcTALE](https://doi.org/10.1186/s12864-023-09228-1)
    - tellTale, an R function similar to AnnoTALE
    - Analysis tools for RVD inventory, repeat lenght


- TALEs classification, phylogeny:
    - Wrappers around [DisTAL](https://doi.org/10.3389/fpls.2015.00545), [FuncTAL](https://doi.org/10.1038/srep21077), AnnoTALE
    - TALE groups inference
    - Easily build Multiple alignments and generate nice plots


- TALE targets predictions:
    - Wrappers around target predictors ([Talvez](https://doi.org/10.1371/journal.pone.0068464) and [PrediTALE](https://doi.org/10.1371/journal.pcbi.1007206))
    - General parser for results aggregation
    - Connector with [daTALbase](https://doi.org/10.1094/MPMI-06-17-0153-FI) (to be done)



**NOTE** :

- For further details, please, take a look at the package [page](https://scunnac.github.io/tantale)

- Install via the `remotes` package

```
remotes::install_github("cunnac/tantale",
                          type = "source",
                          dependencies = TRUE,
                          upgrade = "never"
                          )
```

- This is still a **work in progress** and is not necessarily fully and **properly** implemented!!!
- Documentation could be improved and extended.
- If you feel like contributing, that is great, please send me an email: sebastien.cunnac@ird.fr
- tantale has been written with only Linux systems in mind and will very likely not work on other OS (eg Windows)
- Some of tantale wrappers use code written in other languages. For them to work, you must ideally have **Java and Perl on the PATH** in your system.
- **Conda and Mamba** must be installed as well.
- Perl scripts (DisTAL et al.) used in tantale require a number of perl libraries. Upon experimenting a bit we realized it was more convenient to include those libraries in tantale. The consequence is that those files cause tantale to occupy a quite some disk space...
    
    

