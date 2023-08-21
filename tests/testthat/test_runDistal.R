

talOrfsFile <- system.file("extdata", "SamplePutativeTalOrf.fasta", package = "tantale", mustWork = T)

# run DisTal
distal_output <- runDistal(fasta.file = talOrfsFile, condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")

disTal <- file.path(system.file("tools", "DisTAL1.2_MultipleAlignment", package = "tantale", mustWork = T),
                    "DisTAL_v1.2_matest_M.pl")
lib <- "/home/cunnac/Lab-Related/MyScripts/tantale/inst/tools/DisTAL1.2_MultipleAlignment/lib/Algorithm"
disTalCMD <- glue::glue("perl {disTal} -I {lib} -h")

res <- systemInCondaEnv(envName = "perlforal",
                        condaBinPath = condaBinPath,
                        command = disTalCMD,
                        ignore.stdout = T)
