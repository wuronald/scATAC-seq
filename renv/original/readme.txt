The renv.lock file is from the ArchR website specific for the R version 4.4. 
See below:
download.file(url = "https://pub-9ae435458ecc412abbbc9420a502ec38.r2.dev/renv.lock", destfile = "./renv.lock")
It is missing installation for ArchR version 1.0.3, so that needs to be done manually after renv mediated install of the other dependencies.