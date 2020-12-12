.rs.unloadPackage("raster")
.rs.unloadPackage("sp")


writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")

# uninstall sp

install.packages("sp", type = "source")

# delete and then reinstall, sp, rlang, Rcpp, backports
install.packages('rlang')
install.packages('Rcpp')
install.packages('glue')
install.packages('backports')
install.packages('processx')
install.packages('ps')
install.packages("digest")


devtools::install_github('biomodhub/biomod2')
