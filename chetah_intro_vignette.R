# http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/
# https://rmarkdown.rstudio.com/articles_report_from_r_script.html
# https://sachsmc.github.io/knit-git-markr-guide/knitr/knit.html
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)
knitr::opts_knit$set(
  root.dir = system.file('data', package = 'CHETAH')
)
library(Matrix)
library(CHETAH)
