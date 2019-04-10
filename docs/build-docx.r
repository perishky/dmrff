library(rmarkdown)
library(knitr)

args = commandArgs(trailingOnly=TRUE)

render(args[1], output_file=args[2], output_format="word_document")

