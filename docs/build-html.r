library(knitr)
library(markdown)

knit.report <- function(input.filename, output.filename,
                        stylesheet=getOption("markdown.HTML.stylesheet"), ...) {
    input.filename <- normalizePath(input.filename)
    output.filename <- normalizePath(output.filename)

    output.dir <- dirname(output.filename)
    if (!file.exists(output.dir))
        dir.create(output.dir)

    current.dir <- getwd()
    on.exit(setwd(current.dir))
    setwd(output.dir)

    name <- gsub("\\.[^.]+$", "", basename(output.filename))
    suffix <- gsub(".*\\.([^.]+)$", "\\1", output.filename)
    is.html <- tolower(suffix) %in% c("htm","html")

    if (is.html)
        md.filename <- paste(name, "md", sep=".")
    else
        md.filename <- basename(output.filename)
    
    knit(input.filename, output=md.filename, envir=parent.frame(), ...)

    if (is.html)
        markdownToHTML(md.filename, basename(output.filename), stylesheet=stylesheet)
}

args = commandArgs(trailingOnly=TRUE)
knit.report(args[1], args[2])

