## tutorial.html takes about 3 minutes
## meta-analysis.html takes about 20 minutes

all: tutorial.html meta-analysis.html 

%.docx: %.md
	Rscript --vanilla build-docx.r $< $@

%.html %.md: %.rmd
	Rscript --vanilla build-html.r $< $@

clean: 
	rm -rf *.md *.html *.docx *.pdf *.Rout	


