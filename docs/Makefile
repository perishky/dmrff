all: tutorial.html tutorial.docx

%.docx: %.md
	Rscript --vanilla build-docx.r $< $@

%.html %.md: %.rmd
	Rscript --vanilla build-html.r $< $@

clean: 
	rm -rf *.md *.html *.docx *.pdf *.Rout	


