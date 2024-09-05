
%.pdf: %/slides.md slides.sty
	pandoc -t beamer -V theme:Boadilla --pdf-engine=lualatex -V fontsize=10pt -V colorlinks:true --listings -V 'mainfont:NewComputerModernSans10' -fmarkdown-implicit_figures -H $(word 2,$^) -s $< -o $@

.PHONY: clean

clean:
	rm -f *.pdf

.PHONY: clean
