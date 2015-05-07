docdir = $(pkgdatadir)/doc

doc_DATA = 4ti2_manual.pdf

## Manual sources.
HELPFILES = output-help.txt gensymm-help.txt genmodel-help.txt ppi-help.txt \
	hilbert-help.txt graver-help.txt zsolve-help.txt

EXTRA_DIST =					\
	4ti2_manual.pdf				\
	3by3magic.eps				\
	4ti2_manual.bib				\
	4ti2_manual.tex				\
	4ti2_manual_api.tex			\
	4ti2_manual_advanced.tex		\
	4ti2_manual_beginner.tex		\
	4ti2_manual_ref.tex			\
	$(HELPFILES)				\
	cubepicture.eps				\
	feasible.pdf				\
	feasible.pdf_t				\
	feasible_LP.pdf				\
	feasible_LP.pdf_t			\
	magicsqs.eps				\
	magicsqs.fig				\
	magicsqs.pdf				\
	magicsqs.pdf_t				\
	magicsqs_hil.fig			\
	magicsqs_hil.pdf			\
	magicsqs_hil.pdf_t			\
	magicsqs_ray.fig			\
	magicsqs_ray.pdf			\
	magicsqs_ray.pdf_t			\
	myfancyheadings.sty			\
	solved_IP.pdf				\
	solved_IP.pdf_t				\
	solved_LP.pdf				\
	solved_LP.pdf_t

EXTRA_DIST +=					\
	deloera-hemmecke-koeppe-book.tex	\
	deloera-hemmecke-koeppe-cover.jpg

### We do not attempt to build the manual in Make targets.
### Before making a release, pdflatex the manual; we distribute the PDF.

.PHONY: update-help update-manual 

update-manual: update-help
	-(cd $(srcdir) && pdflatex 4ti2_manual; pdflatex 4ti2_manual)

update-help:
	../src/util/gensymm --help > $(srcdir)/gensymm-help.txt
	../src/util/genmodel --help > $(srcdir)/genmodel-help.txt
	../src/util/output --help > $(srcdir)/output-help.txt
	../src/ppi/ppi --help > $(srcdir)/ppi-help.txt
	../src/zsolve/zsolve --help > $(srcdir)/zsolve-help.txt
	../src/zsolve/hilbert --help > $(srcdir)/hilbert-help.txt
	../src/zsolve/graver --help > $(srcdir)/graver-help.txt
	../src/groebner/qsolve --help > $(srcdir)/qsolve-help.txt
	../src/groebner/rays --help > $(srcdir)/rays-help.txt
	../src/groebner/circuits --help > $(srcdir)/circuits-help.txt
	../src/groebner/groebner --help > $(srcdir)/groebner-help.txt
	../src/groebner/markov --help > $(srcdir)/markov-help.txt
	../src/groebner/minimize --help > $(srcdir)/minimize-help.txt
	../src/groebner/normalform --help > $(srcdir)/normalform-help.txt
	../src/groebner/walk --help > $(srcdir)/walk-help.txt
	../src/groebner/zbasis --help > $(srcdir)/zbasis-help.txt
