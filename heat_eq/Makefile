all: report

report: report.md 
	pandoc report.md -o report.pdf --from markdown --template eisvogel.tex --highlight-style zenburn --variable urlcolor=cyan

clean:
	rm -f report.pdf

remake : clean all

