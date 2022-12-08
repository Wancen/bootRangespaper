#!/bin/bash
for doc in supp main plain
do
  pdflatex $doc
  bibtex $doc
  pdflatex $doc
  pdflatex $doc
done
