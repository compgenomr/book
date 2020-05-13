#!/bin/sh

set -ev

Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book',output_dir='book_pdf')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book',new_session=TRUE)"