#!/bin/sh

set -ev

# render gitbook
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook',new_session=TRUE)"

# render pdf
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book',output_dir='book_pdf')"

# compile each chapter separetely and merge
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book',new_session=TRUE)"

# compile one chapter
#Rscript -e "bookdown::preview_chapter('01-intro2Genomics.Rmd', 'bookdown::pdf_book',new_session=TRUE)"
