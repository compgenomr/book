# works on master branch
# use it after bookdown::render() or knitting index.Rmd
git clone -b gh-pages   https://github.com/compgenomr/book.git book-output
cd book-output
cp -r ../_book/* ./
git add *
git commit -m "Update the book manually 2"
git push origin gh-pages
