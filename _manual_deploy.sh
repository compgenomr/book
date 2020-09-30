# works on master branch
# use it after bookdown::render() or knitting index.Rmd
rm -rf book-output
git clone -b gh-pages   https://github.com/compgenomr/book.git book-output
cd book-output
cp -r ../_book/* ./
git add *
git commit -m "Update the book manually 3"
git push origin gh-pages

