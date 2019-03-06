#!/bin/sh

set -e

[ -z "${GITHUB_PAT}" ] && exit 0
[ "${TRAVIS_BRANCH}" != "master" ] && exit 0

BOOK_DIR=$(pwd)/_book
rm -rf ~/_book
mkdir ~/_book && cd ~/_book
git clone -b gh-pages https://${GITHUB_PAT}@github.com/${TRAVIS_REPO_SLUG}.git .
ls | grep -v ^bookdown[.].* | xargs rm -rf
git ls-files --deleted -z | xargs -0 git rm
cp -r ${BOOK_DIR}/* ./
git add --all *
git commit -m"update homepage (travis build ${TRAVIS_BUILD_NUMBER})"
git push -q -f origin gh-pages
