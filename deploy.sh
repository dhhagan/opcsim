#!/usr/bin/env bash

git checkout master

echo "Building and deploying docs"
cd docs
make clean
make tutorials
make html
cd ..

# commit and push
git add -A
git commit -m "building and pushing docs"
git push origin master

# switch branches
git checkout gh-pages
rm -rf *
touch .nojekyll
git checkout master docs/_build/html
mv ./docs/_build/html/* ./
rm -rf ./docs

git add -A
git commit -m "publish new docs"
git push origin gh-pages

# switch back to docs
git checkout master
