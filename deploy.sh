#!/usr/bin/env bash

git checkout docs

echo "Building and deploying docs"
cd docs
make clean
make notebooks
make html
cd ..

# commit and push
git add -A
git commit -m "building and pushing docs"
git push origin docs

# switch branches
git checkout gh-pages
rm -rf *
touch .nojekyll
git checkout docs/_build/html
mv ./docs/_build/html/* ./
rm -rf ./docs

git add -A
git commit -m "publish new docs"
git push gh-pages

# switch back to docs
git checkout docs
