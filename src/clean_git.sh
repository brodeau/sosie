#!/bin/bash

git filter-branch --index-filter 'git rm --cached --ignore-unmatch *.gz' -- --all
rm -Rf .git/refs/original
rm -Rf .git/logs/
git gc --aggressive --prune=now
git count-objects -v
