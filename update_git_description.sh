#!/bin/bash

git describe --always --tags --dirty > git_description.new
if [[ -f "git_description" ]]; then
   cmp --silent git_description.new git_description && rm git_description.new || mv git_description.new git_description
else
   mv git_description.new git_description
fi
