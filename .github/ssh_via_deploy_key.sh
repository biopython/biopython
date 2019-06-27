#!/bin/bash
# Call ssh using our GitHub repository deploy key (set via -i)
# using -F to make sure this ignores ~/.ssh/config
ssh -i "$HOME/.biopython_doc_deploy.key" -F /dev/null -p 22 $*
