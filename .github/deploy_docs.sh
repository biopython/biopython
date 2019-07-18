#!/bin/bash

if [[ $(git rev-parse --abbrev-ref HEAD) != "master" ]]; then
    echo "Not on master branch, not attempting to deploy documentation"
    exit 0
fi

# Assumes being called from the Biopython repository's root folder,
# (i.e. a clone of https://github.com/biopython/biopython) as part
# of our continuous integration testing to save the compiled docs
# to https://github.com/biopython/docs
#
# i.e. AFTER you have run sphinx-build with:
#
# $ make -C Doc/api/ html
#
# First, we need to enable SSH. First, generate a key locally (You must use an empty passphrase):
#
#     ssh-keygen -t rsa -b 4096 -C "biopython documentation deploy key" -f biopython_doc_key -N ""
#
# Output from the previous command should look something like this:
#
#     Generating public/private rsa key pair.
#     Your identification has been saved in biopython_doc_key
#     Your public key has been saved in biopython_doc_key.pub
#     The key fingerprint is:
#     SHA256:abunchofnumbersandletters biopython documentation deploy key
#     The key's randomart image is:
#     +---[RSA 4096]----+
#     |    xxx x        |
#     |    xxxx x       |
#     | . . .x.x x x    |
#     |  x . xO O o     |
#     |   x  x x o .    |
#     |  .  . + . o .   |
#     | x      E o . .  |
#     |xx.. o xxxx  x   |
#     |x.  . ..+xxx     |
#     +----[SHA256]-----+
#
# Now we're going to add the public key to the docs project. To get the key, run the
# following command:
#
#     cat biopython_doc_key.pub
#
# and copy the output. Go to https://github.com/$USER/docs and
# substitute $USER for the organization your doing this for (biopython for
# the main repository, your GitHub username for your own fork).
# Click on Settings -> Deploy Keys -> Add deploy key
# In the title, pick a name that's descriptive and easy for you to remember.
# I'm going to use "biopython_doc_key". Paste the contents of the public
# key (the one you copied from the cat command a few moments ago) into the
# "Key" text box. Check "Allow write access" and click "Add key".
#
# Now, we have to add the private key to CircleCI. Go to
# https://app.circleci.com/settings/project/github/$USER/biopython
# once again replace $USER with either your github user name or
# biopython. Click "SSH Keys" in the menu on the left and click
# the button "Add SSH Key" at the bottom of the page. Leave the
# host name blank. In your terminal, run the following command:
#     cat biopython_doc_key
# and copy the contents. Back in the CircleCI UI, paste the contents
# you just copied into the "Private Key" text box. Click "Add SSH Key".
#
# We now have to set an environment variable called DOC_KEY
# in the Biopython project settings. To get this value in a format
# we can paste, run the following command and copy the output:
#
#     python -c "print(open('biopython_doc_key').read().strip().replace(' ', r'\ ').replace('\n', r'\\\n'))"
#
# In the CircleCI UI, click on "Project Settings" in the top right, 
# and go to "Environment Variables" on the left. Click "Add Environment
# Variable". In the "Name" input, put "DOC_KEY" and paste the value 
# you copied earlier into the "Value" input.
# 
# Finally, we need to add a User Key to the CircleCI biopython project.
# To do this, go to Biopython project in the CircleCI UI, click
# "project settings" in the top right, click "SSH Keys", and click
# "Add User Key." This allows CircleCI to act as a Git user and
# push the project to the docs repository.

set -e

if [ -z "$DOC_KEY" ]; then
    echo "Missing (secret) environment variable DOC_KEY,"
    echo "which should hold the private SSH deployment key."
    false
fi

set -euo pipefail

DEST_SLUG=biopython/docs
DEST_DIR=`python -c "import Bio; v=Bio.__version__; print('dev' if 'dev' in v else v)"`
SOURCE_DIR=${BUILD_DIR:-$PWD}/Doc/_build/html
WORKING_DIR=/tmp/deploy_biopython_docs
COMMIT_HASH=$(git rev-parse HEAD) # For later when we commit the docs repository

if [ -z "$DEST_DIR" ]; then
   echo "ERROR: Failed to get Biopython version, is it not installed?"
   python -c "import Bio; print(Bio.__version__)"
   false
fi

DEST_DIR=$DEST_DIR/
echo "Aiming to deploy $SOURCE_DIR to $DEST_SLUG branch gh-pages as $DEST_DIR"

# We have to create the SSH key with spaces and new lines, so
# we un-escape whitespace and newline characters to recover
# the SSH deploy key:
python -c "import os; print(os.environ['DOC_KEY'].strip().replace(r'\ ', ' ').replace(r'\n', '\n'))" > $HOME/.biopython_doc_deploy.key
# Check we have a sane looking line structure:
if [ `grep -c "^\-\-\-\-\-" $HOME/.biopython_doc_deploy.key` -ne 2 ]; then
    echo "ERROR: Failed to rebuild the SSH key,"
    wc -l $HOME/.biopython_doc_deploy.key
    md5sum $HOME/.biopython_doc_deploy.key
    false
fi

chmod 600 $HOME/.biopython_doc_deploy.key
mv $HOME/.biopython_doc_deploy.key $HOME/.ssh/id_rsa

# Clone the destination under /tmp (public URL, no key needed)
rm -rf $WORKING_DIR
git clone https://github.com/$DEST_SLUG.git $WORKING_DIR
pushd $WORKING_DIR
git checkout gh-pages
# Switch the git protocol to SSH based so we can use our key
git remote set-url origin --push git@github.com:$DEST_SLUG.git
popd

echo "Copying $SOURCE_DIR/* to $WORKING_DIR/$DEST_DIR/ next"
if [ ! -d $SOURCE_DIR ]; then
    echo "ERROR: Directory $SOURCE_DIR/ does not exist."
    false
fi

# Remove any old files
pushd $WORKING_DIR
if [ -d $DEST_DIR ]; then
    echo "Removing old files"
    git rm -r $DEST_DIR/
fi
mkdir -p $DEST_DIR
echo "Copying files"
cp -R $SOURCE_DIR/* $DEST_DIR/
echo "Staging files in git"
git add $DEST_DIR/

echo "Preparing to deploy documentation"

if [[ -z $(git status --porcelain) ]]; then
    echo "Nothing has changed, nothing needs pushing."
else
    echo "Making commit of new files"
    git config user.email "sphinx@example.org"
    git config user.name "Sphinx"
    git commit -m "Automated update ${COMMIT_HASH}"
    echo "Finally, pushing to $DEST_SLUG gh-pages branch"
    git push origin gh-pages
    echo "Documentation deployed!"
fi

popd
