#!/bin/bash

# Assumes being called from the Biopython repository's root folder,
# (i.e. a clone of https://github.com/biopython/biopython) as part
# of our continuous integration testing to save the compiled docs
# to https://github.com/biopython/docs
#
# In order to have write permissions, we put a private key into the
# TravisCI settings as a secure environment variable, and put the
# matching public key into the GitHub documentation repository's
# settings as a deploy key with write permissions.
#
# Key creation,
#
# $ ssh-keygen -t rsa -b 4096 -C "biopython documentation deploy key" -f biopython_doc_key -N ""
# Generating public/private rsa key pair.
# Your identification has been saved in biopython_doc_key.
# Your public key has been saved in biopython_doc_key.pub.
# The key fingerprint is:
# SHA256:nFfhbwryDLDz8eDEHa4sjdH0gOgwyXGGDUBGfDi5luQ biopython documentation deploy key
# The key's randomart image is:
# +---[RSA 4096]----+
# |===+o       .    |
# |.B.*.. .   . .   |
# |o X . o o . o    |
# | E +   B * o .   |
# |.   . + S *   o  |
# |       X @ . o   |
# |      o * + .    |
# |       .         |
# |                 |
# +----[SHA256]-----+
#
# Next, we add the public key to https://github.com/biopython/docs as
# a deployment key with write permission,
#
# $ cat biopython_doc_key.pub
# ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAACAQDpQ3I6ZpL9cqUpqkHMPALTQg9ya3sL1MVXYjbTnuWnDoRml5UYXVgD8hgOJxwaxDo1BV+fKn68LXPEwlZ5FC6eRSCJz20SvWPkMDhAwChJJ+nE7f/vvK18R3Ge9ksWra8LFSR3EL7joQTN+c1VyaJH22qj1OED3G68Ix+bgvnUpZgeurV8vDV06FVx7H1Q5a5MoTWFdMa9wzJn5o6m7khditOTDKqznFULoOONpw7CsTiJD6drQPk1pwftDxEBMEAG7cKwux/dRWJtzsRQ7IO0d/AhzsqnLJJIgkHzQwmvpGpffWfoomNwF4bWJuWzu6tRcGcX16fLMyGK8kFJaL1zY6gQFkAbfsIdA2G28S79mIC4jT1JtiNYBOV9wIjxyZUyvzSeQGVC7yBafHE5eEb267dgGnDl654XzyIImLSKv/nx8No16UK/e5F+ds3hp0DPTknzeVOGBUEt1k8pEp47J9JVKoeceph0cJbfzFNv9pgOgyaHb1mhs9pI4kIQ3R+ibeAZbPWT709n26Y99Q2MSSZyPuZvX8VBA1NfoENmuTrEn/qqGlvZez3Blh4MIvYg24DFv/rHN92Edk5S7xY0eB7E6D6X/N80ThuBSqxlJpxSQlA+LICcq/EPd37/WT7exiheXysN5oIOvwNgUNNFftDWv2gPBu2bf/foHfAQKQ== biopython documentation deploy key
#
# Finally, we add the private key to TravisCI by going to
# https://travis-ci.org/biopython/biopython/settings or any authorised
# fork like https://travis-ci.org/peterjc/biopython/settings and
# setting DOC_KEY to the following (secret) value:
#
# $ python -c "print(open('biopython_doc_key').read().strip().replace(' ', r'\ ').replace('\n', r'\\\n'))"
# ...
#
# TravisCI requires we escape spaces as '\ ' and newlines as '\\n', and
# we explicitly strip the trailing new line so that we don't get an extra
# one when rebuilding the key later.
#
# Make sure "DISPLAY VALUE IN BUILD LOG" is off (the default).
#
# For testing locally, set local environment $DOC_KEY to this value.
# Thereafter, when ever this script gets run on TravisCI it should
# be able to deplop the HTML documentation to our documentation
# repository (which will dispaly on biopython.org via GitHub pages).

set -e

if [ -z "$DOC_KEY" ]; then
    echo "Missing (secret) environment variable DOC_KEY,"
    echo "which should hold the private SSH deployment key."
    false
fi

set -euo pipefail

DEST_SLUG=biopython/docs
# Could look at $TRAVIS_TAG, e.g. DEST_DIR=${TRAVIS_TAG:-dev}
# However, tags lack the dots in the version number. Since
# Biopython was installed to run Sphinx and build the docs,
# can use this:
DEST_DIR=`python -c "import Bio; v=Bio.__version__; print('dev' if 'dev' in v else v)"`
SOURCE_DIR=${TRAVIS_BUILD_DIR:-$PWD}/Doc/api/_build/html
WORKING_DIR=/tmp/deploy_biopython_docs

if [ -z "$DEST_DIR" ]; then
   echo "ERROR: Failed to get Biopython version, is it not installed?"
   python -c "import Bio; print(Bio.__version__)"
   false
fi
DEST_DIR=$DEST_DIR/api
echo "Aiming to deploy $SOURCE_DIR to $DEST_SLUG branch gh-pages as $DEST_DIR"

# On TravisCI, must create the variable using '\ ' and '\n', so
# here we must unescape the whitespace to recover the SSH deploy key:
python -c "import os; print(os.environ['DOC_KEY'].strip().replace(r'\ ', ' ').replace(r'\n', '\n'))" > $HOME/.biopython_doc_deploy.key
# Check we have a sane looking line structure:
if [ `grep -c "^\-\-\-\-\-" $HOME/.biopython_doc_deploy.key` -ne 2 ]; then
    echo "ERROR: Failed to rebuild the SSH key,"
    wc -l $HOME/.biopython_doc_deploy.key
    md5sum $HOME/.biopython_doc_deploy.key
    false
fi
chmod 600 $HOME/.biopython_doc_deploy.key
export GIT_SSH=${TRAVIS_BUILD_DIR:-$PWD}/.github/ssh_via_deploy_key.sh

if ! [[ -f "$GIT_SSH" ]]; then
    echo "Error, set GIT_SSH="$GIT_SSH" but does not exist"
    false
elif ! [[ -x "$GIT_SSH" ]]; then
    echo "Error, set GIT_SSH="$GIT_SSH" but not executable"
    false
fi;

echo "Setting up clone of $DEST_SLUG locally at $WORKING_DIR"

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

if [[ -z $(git status --porcelain) ]]; then
    echo "Nothing has changed, nothing needs pushing."
else
    echo "Making commit of new files"
    git commit -m "Automated update ${TRAVIS_COMMIT:-}" --author "TravisCI <travisci@example.org>"
    echo "Finally, pushing to $DEST_SLUG gh-pages branch"
    git push origin gh-pages
    echo "Documentation deployed!"
fi

popd
