#!/bin/bash

set -o errexit

# configure git
git config --global user.name "Sourcery-Bot"
git config --global user.email "si-bot@izaakbeekman.com"

# Print diagnostic info
echo "Workflow name: $GITHUB_WORKFLOW"
echo "Action name: $GITHUB_ACTION"
echo "Person/app initiating action: $GITHUB_ACTOR"
echo "Current repository: $GITHUB_REPOSITORY"
echo "Event name triggering workflow: $GITHUB_EVENT_NAME"
echo "Workspace path: $GITHUB_WORKSPACE"
echo "GITHUB_SHA: $GITHUB_SHA"
echo "GitHub branch/tag/ref: $GITHUB_REF"
echo "Current directory: $(pwd)"

echo "Contents of Workspace:"
ls -al "$GITHUB_WORKSPACE"

echo "Git version:"
git --version

echo "Git Status:"
git status

echo "Checking git remotes"
git remote -v
git remote set-url origin https://${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY}.git

echo "Fetch everything and make sure we're up-to-date before mirroring."
git fetch --tags --prune --prune-tags --force --update-head-ok --progress

echo "Branches found:"
git branch -avvv

echo "Seting up the mirror remote..."
git remote set-url --push origin "$MIRROR_URL"

echo "Setting up SSH"
mkdir ~/.ssh
chmod 700 ~/.ssh
echo "$IBB_PWLESS_DEPLOY_KEY" > ~/.ssh/id_ed25519
chmod 600 ~/.ssh/id_ed25519
echo "Checking the sha256 checksum of the ssh key..."
sha256sum ~/.ssh/id_ed25519
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

echo "Configure git for authorized user"
git config --global user.name "Izaak Beekman"
git config --global user.email "ibeekman@paratools.com"
git config --global core.sshCommand "ssh -i ~/.ssh/id_ed25519 -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"
git config --show-origin --list
git remote -v

echo "Testing ssh connection to mirror repo"
ssh -i ~/.ssh/id_ed25519 -T git@github.com || echo "FAILED TO AUTHENTICATE!!!"

echo "Attempting push to MIRROR repository..."
# Push to the mirrored repository
git push --mirror --force --progress
