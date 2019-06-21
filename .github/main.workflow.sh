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


# setup SSH
mkdir ~/.ssh
chmod 700 ~/.ssh
echo "$MIRROR_DEPLOYMENT_KEY" > ~/.ssh/id_ed25519
chmod 600 ~/.ssh/id_ed25519
git config --global core.sshCommand "ssh -i ~/.ssh/id_ed25519 -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"

# Ensure our local repo has EV-ER-Y-THING
git fetch --tags --prune --prune-tags --force --update-head-ok --progress

# Setup the mirror remote
git remote set-url --push origin "$MIRROR_URL"

# # Push to the mirrored repository
# git push --mirror --force --progress
