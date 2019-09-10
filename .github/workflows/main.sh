#!/bin/bash

set -o errexit

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

echo "Branches found:"
git branch -avvv

echo "Setting up SSH"
[ -d ~/.ssh ] || mkdir -p ~/.ssh
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

echo "Verifying ssh access to github"
ssh -i ~/.ssh/id_ed25519 -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -T git@github.com || true

echo "Creating a local mirror of ${GITHUB_REPOSITORY}"
cd ~ || exit 1
git clone --mirror "git@github.com:${GITHUB_REPOSITORY}.git"
pwd
ls -ld "${GITHUB_REPOSITORY#*/}.git"

git clone --verbose --mirror "$MIRROR_URL"


cd "${GITHUB_REPOSITORY#*/}.git" || exit 1

echo "Setting mirror remote url"
git remote set-url origin "${MIRROR_URL}"
git remote -v

echo "Pruning PR refs"
git show-ref | cut -d' ' -f2 | grep 'refs/pull/' | xargs -r -L1 git update-ref -d
git show-ref

git config --show-origin --list

echo "Attempting push to MIRROR repository..."
if ! git push --mirror --force --progress ; then
    sleep 25
    git push --mirror --force --progress || exit 78 # nuetral exit
fi
