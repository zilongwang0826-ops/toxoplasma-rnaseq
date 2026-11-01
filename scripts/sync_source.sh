#!/usr/bin/env bash
# Transfer committed public source only. Requires explicit remote-write approval.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${REMOTE_HOST:?Set an SSH host or user@host}"
: "${REMOTE_WORKSPACE:?Set the absolute dedicated remote project path}"
# SSH identity/options should be supplied through RSYNC_RSH and SSH_ARGS, never key contents.
[[ "$REMOTE_WORKSPACE" =~ ^/home/[A-Za-z0-9_-]+/Projet2-Toxoplasma$ ]] || { echo 'Unexpected remote workspace' >&2; exit 2; }
[[ "$REMOTE_HOST" =~ ^[A-Za-z0-9_.@-]+$ ]] || exit 2
if [[ -n $(git status --porcelain --untracked-files=no) ]]; then echo 'Commit source changes first' >&2; exit 2; fi
revision=$(git rev-parse HEAD)
mkdir -p .cache
bundle=".cache/source-$revision.bundle"
git bundle create "$bundle" HEAD
read -r -a ssh_options <<< "${SSH_ARGS:--o BatchMode=yes -o StrictHostKeyChecking=yes}"
ssh "${ssh_options[@]}" "$REMOTE_HOST" "mkdir -p '$REMOTE_WORKSPACE/.transfer'"
rsync -av --checksum "$bundle" "$REMOTE_HOST:$REMOTE_WORKSPACE/.transfer/"
ssh "${ssh_options[@]}" "$REMOTE_HOST" bash -s -- "$REMOTE_WORKSPACE" "$revision" <<'REMOTE'
set -euo pipefail
cd "$1"
if [[ -d .git ]]; then
    [[ -z $(git status --porcelain --untracked-files=no) ]] || { echo 'Remote source changes detected; reconcile first' >&2; exit 2; }
    git fetch ".transfer/source-$2.bundle" HEAD
    git merge --ff-only FETCH_HEAD
else
    git init -b main
    git fetch ".transfer/source-$2.bundle" HEAD
    git checkout -B main FETCH_HEAD
fi
[[ $(git rev-parse HEAD) == "$2" ]]
printf 'Synchronized commit: %s\n' "$2"
REMOTE
