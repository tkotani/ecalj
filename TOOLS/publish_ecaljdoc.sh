#!/bin/bash
# Publish the documents ecalj/ecaljdoc/ to github.com/ecalj/ecaljdoc (its main is the GitHub Pages site
# https://ecalj.github.io/ecaljdoc/). 2026-10-02: replaces git subtree (procedure: MD/ecaljdoc_publish.md).
#
#   TOOLS/publish_ecaljdoc.sh          copy the committed ecaljdoc/ of HEAD into the publishing clone and commit there; show the change
#   TOOLS/publish_ecaljdoc.sh --push   ... and push it (only when the maintainer says so)
#
# Env: ECALJDOC_PUBLISH_DIR (default ~/work/ecaljdoc_publish, a clone of the publishing repository, made if missing),
#      ECALJDOC_REMOTE (default git@github.com:ecalj/ecaljdoc.git).
# What is published is what HEAD of ecalj has under ecaljdoc/ (git archive), not the working tree: commit first.
# The publishing clone is reset to its origin/main before the copy, so edits made directly on GitHub are overwritten
# (make them in ecalj/ecaljdoc instead).
set -eu
TOP=$(git -C "$(dirname "$(readlink -f "$0")")" rev-parse --show-toplevel)
PUB=${ECALJDOC_PUBLISH_DIR:-$HOME/work/ecaljdoc_publish}
REMOTE=${ECALJDOC_REMOTE:-git@github.com:ecalj/ecaljdoc.git}
PUSH=no; [ "${1:-}" = "--push" ] && PUSH=yes

if [ -n "$(git -C "$TOP" status --porcelain -- ecaljdoc)" ]; then
  echo "publish_ecaljdoc: ecaljdoc/ has uncommitted changes; they would not be published. Commit them first:"
  git -C "$TOP" status --short -- ecaljdoc | head
  exit 1
fi
REV=$(git -C "$TOP" rev-parse --short=9 HEAD)

if [ ! -d "$PUB/.git" ]; then
  git clone -q "$REMOTE" "$PUB"
fi
git -C "$PUB" fetch -q origin
git -C "$PUB" checkout -q main
git -C "$PUB" reset -q --hard origin/main

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
git -C "$TOP" archive HEAD ecaljdoc | tar -x -C "$TMP"
rsync -a --delete --exclude .git "$TMP/ecaljdoc/" "$PUB/"

git -C "$PUB" add -A
if git -C "$PUB" diff --cached --quiet; then
  echo "publish_ecaljdoc: nothing to publish (ecaljdoc/ of ecalj $REV is the same as $REMOTE main)"
  exit 0
fi
git -C "$PUB" commit -q -m "Sync from ecalj $REV ($(date '+%F %H:%M'))"
echo "publish_ecaljdoc: committed in $PUB:"
git -C "$PUB" show --stat --format='%h %s' HEAD | tail -n +1 | head -40
if [ "$PUSH" = yes ]; then
  git -C "$PUB" push origin main
  echo "publish_ecaljdoc: pushed to $REMOTE (the site is rebuilt by its GitHub Actions)"
else
  echo "publish_ecaljdoc: not pushed. To publish: TOOLS/publish_ecaljdoc.sh --push   (or git -C $PUB push origin main)"
fi
