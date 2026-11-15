#!/bin/sh
# Fetch OXSA (Purvis et al. 2017; https://github.com/OXSAtoolbox/OXSA) into ./OXSA,
# alongside this repository's two prior-knowledge files (upstream has no files of those names).
cd "$(dirname "$0")"
tmp=$(mktemp -d)
git clone --depth 1 -q https://github.com/OXSAtoolbox/OXSA "$tmp"
git -C "$tmp" checkout -q 7cf1fdd30d9f4775f52c41af70cf9960d2db6cd7   # version tested
rm -rf "$tmp/.git"
cp -R "$tmp"/. OXSA/
rm -rf "$tmp"
