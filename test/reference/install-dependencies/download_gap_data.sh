#!/bin/bash

INSTALL_DIR="$(pwd)/gap_data"

# Check if gap-data is already downloaded
if [ -d "$INSTALL_DIR" ] && [ -f "$INSTALL_DIR/README.md" ]; then
    echo "gap-data is already cloned at $INSTALL_DIR."
    echo "Skipping cloning. If you want to re-clone, remove the $INSTALL_DIR directory first."
    exit 0
fi

echo "Setting up gap-data directory..."
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "Cloning gap-data..."
git init
git remote add origin https://gitlab.com/acclab/gap-data.git
echo "Fetching specific commit 3072bbea41b636689580ce939ae76785904102ff..."
git fetch origin 3072bbea41b636689580ce939ae76785904102ff
git reset --hard FETCH_HEAD

cd ..

echo "=================================================="
echo "gap-data cloned successfully to $INSTALL_DIR"
echo "=================================================="
