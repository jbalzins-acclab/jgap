#!/bin/bash
set -e

INSTALL_DIR="$(pwd)/tabgap_install"

# Check if tabgap is already installed
if [ -f "$INSTALL_DIR/tabgap/tabgap/tabulate.py" ]; then
    echo "TabGAP is already installed at $INSTALL_DIR/tabgap."
    echo "Skipping installation. If you want to reinstall, remove the $INSTALL_DIR directory first."
    exit 0
fi

# Clean up any partial installations
if [ -d "$INSTALL_DIR" ]; then
    echo "Removing incomplete/old installation at $INSTALL_DIR..."
    rm -rf "$INSTALL_DIR"
fi

echo "Setting up installation directory..."
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "Cloning TabGAP..."
git clone https://gitlab.com/jezper/tabgap.git

cd tabgap
export TABGAP_ROOT="$(pwd)"

echo "Checking out TabGAP to specific commit cbdf7532f3ca29c0aa5a5ae668122ad0a1b7fa1f..."
git fetch origin
git checkout cbdf7532f3ca29c0aa5a5ae668122ad0a1b7fa1f

cd ..

echo "=================================================="
echo "TabGAP cloned successfully to $TABGAP_ROOT"
echo "=================================================="
