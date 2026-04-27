#!/bin/bash
set -e

echo "Starting benchmarking setup..."

# Ensure all install scripts are executable
chmod +x download_gap_data.sh
chmod +x install_quip.sh
chmod +x install_tabgap.sh
chmod +x install_lammps.sh

echo "========================================"
echo "1. Downloading gap-data"
echo "========================================"
./download_gap_data.sh

echo "========================================"
echo "2. Installing QUIP & GAP"
echo "========================================"
./install_quip.sh

echo "========================================"
echo "3. Installing TabGAP"
echo "========================================"
./install_tabgap.sh

echo "========================================"
echo "4. Installing LAMMPS"
echo "========================================"
./install_lammps.sh

echo "========================================"
echo "All installations and downloads completed successfully!"
echo "========================================"

echo "========================================"
echo "Copying jgap/build to jgap_install"
echo "========================================"
mkdir -p jgap_install
cp -r ../build jgap_install/

# Add benchmarking test execution logic here in the future
