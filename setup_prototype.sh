#!/bin/bash
# setup_prototype.sh
# Automated setup script for DCEA prototype (Linux/Mac)

set -e  # Exit on error

echo "======================================"
echo "TCGA-BRCA DCEA Prototype Setup"
echo "======================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Step 1: Create directory structure
echo -e "${YELLOW}[1/5] Creating directory structure...${NC}"
mkdir -p scripts
mkdir -p input/preprocessed/matrices
mkdir -p output/{correlation_matrices,differential_results,hub_subgraphs,visualizations}
touch output/prototype_execution.log
echo -e "${GREEN}✓ Directories created${NC}"
echo ""

# Step 2: Create virtual environment
echo -e "${YELLOW}[2/5] Setting up Python virtual environment...${NC}"
if [ ! -d "venv_proto" ]; then
    python3 -m venv venv_proto
    echo -e "${GREEN}✓ Virtual environment created${NC}"
else
    echo -e "${GREEN}✓ Virtual environment already exists${NC}"
fi
echo ""

# Step 3: Install dependencies
echo -e "${YELLOW}[3/5] Installing Python dependencies...${NC}"
source venv_proto/bin/activate
pip install --quiet --upgrade pip
pip install --quiet -r requirements_proto.txt
echo -e "${GREEN}✓ Dependencies installed${NC}"
echo ""

# Step 4: Verify dependencies
echo -e "${YELLOW}[4/5] Verifying installation...${NC}"
python -c "import numpy, pandas, networkx, matplotlib, scipy, yaml; print('✓ All core modules imported successfully')"
echo ""

# Step 5: Check for input data
echo -e "${YELLOW}[5/5] Checking for input data...${NC}"
if [ -f "input/preprocessed/matrices/tumor_matrix.tsv" ] && [ -f "input/preprocessed/matrices/normal_matrix.tsv" ]; then
    echo -e "${GREEN}✓ Input matrices found${NC}"
    echo ""
    echo -e "${GREEN}======================================"
    echo "Setup complete!"
    echo "======================================"
    echo ""
    echo "Next steps:"
    echo "  1. Activate environment: source venv_proto/bin/activate"
    echo "  2. Run prototype: python scripts/proto_01_prepare_data.py"
    echo ""
else
    echo -e "${RED}⚠ Input matrices not found${NC}"
    echo ""
    echo "Please link or copy your preprocessed data:"
    echo "  Option A (symlink): ln -s /path/to/00_b/matrices input/preprocessed/matrices"
    echo "  Option B (copy): cp /path/to/00_b/matrices/*.tsv input/preprocessed/matrices/"
    echo ""
fi