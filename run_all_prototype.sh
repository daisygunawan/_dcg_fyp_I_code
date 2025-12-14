#!/bin/bash
# run_all_prototype.sh
# Execute complete prototype workflow (Steps 0-4)

set -e  # Exit on error

echo "========================================"
echo "TCGA-BRCA Network Rewiring Prototype"
echo "Complete Workflow Execution"
echo "========================================"
echo ""

# Check if virtual environment is activated
if [[ -z "$VIRTUAL_ENV" ]]; then
    echo "Warning: Virtual environment not detected"
    echo "   Consider activating: source venv_proto/bin/activate"
    echo ""
    read -p "Continue anyway? (y/n): " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Start timing
start_time=$(date +%s)

# Step 0: Create subset matrices
echo "Step 0: Creating subset matrices..."
python scripts/proto_00_create_subset.py
if [ $? -ne 0 ]; then
    echo "Step 0 failed"
    exit 1
fi
echo "Step 0 complete"
echo ""

# Step 1: Data preparation & correlation
echo "Step 1: Computing correlation matrices..."
python scripts/proto_01_prepare_data.py
if [ $? -ne 0 ]; then
    echo "Step 1 failed"
    exit 1
fi
echo "Step 1 complete"
echo ""

# Step 2: Differential analysis
echo "Step 2: Differential co-expression analysis..."
python scripts/proto_02_differential.py
if [ $? -ne 0 ]; then
    echo "Step 2 failed"
    exit 1
fi
echo "Step 2 complete"
echo ""

# Step 3: Build subgraphs
echo "Step 3: Building hub subgraphs..."
python scripts/proto_03_build_subgraphs.py
if [ $? -ne 0 ]; then
    echo "Step 3 failed"
    exit 1
fi
echo "Step 3 complete"
echo ""

# Step 4: Visualize
echo "Step 4: Generating visualizations..."
python scripts/proto_04_visualize.py
if [ $? -ne 0 ]; then
    echo "Step 4 failed"
    exit 1
fi
echo "Step 4 complete"
echo ""

# Step 5: Diagnostic analysis
echo "Step 5: Running diagnostic analysis..."
python scripts/proto_05_diagnostic.py
if [ $? -ne 0 ]; then
    echo "Step 5 failed"
    exit 1
fi
echo "Step 5 complete"
echo ""

# Calculate total time
end_time=$(date +%s)
total_time=$((end_time - start_time))
minutes=$((total_time / 60))
seconds=$((total_time % 60))

# Summary
echo "========================================"
echo "PROTOTYPE COMPLETE!"
echo "========================================"
echo "Total execution time: ${minutes}m ${seconds}s"
echo ""
echo "Results available in:"
echo "   • output/differential_results/"
echo "   • output/visualizations/"
echo ""
echo "Key files for thesis:"
echo "   • output/visualizations/04_hub_comparison_grid.png"
echo "   • output/differential_results/summary_stats.json"
echo ""
echo "Next: Review visualizations and proceed to full pipeline"
echo "========================================"