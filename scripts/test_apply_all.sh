#!/bin/bash
# Test script for --apply_all functionality
# Verifies that the flag is recognized and help text is displayed correctly

set -e

echo "========================================="
echo "Testing --apply_all flag implementation"
echo "========================================="
echo

# Test 1: Check binary exists
echo "[Test 1] Checking if binary exists..."
if [ ! -f "bin/call_features" ]; then
    echo "ERROR: bin/call_features not found. Run 'make all' first."
    exit 1
fi
echo "✓ Binary found"
echo

# Test 2: Check help text includes --apply-all
echo "[Test 2] Verifying --apply-all in help text..."
if ! ./bin/call_features --help | grep -q "apply-all"; then
    echo "ERROR: --apply-all not found in help text"
    exit 1
fi
echo "✓ Help text includes --apply-all flag"
echo

# Test 3: Check detailed help description
echo "[Test 3] Checking help description..."
./bin/call_features --help | grep -A 1 "apply-all"
echo "✓ Help description displayed"
echo

# Test 4: Verify flag is accepted (without required args, should fail gracefully)
echo "[Test 4] Testing flag acceptance..."
if ./bin/call_features --apply-all 2>&1 | grep -q "ERROR.*mtx-dir.*out-prefix"; then
    echo "✓ Flag accepted (correctly requires mtx-dir and out-prefix)"
else
    echo "✗ Unexpected behavior when using --apply-all alone"
    exit 1
fi
echo

echo "========================================="
echo "All basic tests passed!"
echo "========================================="
echo
echo "Next steps:"
echo "1. Test with actual data: provide --mtx-dir, --out-prefix, and --apply-all"
echo "2. Compare output with and without --apply-all on same dataset"
echo "3. Verify extra barcodes are classified in apply-all mode"
echo "4. Check memory usage on large datasets"
echo
echo "Example usage:"
echo "  ./bin/call_features \\"
echo "    --mtx-dir /path/to/matrix \\"
echo "    --out-prefix results/test \\"
echo "    --apply-all \\"
echo "    --tau 0.8 --delta 0.4 --gamma 0.9"

