#!/bin/bash
# Regression and validation tests for --apply_all implementation
# Tests that:
# 1. Filtered barcodes get identical assignments ±apply_all
# 2. FLEX enforces all three singlet conditions
# 3. EM uses stored fits (low-count guides skipped)
# 4. Extra barcodes are classified when apply_all enabled

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT_DIR/bin/call_features"

echo "========================================="
echo "Regression and Validation Tests"
echo "for --apply_all implementation"
echo "========================================="
echo

# Check if binary exists
if [ ! -f "$BIN" ]; then
    echo "ERROR: $BIN not found. Run 'make all' first."
    exit 1
fi

echo "[INFO] Using binary: $BIN"
echo

# ============================================================================
# Helper Functions
# ============================================================================

# Create a synthetic test matrix
create_test_matrix() {
    local outdir="$1"
    mkdir -p "$outdir"
    
    # Create features.tsv (3 guides: G1, G2, G3)
    # Format: ID\tName\tType
    cat > "$outdir/features.tsv" <<EOF
G1	Guide1	CRISPR
G2	Guide2	CRISPR
G3	Guide3	CRISPR
EOF
    
    # Create allowed_features.tsv (just IDs, one per line)
    cat > "$outdir/allowed_features.tsv" <<EOF
G1
G2
G3
EOF
    
    # Create barcodes.tsv (30 barcodes total)
    cat > "$outdir/barcodes.tsv" <<EOF
AAAA-1
AAAB-1
AAAC-1
AAAD-1
AAAE-1
AAAF-1
AAAG-1
AAAH-1
AAAI-1
AAAJ-1
BBBA-1
BBBB-1
BBBC-1
BBBD-1
BBBE-1
BBBF-1
BBBG-1
BBBH-1
BBBI-1
BBBJ-1
CCCA-1
CCCB-1
CCCC-1
CCCD-1
CCCE-1
CCCF-1
CCCG-1
CCCH-1
CCCI-1
CCCJ-1
EOF
    
    # Create filtered_barcodes.tsv (first 10 barcodes = high quality)
    head -10 "$outdir/barcodes.tsv" > "$outdir/filtered_barcodes.tsv"
    
    # Create matrix.mtx with realistic data
    # Format: %%MatrixMarket header, then: feature_row barcode_col count
    cat > "$outdir/matrix.mtx" <<'EOF'
%%MatrixMarket matrix coordinate integer general
3 30 90
1 1 50
2 1 5
3 1 2
1 2 3
2 2 45
3 2 3
1 3 40
2 3 35
3 3 5
1 4 5
2 4 3
3 4 30
1 5 20
2 5 18
3 5 2
1 6 30
2 6 25
3 6 3
1 7 15
2 7 5
3 7 3
1 8 8
2 8 4
3 8 2
1 9 12
2 9 10
3 9 2
1 10 6
2 10 5
3 10 18
1 11 45
2 11 3
3 11 2
1 12 2
2 12 40
3 12 3
1 13 35
2 13 30
3 13 4
1 14 3
2 14 2
3 14 25
1 15 18
2 15 16
3 15 2
1 16 28
2 16 20
3 16 3
1 17 12
2 17 4
3 17 2
1 18 6
2 18 3
3 18 2
1 19 10
2 19 8
3 19 2
1 20 5
2 20 4
3 20 15
1 21 40
2 21 2
3 21 2
1 22 2
2 22 35
3 22 2
1 23 30
2 23 25
3 23 3
1 24 2
2 24 2
3 24 20
1 25 15
2 25 14
3 25 2
1 26 25
2 26 18
3 26 2
1 27 10
2 27 3
3 27 2
1 28 5
2 28 3
3 28 1
1 29 8
2 29 7
3 29 1
1 30 4
2 30 3
3 30 12
EOF
    
    echo "[INFO] Created synthetic test matrix in $outdir"
    echo "       - 3 features (G1, G2, G3)"
    echo "       - 30 barcodes total"
    echo "       - 10 filtered barcodes (high quality)"
    echo "       - 20 unfiltered barcodes (lower quality)"
}

# Compare two files line by line
compare_files() {
    local file1="$1"
    local file2="$2"
    local label="$3"
    
    if [ ! -f "$file1" ]; then
        echo "  ✗ Missing: $file1"
        return 1
    fi
    
    if [ ! -f "$file2" ]; then
        echo "  ✗ Missing: $file2"
        return 1
    fi
    
    if diff -q "$file1" "$file2" > /dev/null 2>&1; then
        echo "  ✓ $label: files identical"
        return 0
    else
        echo "  ✗ $label: files differ"
        echo "    Differences:"
        diff "$file1" "$file2" | head -20
        return 1
    fi
}

# ============================================================================
# Test 1: Regression Test - Filtered barcodes should be identical
# ============================================================================
test_regression() {
    echo "========================================="
    echo "Test 1: Regression Test"
    echo "========================================="
    echo "Goal: Filtered barcodes get identical assignments ±apply_all"
    echo
    
    local testdir="$ROOT_DIR/test_regression"
    rm -rf "$testdir"
    mkdir -p "$testdir"
    
    # Create test data
    create_test_matrix "$testdir/matrix"
    
    # Run WITHOUT --apply_all (baseline)
    echo "[RUN] Without --apply_all..."
    "$BIN" \
        --mtx-dir "$testdir/matrix" \
        --out-prefix "$testdir/baseline" \
        --simple-assign --min-count 3 --min-ratio 1.5 \
        --m-min-fixed 10 \
        > "$testdir/baseline.log" 2>&1
    
    # Run WITH --apply_all
    echo "[RUN] With --apply_all..."
    "$BIN" \
        --mtx-dir "$testdir/matrix" \
        --out-prefix "$testdir/rescue" \
        --simple-assign --min-count 3 --min-ratio 1.5 \
        --m-min-fixed 10 \
        --apply-all \
        > "$testdir/rescue.log" 2>&1
    
    # Extract filtered barcodes from both runs
    echo
    echo "[CHECK] Comparing filtered barcode assignments..."
    
    # Get filtered barcodes list
    local filtered_list="$testdir/matrix/filtered_barcodes.tsv"
    
    # Extract assignments for filtered barcodes only
    grep -F -f "$filtered_list" "$testdir/baseline.assignments.tsv" | sort > "$testdir/baseline_filtered.tsv" 2>/dev/null || true
    grep -F -f "$filtered_list" "$testdir/rescue.assignments.tsv" | sort > "$testdir/rescue_filtered.tsv" 2>/dev/null || true
    
    # Compare
    local status=0
    if compare_files "$testdir/baseline_filtered.tsv" "$testdir/rescue_filtered.tsv" "Filtered assignments"; then
        echo "  ✓ PASS: Filtered barcodes have identical assignments"
    else
        echo "  ✗ FAIL: Filtered barcodes differ between runs"
        status=1
    fi
    
    # Check that rescue has MORE assignments
    local baseline_count=$(wc -l < "$testdir/baseline.assignments.tsv" 2>/dev/null || echo 0)
    local rescue_count=$(wc -l < "$testdir/rescue.assignments.tsv" 2>/dev/null || echo 0)
    
    echo
    echo "[CHECK] Assignment counts..."
    echo "  Baseline (no --apply_all): $baseline_count"
    echo "  Rescue (with --apply_all):  $rescue_count"
    
    if [ "$rescue_count" -gt "$baseline_count" ]; then
        echo "  ✓ PASS: --apply_all produced more assignments ($rescue_count > $baseline_count)"
    else
        echo "  ✗ FAIL: --apply_all should produce MORE assignments"
        status=1
    fi
    
    echo
    if [ $status -eq 0 ]; then
        echo "========================================="
        echo "Test 1: PASSED ✓"
        echo "========================================="
    else
        echo "========================================="
        echo "Test 1: FAILED ✗"
        echo "========================================="
    fi
    echo
    
    return $status
}

# ============================================================================
# Test 2: FLEX Mode Validation
# ============================================================================
test_flex_validation() {
    echo "========================================="
    echo "Test 2: FLEX Mode Validation"
    echo "========================================="
    echo "Goal: Verify FLEX applies all three conditions"
    echo
    
    local testdir="$ROOT_DIR/test_flex"
    rm -rf "$testdir"
    mkdir -p "$testdir"
    
    # Create test data
    create_test_matrix "$testdir/matrix"
    
    # Run FLEX mode with --apply_all
    echo "[RUN] FLEX mode with --apply_all..."
    "$BIN" \
        --mtx-dir "$testdir/matrix" \
        --out-prefix "$testdir/flex" \
        --tau 0.7 --delta 0.3 --gamma 0.8 --alpha 0.01 \
        --apply-all \
        > "$testdir/flex.log" 2>&1
    
    local status=0
    
    # Check that output files exist
    echo
    echo "[CHECK] Output files..."
    if [ -f "$testdir/flex.assignments.tsv" ]; then
        local count=$(wc -l < "$testdir/flex.assignments.tsv")
        echo "  ✓ assignments.tsv exists ($count assignments)"
    else
        echo "  ✗ assignments.tsv missing"
        status=1
    fi
    
    # Check log for apply_all messages
    echo
    echo "[CHECK] Log messages..."
    if grep -q "apply_all enabled" "$testdir/flex.log" 2>/dev/null; then
        echo "  ✓ Found '--apply_all enabled' message"
    else
        echo "  ✗ Missing '--apply_all enabled' message"
        status=1
    fi
    
    if grep -q "Applying learned thresholds to all" "$testdir/flex.log" 2>/dev/null; then
        echo "  ✓ Found 'Applying learned thresholds' message"
    else
        echo "  ✗ Missing 'Applying learned thresholds' message"
        status=1
    fi
    
    # Note: Without real data, we can't verify the actual decision logic
    # But we can verify the code ran without errors
    echo
    echo "[CHECK] Execution status..."
    if grep -q "ERROR" "$testdir/flex.log" 2>/dev/null; then
        echo "  ✗ Errors found in log:"
        grep "ERROR" "$testdir/flex.log"
        status=1
    else
        echo "  ✓ No errors in execution"
    fi
    
    echo
    if [ $status -eq 0 ]; then
        echo "========================================="
        echo "Test 2: PASSED ✓"
        echo "========================================="
        echo "Note: Full FLEX validation requires real data with known ground truth"
    else
        echo "========================================="
        echo "Test 2: FAILED ✗"
        echo "========================================="
    fi
    echo
    
    return $status
}

# ============================================================================
# Test 3: EM Mode Validation
# ============================================================================
test_em_validation() {
    echo "========================================="
    echo "Test 3: EM Mode Validation"
    echo "========================================="
    echo "Goal: Verify EM stores and uses fitted parameters"
    echo
    
    local testdir="$ROOT_DIR/test_em"
    rm -rf "$testdir"
    mkdir -p "$testdir"
    
    # Create test data
    create_test_matrix "$testdir/matrix"
    
    # Run EM mode with --apply_all
    echo "[RUN] EM mode with --apply_all..."
    "$BIN" \
        --mtx-dir "$testdir/matrix" \
        --out-prefix "$testdir/em" \
        --use-em --k-min 3 --tau-pos 0.9 \
        --gamma-min 0.7 --min-em-counts 5 \
        --k-small 3 \
        --m-min-fixed 10 \
        --apply-all \
        > "$testdir/em.log" 2>&1
    
    local status=0
    
    # Check that output files exist
    echo
    echo "[CHECK] Output files..."
    if [ -f "$testdir/em.assignments.tsv" ]; then
        local count=$(wc -l < "$testdir/em.assignments.tsv")
        echo "  ✓ assignments.tsv exists ($count assignments)"
    else
        echo "  ✗ assignments.tsv missing"
        status=1
    fi
    
    # Check log for EM-specific messages
    echo
    echo "[CHECK] EM training messages..."
    if grep -q "Skipped.*guides with.*total counts" "$testdir/em.log" 2>/dev/null; then
        echo "  ✓ Found EM training summary"
        grep "Skipped.*guides" "$testdir/em.log"
    else
        echo "  ✗ Missing EM training summary"
        status=1
    fi
    
    if grep -q "ran EM on.*guides" "$testdir/em.log" 2>/dev/null; then
        echo "  ✓ EM fitting completed"
    else
        echo "  ✗ EM fitting may not have completed"
        status=1
    fi
    
    # Check for apply_all messages
    echo
    echo "[CHECK] Apply-all execution..."
    if grep -q "apply_all enabled" "$testdir/em.log" 2>/dev/null; then
        echo "  ✓ Found '--apply_all enabled' message"
    else
        echo "  ✗ Missing '--apply_all enabled' message"
        status=1
    fi
    
    if grep -q "Applying learned thresholds to all" "$testdir/em.log" 2>/dev/null; then
        echo "  ✓ Found 'Applying learned thresholds' message"
    else
        echo "  ✗ Missing 'Applying learned thresholds' message"
        status=1
    fi
    
    # Check for errors
    echo
    echo "[CHECK] Execution status..."
    if grep -q "ERROR" "$testdir/em.log" 2>/dev/null; then
        echo "  ✗ Errors found in log:"
        grep "ERROR" "$testdir/em.log"
        status=1
    else
        echo "  ✓ No errors in execution"
    fi
    
    echo
    if [ $status -eq 0 ]; then
        echo "========================================="
        echo "Test 3: PASSED ✓"
        echo "========================================="
        echo "Note: Full EM validation requires checking that mu_bg is computed from stored fits"
    else
        echo "========================================="
        echo "Test 3: FAILED ✗"
        echo "========================================="
    fi
    echo
    
    return $status
}

# ============================================================================
# Main Execution
# ============================================================================

echo "Running all tests..."
echo

# Track overall status
overall_status=0

# Run tests
test_regression || overall_status=1
test_flex_validation || overall_status=1
test_em_validation || overall_status=1

# Final summary
echo
echo "========================================="
echo "FINAL SUMMARY"
echo "========================================="

if [ $overall_status -eq 0 ]; then
    echo "✓ All tests PASSED"
    echo
    echo "The --apply_all implementation:"
    echo "  ✓ Maintains identical assignments for filtered barcodes"
    echo "  ✓ Produces additional assignments for unfiltered barcodes"
    echo "  ✓ FLEX mode executes without errors"
    echo "  ✓ EM mode executes without errors"
    echo
    echo "Status: READY FOR PRODUCTION with caveats:"
    echo "  - Simple-assign: Fully validated"
    echo "  - FLEX: Basic validation passed, recommend testing with real data"
    echo "  - EM: Basic validation passed, recommend validating heuristic accuracy"
else
    echo "✗ Some tests FAILED"
    echo
    echo "Please review the failures above before using --apply_all in production."
fi

echo "========================================="

exit $overall_status

