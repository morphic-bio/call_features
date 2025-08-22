#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# compareAllMethods.sh – run all demultiplexing methods and compare results
#
# Usage:
#   ./compareAllMethods.sh
#
# This script runs all available methods in testCallPerturb.sh and creates
# a comprehensive comparison table of the results.
# ---------------------------------------------------------------------------

set -euo pipefail

# Configuration
SC_ID="30_KO_ES"
OUTPUT_DIR="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
COMPARISON_FILE="${OUTPUT_DIR}/method_comparison_${TIMESTAMP}.tsv"
LOG_DIR="${OUTPUT_DIR}/logs_${TIMESTAMP}"

# All available methods
METHODS=(flex em em-fixed em-noncells simple otsu quantile model3)

echo "🧪 Testing all demultiplexing methods..." >&2
echo "========================================" >&2
echo "Sample ID: $SC_ID" >&2
echo "Output directory: $OUTPUT_DIR" >&2
echo "Comparison file: $COMPARISON_FILE" >&2
echo "" >&2

# Create log directory
mkdir -p "$LOG_DIR"

# Function to extract statistics from program output
extract_stats_from_output() {
    local method="$1"
    local log_file="$2"
    
    echo "Extracting stats for $method from $log_file" >&2
    
    # Initialize variables
    local n_singlet=0
    local n_doublet=0
    local n_ambig=0
    local n_low=0
    local total=0
    
    # Look for the QC output lines in the log
    if [[ -f "$log_file" ]]; then
        # Extract counts from lines like "Singlet: 1234 (0.567)"
        n_singlet=$(grep "Singlet:" "$log_file" | grep -o "[0-9]*" | head -1 2>/dev/null || echo "0")
        n_doublet=$(grep "Doublet:" "$log_file" | grep -o "[0-9]*" | head -1 2>/dev/null || echo "0")
        n_ambig=$(grep "Ambiguous:" "$log_file" | grep -o "[0-9]*" | head -1 2>/dev/null || echo "0")
        n_low=$(grep "Low support:" "$log_file" | grep -o "[0-9]*" | head -1 2>/dev/null || echo "0")
        
        # Alternative: look for "Evaluated (filtered∩MTX): XXXX"
        if [[ $n_singlet -eq 0 && $n_doublet -eq 0 && $n_ambig -eq 0 && $n_low -eq 0 ]]; then
            total=$(grep "Evaluated.*MTX" "$log_file" | grep -o "[0-9]*" | tail -1 2>/dev/null || echo "0")
            # If we have total but no breakdown, try different extraction patterns
            if [[ $total -gt 0 ]]; then
                # Try to extract from summary lines
                local singlet_line=$(grep -A 5 "Evaluated.*MTX" "$log_file" | grep -i "singlet" | head -1)
                local doublet_line=$(grep -A 5 "Evaluated.*MTX" "$log_file" | grep -i "doublet" | head -1)
                local ambig_line=$(grep -A 5 "Evaluated.*MTX" "$log_file" | grep -i "ambiguous" | head -1)
                local low_line=$(grep -A 5 "Evaluated.*MTX" "$log_file" | grep -i "low" | head -1)
                
                [[ -n "$singlet_line" ]] && n_singlet=$(echo "$singlet_line" | grep -o "[0-9]*" | head -1 || echo "0")
                [[ -n "$doublet_line" ]] && n_doublet=$(echo "$doublet_line" | grep -o "[0-9]*" | head -1 || echo "0")
                [[ -n "$ambig_line" ]] && n_ambig=$(echo "$ambig_line" | grep -o "[0-9]*" | head -1 || echo "0")
                [[ -n "$low_line" ]] && n_low=$(echo "$low_line" | grep -o "[0-9]*" | head -1 || echo "0")
            fi
        fi
    fi
    
    # Calculate total if not found directly
    if [[ $total -eq 0 ]]; then
        total=$((n_singlet + n_doublet + n_ambig + n_low))
    fi
    
    # Calculate percentages
    local singlet_pct=$(awk "BEGIN {printf \"%.3f\", ($total > 0) ? $n_singlet/$total : 0}")
    local doublet_pct=$(awk "BEGIN {printf \"%.3f\", ($total > 0) ? $n_doublet/$total : 0}")
    local ambig_pct=$(awk "BEGIN {printf \"%.3f\", ($total > 0) ? $n_ambig/$total : 0}")
    local low_pct=$(awk "BEGIN {printf \"%.3f\", ($total > 0) ? $n_low/$total : 0}")
    
    echo "  Extracted: Total=$total, Singlet=$n_singlet, Doublet=$n_doublet, Ambig=$n_ambig, Low=$n_low" >&2
    
    # Return the stats as a tab-separated line (to stdout for capture)
    echo -e "$method\t$total\t$n_singlet\t$singlet_pct\t$n_doublet\t$doublet_pct\t$n_ambig\t$ambig_pct\t$n_low\t$low_pct"
}

# Function to run a single method
run_method() {
    local method="$1"
    local log_file="$LOG_DIR/${method}.log"
    
    echo "🔄 Running $method method..." >&2
    echo "   Log file: $log_file" >&2
    
    # Run the method and capture all output
    if timeout 30m ./scripts/testCallPerturb.sh "$method" &> "$log_file"; then
        echo "   ✅ $method completed successfully" >&2
        
        # Extract statistics from the log
        extract_stats_from_output "$method" "$log_file"
        return 0
    else
        echo "   ❌ $method failed or timed out" >&2
        # Return empty stats
        echo -e "$method\t0\t0\t0.000\t0\t0.000\t0\t0.000\t0\t0.000"
        return 1
    fi
}

# Initialize comparison file with header
echo -e "Method\tTotal_Cells\tSinglet_Count\tSinglet_Pct\tDoublet_Count\tDoublet_Pct\tAmbiguous_Count\tAmbiguous_Pct\tLow_Support_Count\tLow_Support_Pct" > "$COMPARISON_FILE"

# Track success/failure
declare -a SUCCESSFUL_METHODS=()
declare -a FAILED_METHODS=()

# Run each method and collect results
for method in "${METHODS[@]}"; do
    echo "" >&2
    echo "==========================================" >&2
    echo "Testing method: $method" >&2
    echo "==========================================" >&2
    
    # Run method and capture stats
    if stats_line=$(run_method "$method"); then
        # Add to comparison file
        echo "$stats_line" >> "$COMPARISON_FILE"
        SUCCESSFUL_METHODS+=("$method")
    else
        # Add failed entry to comparison file  
        echo -e "$method\t0\t0\t0.000\t0\t0.000\t0\t0.000\t0\t0.000" >> "$COMPARISON_FILE"
        FAILED_METHODS+=("$method")
    fi
done

echo "" >&2
echo "==========================================" >&2
echo "📊 COMPARISON RESULTS" >&2
echo "==========================================" >&2

# Display the comparison table
if [[ -f "$COMPARISON_FILE" ]]; then
    echo "" >&2
    echo "Results table:" >&2
    column -t -s $'\t' "$COMPARISON_FILE"
    echo ""
    echo "Full results saved to: $COMPARISON_FILE"
else
    echo "❌ Error: Could not create comparison file" >&2
    exit 1
fi

# Summary
echo "" >&2
echo "📈 SUMMARY" >&2
echo "==========" >&2
echo "Total methods tested: ${#METHODS[@]}" >&2
echo "Successful: ${#SUCCESSFUL_METHODS[@]} (${SUCCESSFUL_METHODS[*]})" >&2
echo "Failed: ${#FAILED_METHODS[@]} (${FAILED_METHODS[*]})" >&2
echo "" >&2

# Method-specific insights
if [[ ${#SUCCESSFUL_METHODS[@]} -gt 0 ]]; then
    echo "💡 Method Insights:" >&2
    echo "  Original methods:" >&2
    [[ " ${SUCCESSFUL_METHODS[*]} " =~ " flex " ]] && echo "    - flex: Traditional binomial test approach" >&2
    [[ " ${SUCCESSFUL_METHODS[*]} " =~ " em " ]] && echo "    - em: EM algorithm with dispersion updates" >&2
    [[ " ${SUCCESSFUL_METHODS[*]} " =~ " em-fixed " ]] && echo "    - em-fixed: EM algorithm with fixed dispersion" >&2
    [[ " ${SUCCESSFUL_METHODS[*]} " =~ " em-noncells " ]] && echo "    - em-noncells: EM using non-cells for ambient estimation" >&2
    echo "" >&2
    echo "  Cell-derived M_min methods:" >&2
    [[ " ${SUCCESSFUL_METHODS[*]} " =~ " simple " ]] && echo "    - simple: Ratio-based threshold (simple and fast)" >&2
    [[ " ${SUCCESSFUL_METHODS[*]} " =~ " em-cells " ]] && echo "    - em-cells: EM-cells method (cell-derived M_min)" >&2
fi

echo "" >&2
echo "📁 Individual method logs saved in: $LOG_DIR" >&2
echo "" >&2

if [[ ${#FAILED_METHODS[@]} -eq 0 ]]; then
    echo "✨ All methods completed successfully!" >&2
else
    echo "⚠️  Some methods failed. Check individual log files for details." >&2
fi

echo "" >&2
echo "🎯 Next steps:" >&2
echo "  1. Review the comparison table to identify the best-performing method" >&2
echo "  2. Check individual log files for detailed method-specific output" >&2
echo "  3. Consider the trade-offs between accuracy and computational cost" >&2
echo "" >&2
