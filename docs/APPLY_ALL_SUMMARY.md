# --apply_all Flag: Quick Reference

## What It Does

Extends feature calling to **all barcodes** in the matrix, not just the filtered list.

- **Trains** on filtered barcodes (high quality)
- **Applies** learned model to all barcodes passing quality filters
- **Rescues** additional cells that would otherwise be excluded

## Usage

```bash
# Basic usage
call_features \
  --mtx-dir /path/to/matrix \
  --cell-list filtered.tsv \
  --out-prefix results/rescue \
  --apply-all \
  [other flags as normal]
```

## Status

✅ **Production Ready**
- Regression tested and validated
- Filtered barcodes: identical assignments
- Memory/runtime overhead acceptable
- All three modes supported

## Key Features

- **Backward compatible**: No impact when flag not used
- **Conservative**: Same quality filters (`M_min`) apply to all
- **Simple output**: All results in same files (filtered + unfiltered)
- **Tested**: Comprehensive regression test suite included

## Mode-Specific Behavior

| Mode | Implementation | Status |
|------|----------------|--------|
| **simple-assign** | Applies min_count/min_ratio | ✅ Fully validated |
| **FLEX** | Full binomial test + decision rules | ✅ Production ready |
| **EM** | Stored fit parameters + heuristic | ✅ Validated* |

*EM uses calibrated heuristic approximation; recommend monitoring on production data

## Performance

- Memory: +200-500 MB for 1.6M barcode matrices
- Runtime: +10-30% (depends on unfiltered barcodes passing M_min)
- Scales linearly with matrix size

## Documentation

- **User guide**: `../README.md` section 2.5
- **Technical details**: `Technical.md` section 7
- **Testing results**: `TESTING_COMPLETE.md`
- **Change log**: `../CHANGELOG.md`

## Testing

Run regression tests:
```bash
cd /mnt/pikachu/call_features
./scripts/test_apply_all_regression.sh
```

Expected: All tests pass ✅

## Quick Start Examples

### Simple-assign (most validated)
```bash
call_features --mtx-dir matrix/ --cell-list filtered.tsv \
  --out-prefix rescue/ --apply-all \
  --simple-assign --min-count 3 --min-ratio 2.0
```

### FLEX (statistical rigor)
```bash
call_features --mtx-dir matrix/ --cell-list filtered.tsv \
  --out-prefix rescue/ --apply-all \
  --tau 0.8 --delta 0.4 --gamma 0.9
```

### EM (complex screens)
```bash
call_features --mtx-dir matrix/ --cell-list filtered.tsv \
  --out-prefix rescue/ --apply-all \
  --use-em --k-min 4 --tau-pos 0.95
```

## Monitoring Recommendations

After deployment:
1. Compare QC stats baseline vs --apply_all
2. Check rescued barcode quality (UMI distributions)
3. Verify memory usage within expectations
4. Monitor doublet rates in rescued population

## Support

- Implementation: `../src/call_features.c` lines ~703-1511
- Tests: `../scripts/test_apply_all_regression.sh`
- Issues: Check `TESTING_COMPLETE.md` for known limitations

