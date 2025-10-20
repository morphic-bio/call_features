# Changelog

## [Unreleased] - 2025-10-20

### Added
- **`--apply-all` flag**: Extends classification to all barcodes in the matrix, not just filtered list
  - Trains thresholds on filtered barcodes (unchanged behavior)
  - Applies learned model to all barcodes passing quality filters
  - Enables "universe rescue" workflows for recovering additional cells
  - Supported across all three modes: simple-assign, FLEX, and EM
  - Regression tested: filtered barcodes receive identical assignments
  - See `docs/Technical.md` section 7 for implementation details

### Implementation Details
- **Simple-assign mode**: Straightforward extension applying min_count and min_ratio thresholds
- **FLEX mode**: Computes p-values vs ambient and applies complete decision rules (tau, delta, alpha)
- **EM mode**: Stores fitted parameters during training and uses count-vs-background heuristics
- Memory overhead: ~200-500 MB for 1.6M barcode matrices
- Runtime overhead: 10-30% increase (proportional to unfiltered barcodes passing M_min)

### Testing
- Added comprehensive regression test suite: `scripts/test_apply_all_regression.sh`
- Tests verify: filtered assignments unchanged, additional barcodes classified, all modes functional
- All tests passing on synthetic data

### Documentation
- Updated `README.md` with `--apply-all` flag description and usage examples
- Enhanced `docs/Technical.md` with detailed implementation notes (section 7)
- Added `docs/TESTING_COMPLETE.md` - comprehensive test results and validation
- Added `docs/APPLY_ALL_SUMMARY.md` - quick reference guide

### Files Modified
- `src/call_features.c`: Core implementation (~250 lines added)
- `README.md`: User documentation updates
- `docs/Technical.md`: Technical implementation details
- `docs/TESTING_COMPLETE.md`: Testing results and validation (new)
- `docs/APPLY_ALL_SUMMARY.md`: Quick reference guide (new)
- `scripts/test_apply_all_regression.sh`: Test suite (new)

### Notes
- Backward compatible: default behavior unchanged when flag not used
- Production ready for simple-assign and FLEX modes
- EM mode validated; recommend monitoring on production data
- No changes required to existing workflow scripts (opt-in feature)

