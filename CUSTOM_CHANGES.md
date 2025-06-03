# Custom Modifications to NIST SP800-90B

## Compression-Only Test Tool

This fork adds a fast compression-only entropy assessment tool for preliminary testing.

### New Files:
- `cpp/compression_only_main.cpp` - Main program for compression-only testing
- Modified `cpp/Makefile` - Added compression_only target

### Usage:
```bash
make compression_only
./ea_compression_only your_data_file.bin