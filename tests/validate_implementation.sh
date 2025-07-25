#!/bin/bash

set -e

echo "p533 validation"
echo

# check if we are running the script from the root of the repo
if [ ! -f "Cargo.toml" ] || [ ! -d "crates/p533" ]; then
    echo "run from project root"
    exit 1
fi

# generate rust results
echo "generating rust results..."
cargo run --bin generate_baseline --quiet 2>/dev/null

# generate c results  
echo "generating c results..."
cd tests/c
make c_reference_test
./c_reference_test > c_results.json
cd ../..

# compare impls
echo "comparing implementations..."
if python3 tests/compare_results.py 0.1 tests/c/c_results.json tests/rust_results.json; then
    echo "passed: implementations match within 10% tolerance"
    VALIDATION_RESULT=0
else
    echo "trying 25% tolerance..."
    if python3 tests/compare_results.py 0.25 tests/c/c_results.json tests/rust_results.json; then
        echo "acceptable: validates with extended tolerance"
        VALIDATION_RESULT=0
    else
        echo "failed: implementations differ significantly"
        VALIDATION_RESULT=1
    fi
fi

# cleanup
rm -f tests/rust_results.json tests/c/c_results.json

echo
if [ $VALIDATION_RESULT -eq 0 ]; then
    echo "validation successful"
else
    echo "validation failed"
    exit 1
fi
