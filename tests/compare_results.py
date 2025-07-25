#!/usr/bin/env python3

import json
import sys

def load_results(filename):
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        sys.exit(1)

def compare_values(a, b):
    if a == 0 and b == 0:
        return 0.0
    if a == 0:
        return 1.0
    return abs((b - a) / a)

def compare_results(c_results, rust_results, tolerance):
    failures = []
    
    # Create lookup dict for rust results
    rust_dict = {r['name']: r for r in rust_results}
    
    for c_result in c_results:
        name = c_result['name']
        
        if 'error' in c_result:
            failures.append(f"{name}: C calculation failed")
            continue
            
        if name not in rust_dict:
            failures.append(f"{name}: missing in Rust results")
            continue
            
        rust_result = rust_dict[name]
        
        # Compare key fields
        fields = ['distance', 'muf50', 'bcr', 'snr']
        
        for field in fields:
            if field in c_result and field in rust_result:
                diff = compare_values(c_result[field], rust_result[field])
                if diff > tolerance:
                    failures.append(f"{name}.{field}: {diff*100:.1f}% (C={c_result[field]:.3f}, Rust={rust_result[field]:.3f})")
    
    return failures

def main():
    if len(sys.argv) != 4:
        print("Usage: compare_results.py <tolerance> <c_results.json> <rust_results.json>")
        sys.exit(1)
    
    tolerance = float(sys.argv[1])
    c_file = sys.argv[2]
    rust_file = sys.argv[3]
    
    print(f"Comparing with {tolerance*100:.1f}% tolerance")
    
    c_results = load_results(c_file)
    rust_results = load_results(rust_file)
    
    failures = compare_results(c_results, rust_results, tolerance)
    
    if not failures:
        print(f"PASS: All {len(c_results)} scenarios within tolerance")
        return 0
    else:
        print(f"FAIL: {len(failures)} issues found:")
        for failure in failures:
            print(f"  {failure}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
