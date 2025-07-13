# P533 HF Propagation Prediction - Rust Implementation

This is a Rust implementation of ITU-R P.533-12 for HF propagation prediction, integrating ITU-R P.372-16 for atmospheric noise calculations.

It's heavily inspired by the [C implementation](https://github.com/ITU-R-Study-Group-3/ITU-R-HF) of the same recommendation, created by a group of researchers.

## Features

- **Complete P.533-12 implementation**: Maximum Usable Frequency (MUF), field strength, and circuit reliability calculations
- **P.372-16 atmospheric noise integration**: Accurate noise modeling for different terrain types and conditions
- **Multiple propagation modes**: E-layer and F2-layer propagation mode analysis

## Crates
- `p372` - **Fully implemented** - Radio noise characteristics and calculation methods
- `p533` - **Fully implemented** - HF propagation prediction method

## Quick Start

### Basic Usage

```rust
use p533::*;

// Create a path prediction from Dublin to New York
let path = create_path_prediction(
    53.3, -6.2,   // Dublin: 53.3°N, 6.2°W
    40.7, -74.0,  // New York: 40.7°N, 74.0°W
    2024,         // Year
    5,            // Month (May)
    12,           // Hour (12:00 UTC)
    100,          // Solar Sunspot Number
    14.0,         // Frequency (MHz)
    100.0,        // TX Power (W)
    3.0,          // TX Antenna Gain (dBi)
    2.0           // RX Antenna Gain (dBi)
);

// Access the results
println!("MUF50: {:.2} MHz", path.muf50);
println!("Field Strength: {:.2} dB(1μV/m)", path.es);
println!("Circuit Reliability: {:.1}%", path.bcr);
```

### Convenience Functions

```rust
// Get MUF information only
let (muf50, muf90, muf10, opmuf, distance) = get_muf_info(
    53.3, -6.2, 40.7, -74.0,
    2024, 5, 12, 100
);

// Calculate circuit reliability
let reliability = calculate_circuit_reliability(
    53.3, -6.2, 40.7, -74.0,
    2024, 5, 12, 100,
    14.0, 20.0, // 14 MHz, 20 dB required SNR
    100.0, 3.0, 2.0
);
```

## Examples

The crate includes comprehensive examples:

### 1. Simple Example
```bash
cargo run -p p533 --example simple
```
Basic usage with a single path prediction and interpretation.

### 2. Comprehensive Example
```bash
cargo run -p p533 --example comprehensive
```
Multiple analysis types including:
- Frequency analysis across different bands
- 24-hour time-based analysis
- Different path examples
- Solar activity impact analysis

### 3. Technical Example
```bash
cargo run -p p533 --example technical
```
Detailed technical analysis showing:
- Control point information
- Propagation mode details
- Noise analysis
- Complete parameter breakdown

## What is P.372-16 & P.533-12?
- P.372-16 is a recommendation from ITU-R that describes a method to predict the radio noise floor in the HF band. This implementation shouldn't be used for any other purpose than study and research.
- P.533-12 is a recommendation from ITU-R that describes a method to predict HF propagation characteristics. This implementation shouldn't be used for any other purpose than study and research.

PDFs can be found on their official library:
- [https://www.itu.int/rec/R-REC-P/en](https://www.itu.int/rec/R-REC-P/en)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
