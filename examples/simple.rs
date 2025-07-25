use p533::*;

fn main() {
    println!("P533 HF Propagation Prediction - Simple Example");
    println!("===============================================");

    // Create a prediction for Dublin to New York on 14 MHz
    let path = create_path_prediction(
        53.3, -6.2, // Dublin: 53.3°N, 6.2°W
        40.7, -74.0, // New York: 40.7°N, 74.0°W
        2024,  // Year
        5,     // Month (May)
        12,    // Hour (12:00 UTC)
        100,   // Solar Sunspot Number
        14.0,  // Frequency (MHz)
        100.0, // TX Power (W)
        3.0,   // TX Antenna Gain (dBi)
        2.0,   // RX Antenna Gain (dBi)
    );

    // Display the results
    println!("\nPath: Dublin to New York");
    println!("Distance: {:.1} km", path.distance);
    println!("Frequency: {:.1} MHz", path.frequency);
    println!(
        "Date/Time: {}-{:02} {:02}:00 UTC",
        path.year, path.month, path.hour
    );
    println!();

    println!("Maximum Usable Frequencies:");
    println!("  MUF50 (50%): {:.2} MHz", path.muf50);
    println!("  MUF90 (90%): {:.2} MHz", path.muf90);
    println!("  MUF10 (10%): {:.2} MHz", path.muf10);
    println!();

    println!("Signal Performance:");
    println!("  Field Strength: {:.2} dB(1μV/m)", path.es);
    println!("  Received Power: {:.2} dBW", path.pr);
    println!("  Signal-to-Noise Ratio: {:.2} dB", path.snr);
    println!("  Circuit Reliability: {:.1}%", path.bcr);
    println!();

    // Provide interpretation
    if path.frequency <= path.muf90 {
        println!("✓ Excellent propagation conditions");
    } else if path.frequency <= path.muf50 {
        println!("✓ Good propagation conditions");
    } else if path.frequency <= path.muf10 {
        println!("⚠ Marginal propagation conditions");
    } else {
        println!("✗ Poor propagation conditions");
    }

    if path.bcr >= 90.0 {
        println!("✓ Excellent circuit reliability");
    } else if path.bcr >= 70.0 {
        println!("✓ Good circuit reliability");
    } else if path.bcr >= 50.0 {
        println!("⚠ Moderate circuit reliability");
    } else {
        println!("✗ Poor circuit reliability");
    }
}
