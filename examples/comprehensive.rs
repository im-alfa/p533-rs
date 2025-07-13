use p533::*;

fn main() {
    println!("P533 HF Propagation Prediction - Comprehensive Example");
    println!("======================================================");
    
    // Example 1: Basic HF propagation prediction
    println!("\n1. Basic HF Propagation Prediction");
    println!("----------------------------------");
    
    // Dublin to New York path
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
    
    println!("Path: Dublin to New York");
    println!("Distance: {:.1} km", path.distance);
    println!("Frequency: {:.1} MHz", path.frequency);
    println!("Date/Time: {}-{:02} {:02}:00 UTC", path.year, path.month, path.hour);
    println!("Solar Activity (SSN): {}", path.ssn);
    println!();
    
    // Display key results
    println!("Propagation Results:");
    println!("  MUF50 (Median): {:.2} MHz", path.muf50);
    println!("  MUF90 (90%):    {:.2} MHz", path.muf90);
    println!("  MUF10 (10%):    {:.2} MHz", path.muf10);
    println!("  OPMUF:          {:.2} MHz", path.opmuf);
    println!();
    
    println!("Signal Strength:");
    println!("  Field Strength: {:.2} dB(1μV/m)", path.es);
    println!("  Received Power: {:.2} dBW", path.pr);
    println!("  Signal-to-Noise Ratio: {:.2} dB", path.snr);
    println!();
    
    println!("Circuit Performance:");
    println!("  Basic Circuit Reliability: {:.1}%", path.bcr);
    println!("  Required SNR: {:.1} dB", path.snrr);
    println!();
    
    // Example 2: Multiple frequency analysis
    println!("\n2. Frequency Analysis");
    println!("---------------------");
    
    let frequencies = vec![7.0, 10.0, 14.0, 18.0, 21.0, 28.0];
    println!("Testing different frequencies for Dublin to New York path:");
    println!("Freq (MHz) | MUF50 | BCR (%) | SNR (dB) | Status");
    println!("-----------|-------|---------|----------|--------");
    
    for freq in frequencies {
        let path = create_path_prediction(
            53.3, -6.2, 40.7, -74.0,
            2024, 5, 12, 100,
            freq, 100.0, 3.0, 2.0
        );
        
        let status = if freq <= path.muf90 {
            "Excellent"
        } else if freq <= path.muf50 {
            "Good"
        } else if freq <= path.muf10 {
            "Marginal"
        } else {
            "Poor"
        };
        
        println!("{:8.1}   | {:5.1} | {:7.1} | {:8.1} | {}", 
                freq, path.muf50, path.bcr, path.snr, status);
    }
    
    // Example 3: Time-based analysis
    println!("\n3. Time-Based Analysis (24-hour)");
    println!("--------------------------------");
    
    println!("Analyzing 14 MHz propagation over 24 hours:");
    println!("Hour (UTC) | MUF50 | BCR (%) | SNR (dB)");
    println!("-----------|-------|---------|----------");
    
    for hour in (0..24).step_by(3) {
        let path = create_path_prediction(
            53.3, -6.2, 40.7, -74.0,
            2024, 5, hour, 100,
            14.0, 100.0, 3.0, 2.0
        );
        
        println!("{:8}   | {:5.1} | {:7.1} | {:8.1}", 
                hour, path.muf50, path.bcr, path.snr);
    }
    
    // Example 4: Different paths
    println!("\n4. Different Path Examples");
    println!("-------------------------");
    
    let paths = vec![
        ("Dublin to New York", 53.3, -6.2, 40.7, -74.0),
        ("Los Angeles to Tokyo", 34.0, -118.0, 35.7, 139.7),
        ("Sydney to San Francisco", -33.9, 151.2, 37.8, -122.4),
        ("Cape Town to Buenos Aires", -33.9, 18.4, -34.6, -58.4),
    ];
    
    println!("14 MHz propagation for different paths:");
    println!("Path                      | Distance | MUF50 | BCR (%) | SNR (dB)");
    println!("--------------------------|----------|-------|---------|----------");
    
    for (name, tx_lat, tx_lng, rx_lat, rx_lng) in paths {
        let path = create_path_prediction(
            tx_lat, tx_lng, rx_lat, rx_lng,
            2024, 5, 12, 100,
            14.0, 100.0, 3.0, 2.0
        );
        
        println!("{:24}  | {:6.0} km | {:5.1} | {:7.1} | {:8.1}", 
                name, path.distance, path.muf50, path.bcr, path.snr);
    }
    
    // Example 5: Solar activity impact
    println!("\n5. Solar Activity Impact");
    println!("-----------------------");
    
    let solar_conditions = vec![
        ("Solar Minimum", 10),
        ("Low Activity", 50),
        ("Moderate Activity", 100),
        ("High Activity", 150),
        ("Solar Maximum", 200),
    ];
    
    println!("Impact of solar activity on 14 MHz (Dublin to New York):");
    println!("Condition         | SSN | MUF50 | BCR (%) | SNR (dB)");
    println!("------------------|-----|-------|---------|----------");
    
    for (condition, ssn) in solar_conditions {
        let path = create_path_prediction(
            53.3, -6.2, 40.7, -74.0,
            2024, 5, 12, ssn,
            14.0, 100.0, 3.0, 2.0
        );
        
        println!("{:16}  | {:3} | {:5.1} | {:7.1} | {:8.1}", 
                condition, ssn, path.muf50, path.bcr, path.snr);
    }
    
    // Example 6: Using the convenience functions
    println!("\n6. Convenience Functions");
    println!("-----------------------");
    
    // Get MUF information
    let (muf50, muf90, muf10, opmuf, distance) = get_muf_info(
        53.3, -6.2, 40.7, -74.0,
        2024, 5, 12, 100
    );
    
    println!("MUF Information for Dublin to New York:");
    println!("  Distance: {:.1} km", distance);
    println!("  MUF50: {:.2} MHz", muf50);
    println!("  MUF90: {:.2} MHz", muf90);
    println!("  MUF10: {:.2} MHz", muf10);
    println!("  OPMUF: {:.2} MHz", opmuf);
    println!();
    
    // Calculate circuit reliability
    let reliability = calculate_circuit_reliability(
        53.3, -6.2, 40.7, -74.0,
        2024, 5, 12, 100,
        14.0, 20.0, // 14 MHz, 20 dB required SNR
        100.0, 3.0, 2.0
    );
    
    println!("Circuit Reliability for 14 MHz:");
    println!("  Required SNR: 20 dB");
    println!("  Estimated Reliability: {:.1}%", reliability);
    
    println!("\n======================================================");
    println!("Example completed successfully!");
    println!("For more detailed analysis, examine the PathData structure");
    println!("returned by create_path_prediction().");
}
