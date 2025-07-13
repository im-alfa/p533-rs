use p533::*;

fn main() {
    println!("P533 HF Propagation Prediction - Technical Example");
    println!("==================================================");
    
    // Create a prediction for Dublin to New York on 14 MHz
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
    
    println!("Basic Path Information:");
    println!("======================");
    println!("Path distance: {:.1} km", path.distance);
    println!("Frequency: {:.1} MHz", path.frequency);
    println!("TX Power: {:.1} dBW", path.txpower);
    println!("Bandwidth: {:.0} Hz", path.bw);
    println!();
    
    println!("Ionospheric Conditions:");
    println!("=======================");
    println!("Year: {}, Month: {}, Hour: {} UTC", path.year, path.month, path.hour);
    println!("Solar Sunspot Number: {}", path.ssn);
    println!();
    
    println!("Control Points:");
    println!("===============");
    for (i, cp) in path.cp.iter().enumerate().take(5) {
        if cp.l.lat != 0.0 || cp.l.lng != 0.0 {
            println!("  CP {}: Lat={:.2}°, Lng={:.2}°", i + 1, cp.l.lat, cp.l.lng);
            println!("        foF2={:.2} MHz, foE={:.2} MHz, M3kF2={:.2}", 
                     cp.fof2, cp.foe, cp.m3kf2);
        }
    }
    println!();
    
    println!("Maximum Usable Frequencies:");
    println!("===========================");
    println!("  Basic MUF: {:.2} MHz", path.bmuf);
    println!("  MUF50 (50%): {:.2} MHz", path.muf50);
    println!("  MUF90 (90%): {:.2} MHz", path.muf90);
    println!("  MUF10 (10%): {:.2} MHz", path.muf10);
    println!("  Operational MUF: {:.2} MHz", path.opmuf);
    println!("  OPMUF90: {:.2} MHz", path.opmuf90);
    println!("  OPMUF10: {:.2} MHz", path.opmuf10);
    println!();
    
    println!("Propagation Modes:");
    println!("==================");
    println!("Note: Detailed mode information available in PathData structure");
    println!("F2 layer modes available for analysis");
    println!("E layer modes available for analysis");
    println!();
    
    println!("Signal Analysis:");
    println!("================");
    println!("  Field strength (Es): {:.2} dB(1μV/m)", path.es);
    println!("  Field strength (Ep): {:.2} dB(1μV/m)", path.ep);
    println!("  Received power: {:.2} dBW", path.pr);
    println!("  Signal-to-noise ratio: {:.2} dB", path.snr);
    println!("  Required SNR: {:.1} dB", path.snrr);
    println!("  Required SIR: {:.1} dB", path.sirr);
    println!();
    
    println!("Noise Analysis:");
    println!("===============");
    println!("  Frequency dispersion: {:.2} Hz", path.f0);
    println!("  Time spread: {:.2} ms", path.t0);
    println!();
    
    println!("Reliability Analysis:");
    println!("====================");
    println!("  Basic Circuit Reliability: {:.1}%", path.bcr);
    println!("  Upper decile deviation (SNR): {:.2} dB", path.du_sn);
    println!("  Lower decile deviation (SNR): {:.2} dB", path.dl_sn);
    println!();
    
    println!("Antenna Information:");
    println!("===================");
    println!("  TX Antenna: {}", path.a_tx.name);
    println!("  RX Antenna: {}", path.a_rx.name);
    println!("  Receiving antenna gain: {:.1} dBi", path.grw);
    println!("  EIRP: {:.1} dBW", path.eirp);
    println!();
    
    println!("Technical Parameters:");
    println!("====================");
    println!("  Modulation: {}", if path.modulation == 0 { "Analog" } else { "Digital" });
    println!("  Path type: {}", if path.sor_l == 0 { "Short" } else { "Long" });
    println!("  SNR percentage: {}%", path.snrxxp);
    println!();
    
    println!("Interpretation:");
    println!("===============");
    if path.frequency <= path.muf90 {
        println!("  Frequency suitability: Excellent - Frequency well below MUF90");
    } else if path.frequency <= path.muf50 {
        println!("  Frequency suitability: Good - Frequency is below median MUF");
    } else if path.frequency <= path.muf10 {
        println!("  Frequency suitability: Marginal - Frequency near upper MUF limit");
    } else {
        println!("  Frequency suitability: Poor - Frequency above MUF10");
    }
    
    if path.snr >= path.snrr + 10.0 {
        println!("  Signal quality: Excellent - SNR well above requirement");
    } else if path.snr >= path.snrr {
        println!("  Signal quality: Good - SNR meets requirement");
    } else {
        println!("  Signal quality: Poor - SNR below requirement");
    }
    
    if path.bcr >= 90.0 {
        println!("  Circuit reliability: Excellent - Very high reliability");
    } else if path.bcr >= 70.0 {
        println!("  Circuit reliability: Good - Acceptable reliability");
    } else if path.bcr >= 50.0 {
        println!("  Circuit reliability: Marginal - Moderate reliability");
    } else {
        println!("  Circuit reliability: Poor - Low reliability");
    }
    
    println!();
    println!("Recommendations:");
    println!("================");
    println!("  - Optimal frequency range: {:.1} - {:.1} MHz", 
             path.muf90 * 0.8, path.muf50);
}
