pub mod constants;
pub mod path_data;
pub mod geometry;
pub mod initialize_path_util;
pub mod muf_basic;
pub mod muf_variability;
pub mod muf_operational;
pub mod e_layer_screening_frequency;
pub mod median_skywave_field_strength_short_util;
pub mod median_skywave_field_strength_long_util;
pub mod between_7000km_and_9000km_util;
pub mod median_available_receiver_power_util;
pub mod circuit_reliability_util;
pub mod read_ion_parameters;
pub mod read_p1239_util;
pub mod calculate_cp_parameters;

pub use constants::*;
pub use path_data::*;

/// This program provides methods for the prediction of available frequencies, signal levels, and the predicted reliability
/// for analogue and digital-modulated HF systems, taking into account not only the signal-to-noise ratio but also the expected
/// time and frequency spreads of the channel. This program calculates the HF path parameters that appear in ITU-R P.533-12.
/// This program is loosely based on the program REC533. This model uses control points to determine the modes of propagation of
/// HF signals in the ionosphere.
pub fn calculate_p533(path: &mut PathData) {
    // Calculate the distances between rx and tx, find the midpoint of the path, find the midpoint distance and initialize the path 
    // This will also determine the ionospheric parameters for 3 of the potential 5 control points.
    initialize_path_util::initialize_path(path);

    // Load static data files after path initialization
    let _ = read_ion_parameters::initialize_ionospheric_data(path);
    let _ = read_p1239_util::read_p1239(path);

    // Calculate control point parameters for all control points
    for i in 0..5 {
        let mut cp = path.cp[i].clone();
        calculate_cp_parameters::calculate_cp_parameters(path, &mut cp);
        path.cp[i] = cp;
    }

    /************************************************************/
    /* Part 1 – Frequency availability                          */
    /************************************************************/
    /*
       The following routines in Part 1:
            i) MUFBasic()
            ii) MUFVariability()
            iii) MUFOperational()
       Have to be executed in the order above because they slowly populate the path structure
       as the calculation progresses. Typically they follow the flow of ITU-R P533-11
       since the standard is not necessarily a description of a computer algorithm some of the 
       calculation is out of sequence with ITU-R P.533-12.

       Note: The following 4 subroutines are only applicable to paths less than or equal to 9000 km
     */

    // Determine the basic MUF (BMUF) This will also determine R - d0/2 and T - d0/2
    // Control points if necessary.
    muf_basic::calculate_muf_basic(path);

    // Determine Fprob for each mode and the path the 50% MUF (MUF50), 90% MUF (MUF90) and the 10% MUF (MUF10)
    muf_variability::calculate_muf_variability(path);

    // Determine the for each mode and the path the operational MUF (OPMUF), 90% OPMUF (OPMUF90) and the 10% OPOMUF (OPMUF10)
    muf_operational::calculate_muf_operational(path);

    // E Layer Screening Frequency is determine contingent on the path length
    e_layer_screening_frequency::calculate_e_layer_screening_frequency(path);

    /************************************************************/
    /* Part 2 – Median sky-wave field strength                  */
    /************************************************************/
    /*
     * Each of the routines below will initially check the path->distance to determine if the calculation should proceed.
     * As in Part 1 above these routines are designed to be executed in the following order 
     * 		i)		MedianSkywaveFieldStrengthShort()	Calculation for path->distance < 9000 km
     *		ii)		MedianSkywaveFieldStrengthLong()	Calculation for path->distance > 9000 km
     *		iii)	Between7000kmand9000km()			Interpolation for path->distance between 7000 and 9000 km
     */

    median_skywave_field_strength_short_util::median_skywave_field_strength_short(path);

    median_skywave_field_strength_long_util::median_skywave_field_strength_long(path);

    between_7000km_and_9000km_util::between_7000km_and_9000km(path);

    median_available_receiver_power_util::median_available_receiver_power(path);

    /************************************************************/
    /* Part 3 – The prediction of system performance            */
    /************************************************************/

    // Call noise from the P372 crate 
    // Note: Using a default terrain category for now - this should be configurable
    path.noise_p = p372::noise(path.hour as u32, path.l_rx.lng, path.l_rx.lat, p372::TerrainCategory::Rural, path.frequency);

    circuit_reliability_util::circuit_reliability(path);
}

/// Create and calculate a new P533 path prediction
/// This is a high-level convenience function for creating path predictions
pub fn create_path_prediction(
    tx_lat: f64, tx_lng: f64, rx_lat: f64, rx_lng: f64,
    year: i32, month: i32, hour: i32, ssn: i32,
    frequency: f64, txpower: f64,
    tx_antenna_gain: f64, rx_antenna_gain: f64
) -> PathData {
    use crate::path_data::{PathData, Location, Antenna};

    // Initialize path data structure
    let mut path = PathData {
        year,
        month,
        hour,
        ssn,
        modulation: ANALOG,
        sor_l: SHORT_PATH,
        frequency,
        bw: 3000.0, // Default 3 kHz bandwidth
        txpower,
        snrxxp: 90,
        snrr: 20.0,
        sirr: 20.0,
        f0: 0.0,
        t0: 0.0,
        a: 0.0,
        tw: 0.0,
        fw: 0.0,
        l_tx: Location { lat: tx_lat * std::f64::consts::PI / 180.0, lng: tx_lng * std::f64::consts::PI / 180.0 },
        l_rx: Location { lat: rx_lat * std::f64::consts::PI / 180.0, lng: rx_lng * std::f64::consts::PI / 180.0 },
        a_tx: Antenna {
            name: "TX Antenna".to_string(),
            freqn: 1,
            freqs: vec![frequency],
            pattern: vec![vec![vec![tx_antenna_gain; 91]; 360]; 1],
        },
        a_rx: Antenna {
            name: "RX Antenna".to_string(),
            freqn: 1,
            freqs: vec![frequency],
            pattern: vec![vec![vec![rx_antenna_gain; 91]; 360]; 1],
        },
        ..Default::default()
    };

    // Run the full calculation
    calculate_p533(&mut path);

    path
}

/// Calculate circuit reliability percentage
/// This function provides a high-level interface for circuit reliability calculations
pub fn calculate_circuit_reliability(
    tx_lat: f64, tx_lng: f64, rx_lat: f64, rx_lng: f64,
    year: i32, month: i32, hour: i32, ssn: i32,
    frequency: f64, required_snr: f64,
    tx_power: f64, tx_antenna_gain: f64, rx_antenna_gain: f64
) -> f64 {
    let mut path = create_path_prediction(
        tx_lat, tx_lng, rx_lat, rx_lng,
        year, month, hour, ssn,
        frequency, tx_power,
        tx_antenna_gain, rx_antenna_gain
    );
    
    path.snrr = required_snr;
    
    // Return basic estimate based on frequency vs MUF
    if path.muf90 > 0.0 {
        if frequency <= path.muf90 {
            90.0 // 90% reliability if frequency is below 90% MUF
        } else if frequency <= path.muf50 {
            50.0 // 50% reliability if frequency is below median MUF
        } else if frequency <= path.muf10 {
            10.0 // 10% reliability if frequency is below 10% MUF
        } else {
            0.0 // Very low reliability above MUF
        }
    } else {
        0.0
    }
}

/// Get MUF information for a path
/// This provides access to the Maximum Usable Frequency calculations
pub fn get_muf_info(
    tx_lat: f64, tx_lng: f64, rx_lat: f64, rx_lng: f64,
    year: i32, month: i32, hour: i32, ssn: i32
) -> (f64, f64, f64, f64, f64) {
    let path = create_path_prediction(
        tx_lat, tx_lng, rx_lat, rx_lng,
        year, month, hour, ssn,
        10.0, 100.0, // Default frequency and power
        0.0, 0.0     // Default antenna gains
    );
    
    (path.muf50, path.muf90, path.muf10, path.opmuf, path.distance)
}
