use crate::constants::ControlPointIndexes;
use crate::geometry::great_circle_point;
use std::f64::consts::PI;

mod constants;
mod geometry;

#[derive(Clone)]
pub struct Location {
    pub latitude: f64,
    pub longitude: f64,
}

pub struct Antenna {
    pub freqn: u32,
    pub frequencies: Vec<f64>,
    pub pattern: Vec<Vec<f64>>,
}

pub struct P533Input {
    pub year: u32,
    pub month: u32,
    pub hour: u32,

    pub ssn: u32,

    pub modulation: u32,

    pub short_or_long_path_switch: ShortOrLongPathSwitch,

    // Mhz
    pub frequency: f64,
    // Hz
    pub bandwidth: f64,

    pub tx_power: f64,

    pub snr_xxp: u32,
    pub snr_r: f64,
    pub sir_r: f64,

    pub f0: f64,
    pub t0: f64,

    pub a: f64,
    pub tw: f64,
    pub fw: f64,

    pub tx_location: Location,
    pub rx_location: Location,

    pub tx_antenna: Antenna,
    pub rx_antenna: Antenna,
}

#[derive(Clone)]
struct SolarParameters {
    ha: f64,   // hour angle (radians)
    sha: f64,  // Sunrise/sunset hour angle (radians)
    sza: f64,  // Solar zenith angle (radians)
    decl: f64, // Solar declination (radians)
    eot: f64,  // Equation of time (minutes)
    lsr: f64,  // local sunrise (hours)
    lsn: f64,  // local solar noon (hours)
    lss: f64,  // local sunset (hours)
}

#[derive(Clone)]
struct ControlPt {
    l: Location,
    distance: f64, // This is the distance (km) from the transmitter to the CP and not the hop range
    fo_e: f64,     // E layer critical frequency (MHz)
    fo_f2: f64,    // F2 layer critical frequency (MHz)
    m3k_f2: f64,   // F2 layer critical frequency @ 3000 km (MHz)
    dip: Vec<f64>, // Magnetic dip (radians)
    fh: Vec<f64>,  // Gyrofrequency (MHz)
    local_time: f64, // Local time (hours)
    hr: f64,       // Mirror reflection point (km)
    x: f64,        // foE/foF2 ratio used in the calculation of the F2MUF

    sun: SolarParameters, // Solar parameters
}

impl Default for ControlPt {
    fn default() -> Self {
        ControlPt {
            l: Location {
                latitude: 0.0,
                longitude: 0.0,
            },
            distance: 0.0,
            fo_e: 0.0,
            fo_f2: 0.0,
            m3k_f2: 0.0,
            dip: vec![0.0; 2],
            fh: vec![0.0; 2],
            local_time: 0.0,
            hr: 0.0,
            x: 0.0,
            sun: SolarParameters {
                ha: 0.0,
                sha: 0.0,
                sza: 0.0,
                decl: 0.0,
                eot: 0.0,
                lsr: 0.0,
                lsn: 0.0,
                lss: 0.0,
            },
        }
    }
}

#[derive(PartialEq)]
pub enum ShortOrLongPathSwitch {
    Short,
    Long,
}

fn initialise_control_points(
    tx_location: &Location,
    rx_location: &Location,
    distance: f64,
) -> Vec<ControlPt> {
    let mut control_points = vec![ControlPt::default(); 5];

    let fracd = 0.5;
    let mddl = great_circle_point(tx_location, rx_location, distance, fracd);
    control_points[ControlPointIndexes::MP as usize].distance = mddl.distance;
    control_points[ControlPointIndexes::MP as usize].l = mddl.location;

    control_points
}

pub fn calculate_p533(input: P533Input) {
    let mut distance = match geometry::great_circle_distance(&input.tx_location, &input.rx_location)
    {
        val if val == 0.0 => constants::DBL_EPSILON,
        distance => distance,
    };

    if input.short_or_long_path_switch == ShortOrLongPathSwitch::Long {
        distance = constants::EARTH_RADIUS_0 * PI * 2.0 - distance;
    }

    initialise_control_points(&input.tx_location, &input.rx_location, distance);
}
