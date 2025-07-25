/// Constants for P533 calculations

pub const TRUE: i32 = 1;
pub const FALSE: i32 = 0;
pub const PI: f64 = 3.14159265358979323846;
pub const R0: f64 = 6371.009; // km International Union of Geodesy and Geophysics mean Earth radius
pub const D2R: f64 = 0.0174532925; // PI/180
pub const R2D: f64 = 57.2957795; // 180/PI
pub const VOF_L: f64 = 299792458.0; // Velocity of light (m/s)

pub const TINY_DB: i32 = -307;
pub const TOO_BIG: f64 = 1.7976931348623157E+308;

pub const DBL_MAX: f64 = 1.7976931348623157E+308;
pub const DBL_MIN: f64 = 2.2250738585072014E-308;
pub const DBL_EPSILON: f64 = 2.2204460492503131E-016;

// Control point index names for readability
// These are defined from the sense of the short model
// Please note these change meaning when the long model is exclusively used
// i.e when the path->distance is > 9000 km. This is done for diagnostic purposes.
pub const T1K: usize = 0; // T + 1000 (km)
                          // Note: Alternative use in long model penetration point closest to the transmitter at the current hour
pub const TD02: usize = 1; // T + d0/2 (km)
pub const MP: usize = 2; // path mid-path (km);
pub const RD02: usize = 3; // R - d0/2 (km)
                           // Note: Alternative use in lone model as R - dM/2
pub const R1K: usize = 4; // R - 1000 (km)
                          // Note: Alternative use in long model at last penetration point, 2*nL, at current hour

// foF2 variability index names for readability
pub const WINTER: usize = 0;
pub const EQUINOX: usize = 1;
pub const SUMMER: usize = 2;

// Decile flags
pub const DL: usize = 0; // Lower decile
pub const DU: usize = 1; // Upper decile

pub const DAY: usize = 0; // DAY index for the rop array in CalculateMUFOperational()
pub const NIGHT: usize = 1; // NIGHT index for the rop array in CalculateMUFOperational()

pub const JAN: i32 = 0;
pub const FEB: i32 = 1;
pub const MAR: i32 = 2;
pub const APR: i32 = 3;
pub const MAY: i32 = 4;
pub const JUN: i32 = 5;
pub const JUL: i32 = 6;
pub const AUG: i32 = 7;
pub const SEP: i32 = 8;
pub const OCT: i32 = 9;
pub const NOV: i32 = 10;
pub const DEC: i32 = 11;

// For the determination of the lowest order E and F2 mode
pub const NO_LOWEST_MODE: i32 = 99;

// Modulation flags
pub const ANALOG: i32 = 0;
pub const DIGITAL: i32 = 1;

// Maximum Sun Spot Number
pub const MAX_SSN: i32 = 160;

// Maximum number of F2 modes
pub const MAX_F2_MDS: usize = 6;

// Maximum number of E modes
pub const MAX_E_MDS: usize = 3;

// Maximum number of modes
pub const MAX_MDS: usize = MAX_E_MDS + MAX_F2_MDS;

// Direction of the AntennaGain()
pub const TX_TO_RX: i32 = 1;
pub const RX_TO_TX: i32 = 2;
pub const TXTORX: i32 = 1;
pub const RXTOTX: i32 = 2;

// Short or long path flags
pub const SHORT_PATH: i32 = 0;
pub const LONG_PATH: i32 = 1;

// Minimum Elevation Angle (degrees) for the Short model
pub const MIN_ELE_ANGLE_S: f64 = 3.0;
// Minimum Elevation Angle (degree) for the Long model
pub const MIN_ELE_ANGLE_L: f64 = 3.0;

// Indices for magfit() The gyrofrequency and magnetic dip are calculated at a height
//      of either 100 or 300 km. For absorption calculation 100 km is typically used
//      while 300 km is used for other calculations.
pub const HR_100_KM: usize = 0; // height = 100 (km)
pub const HR_300_KM: usize = 1; // height = 300 (km)

// Note: Alternative use in long model as T + dM/2
