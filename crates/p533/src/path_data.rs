use p372::NoiseParams;
use crate::constants::*;

#[derive(Debug, Clone, Copy)]
pub struct Location {
    pub lat: f64,
    pub lng: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct SolarParameters {
    pub ha: f64,      // hour angle (radians)
    pub sha: f64,     // Sunrise/sunset hour angle (radians)
    pub sza: f64,     // Solar zenith angle (radians)
    pub decl: f64,    // Solar declination (radians)
    pub eot: f64,     // Equation of time (minutes)
    pub lsr: f64,     // local sunrise (hours)
    pub lsn: f64,     // local solar noon (hours)
    pub lss: f64,     // local sunset (hours)
}

impl Default for SolarParameters {
    fn default() -> Self {
        SolarParameters {
            ha: 0.0,
            sha: 0.0,
            sza: 0.0,
            decl: 0.0,
            eot: 0.0,
            lsr: 0.0,
            lsn: 0.0,
            lss: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ControlPt {
    pub l: Location,
    pub distance: f64,    // This is the distance (km) from the transmitter to the CP and not the hop range
    pub foe: f64,         // E layer critical frequency (MHz)
    pub fof2: f64,        // F2 layer critical frequency (MHz)
    pub m3kf2: f64,       // F2 layer critical frequency @ 3000 km (MHz)
    pub dip: [f64; 2],    // Magnetic dip (radians)
    pub fh: [f64; 2],     // Gyrofrequency (MHz)
    pub ltime: f64,       // Local time (hours)
    pub hr: f64,          // Mirror reflection point (km)
    pub x: f64,           // foE/foF2 ratio used in the calculation of the F2MUF
    pub sun: SolarParameters, // Solar parameters
}

impl Default for ControlPt {
    fn default() -> Self {
        ControlPt {
            l: Location { lat: 0.0, lng: 0.0 },
            distance: 0.0,
            foe: 0.0,
            fof2: 0.0,
            m3kf2: 0.0,
            dip: [0.0; 2],
            fh: [0.0; 2],
            ltime: 0.0,
            hr: 0.0,
            x: 0.0,
            sun: SolarParameters::default(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Mode {
    // Define the myriad of MUFs
    pub bmuf: f64,    // Basic MUF (MHz). Typically there is no difference between the basic and the 50% MUF
    // The BMUF is checked to see if it is != 0.0 to determine if the mode exists
    pub muf90: f64,   // MUF exceeded for 90% of the days of the month (MHz)
    pub muf50: f64,   // MUF exceeded for 50% of the days of the month(MHz)
    pub muf10: f64,   // MUF exceeded for 10% of the days of the month(MHz)
    pub opmuf: f64,   // Operation MUF(MHz)
    pub opmuf10: f64, // Operation MUF exceeded 10% of the days of the month(MHz)
    pub opmuf90: f64, // Operation MUF exceeded 90% of the days of the month(MHz)
    pub fprob: f64,   // Probability that the mode is supported at the frequency of interest
    pub deltal: f64,  // Lower decile for the MUF calculations
    pub deltau: f64,  // Upper decile for the MUF calculations
    // Other parameters associated with the mode
    pub hr: f64,      // Reflection height for the mode
    pub fs: f64,      // E-Layer screening frequency for F2 modes only(MHz)
    pub lb: f64,      // < 9000 km path basic loss
    pub ew: f64,      // < 9000 km field strength(dB(1 uV/m))
    pub ele: f64,     // Elevation angle
    pub prw: f64,     // Receiver power (dBW)
    pub grw: f64,     // Receive antenna gain (dBi)
    pub tau: f64,     // Time delay
    pub mc: i32,
}

impl Default for Mode {
    fn default() -> Self {
        Mode {
            bmuf: 0.0,
            muf90: 0.0,
            muf50: 0.0,
            muf10: 0.0,
            opmuf: 0.0,
            opmuf10: 0.0,
            opmuf90: 0.0,
            fprob: 0.0,
            deltal: 0.0,
            deltau: 0.0,
            hr: 0.0,
            fs: 0.0,
            lb: 0.0,
            ew: 0.0,
            ele: 0.0,
            prw: 0.0,
            grw: 0.0,
            tau: 0.0,
            mc: 0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Beam {
    pub azm: f64,     // Azimuth
    pub ele: f64,     // Elevation angle
    pub g: f64,       // Gain for the azimuth and elevation
}

#[derive(Debug, Clone)]
pub struct Antenna {
    pub name: String,

    /*
     * Int used to track the number of frequencies for which we have pattern data
     * (e.g. the size of the freqs array).
     */
    pub freqn: i32,

    /*
     * An array to store the frequencies we have pattern data for.
     */
    pub freqs: Vec<f64>,

    // 3D double pointer to the antenna pattern data
    // [freq_index][azimuth][elevation]
    // The following is assumed about the antenna pattern when the program is run:
    //      i) The orientation is correct. The antenna pattern is in the orientation as it would be on the Earth.
    //      ii) The data is valid. It is the responsibility of the calling program to ensure this.
    pub pattern: Vec<Vec<Vec<f64>>>,
}

impl Default for Antenna {
    fn default() -> Self {
        Antenna {
            name: String::new(),
            freqn: 0,
            freqs: Vec::new(),
            pattern: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct PathData {
    // User-provided Input ************************************************************************

    pub year: i32,
    pub month: i32,          // Note: This is 0 - 11
    pub hour: i32,           // Note: This is an hour index 0 - 23
    //       Where 1 - 24 UTC is required add one and rollover
    pub ssn: i32,            // 12-month smoothed sun sport number a.k.a. R12

    pub modulation: i32,     // Modulation flag

    pub sor_l: i32,          //  Short or long path switch

    pub frequency: f64,      // Frequency (MHz)
    pub bw: f64,             // Bandwidth (Hz)

    pub txpower: f64,        // Transmitter power (dB(1 kW))

    pub snrxxp: i32,         // Required signal-to-noise ration (%) of the time (1 to 99)
    pub snrr: f64,           // Required signal-to-noise ratio (dB)
    pub sirr: f64,           // Required signal-to-interference ratio (dB)

    // Parameters for approximate basic circuit reliability for digital modulation
    pub f0: f64,             // Frequency dispersion at a level -10 dB relative to the peak signal amplitude
    pub t0: f64,             // Time spread at a level -10 dB relative to the peak signal amplitude

    // Parameters for digital modulation performance
    pub a: f64,              // Required A ratio (dB)
    pub tw: f64,             // Time window (msec)
    pub fw: f64,             // Frequency window (Hz)

    pub l_tx: Location,
    pub l_rx: Location,
    pub a_tx: Antenna,
    pub a_rx: Antenna,

    // End User Provided Input *********************************************************************

    // Array pointers ******************************************************************************
    // The advantage of having these pointers in the PathData structure is that p533() can be
    // re-entered with the data allocations intact since they are determined and loaded externally
    // to p533(). This is done to make area coverage calculations, multiple hours and/or
    // any calculations that require the path be examined for another location or time within the
    // current month. If the month changes foF2 and M3kF2 will have to be reloaded, while the pointer
    // foF2var does not since it is for the entire year
    // Pointers to array extracted from the coefficients in ~/IonMap directory
    pub fof2: Vec<Vec<Vec<Vec<f32>>>>,         // foF2
    pub m3kf2: Vec<Vec<Vec<Vec<f32>>>>,        // M(3000)F2
    // Pointer to array extracted from the file "P1239-2 Decile Factors.txt"
    pub fof2var: Vec<Vec<Vec<Vec<Vec<f64>>>>>, // foF2 Variability from ITU-R P.1239-2 TABLE 2 and TABLE 3

    // End Array Pointers *************************************************************************

    // Calculated Parameters **********************************************************************
    pub season: i32,         // This is used for MUF calculations
    pub distance: f64,       // This is the great circle distance (km) between the rx and tx
    pub ptick: f64,          // Slant range
    pub dmax: f64,           // d sub max (km) determined as a function of the midpoint of the path and other parameter
    pub b: f64,              // Intermediate value when calculating dmax also determined at midpoint of the path
    pub ele: f64,            // For paths that are longer than 9000 km this is the composite elevation angle

    // MUFs
    pub bmuf: f64,           // Basic MUF (MHz)
    pub muf50: f64,          // MUF exceeded for 50% of the days of the month (MHz)
    pub muf90: f64,          // MUF exceeded for 90% of the days of the month (MHz)
    pub muf10: f64,          // MUF exceeded for 10% of the days of the month (MHz)
    pub opmuf: f64,          // Operation MUF (MHz)
    pub opmuf90: f64,        // OPMUF exceeded for 90% of the days of the month (MHz)
    pub opmuf10: f64,        // OPMUF exceeded for 10% of the days of the month (MHz)
    // Highest probable frequency, HPF, is 10% MUF (MHz)
    // Optimum working frequency, FOT, is 90% MUF (MHz)

    pub n0_f2: i32,          // Lowest order F2 mode ( 0 to MAXF2MODES )
    pub n0_e: i32,           // Lowest order E mode ( 0 to 2 )

    // Signal powers
    pub es: f64,  // The overall resultant equivalent median sky-wave field strength for path->distance < 7000 km
    pub el: f64,  // The overall resultant median field strength for paths->distance > 9000 km
    pub ei: f64,  // For paths->distance between 7000 and 9000 km the interpolated resultant median field strength
    pub ep: f64,  // The Path Field Strength (dBu) Depending on the path distance this is either Es, El or Ei.
    pub pr: f64,  // Median available receiver power

    // Short path (< 7000 km) parameters
    pub lz: f64,             // "Not otherwise included" loss

    // Long path (> 9000 km) parameters
    pub e0: f64,             // The free-space field strength for 3 MW EIRP
    pub gap: f64,            // Focusing on long distance gain (dB)
    pub ly: f64,             // "Not otherwise included" loss
    pub fm: f64,             // Upper reference frequency
    pub fl: f64,             // Lower reference frequency
    pub f: f64,              // f(f, fH, fL, fM) in eqn 28 P.533-12
    pub fh: f64,             // Mean gyrofrequency
    pub gtl: f64,            // Largest antenna gain in the range 0 to 8 degrees
    pub k: [f64; 2],         // Correction factor

    // Signal-to-noise ratio
    pub snr: f64,            // Median resultant signal-to-noise ratio (dB) for bandwidth b (Hz)
    pub du_sn: f64,          // Upper decile deviation of the signal-to-noise ratio (dB)
    pub dl_sn: f64,          // Lower decile deviation of the signal-to-noise ratio (dB)

    // Signal-to-noise at the required reliability
    pub snrxx: f64,          //

    // Digitally modulated system stats
    pub sir: f64,   // Signal-to-interference ratio (db)
    pub du_si: f64, // Upper decile deviation of the signal-to-interference ratio (db)
    pub dl_si: f64, // Lower decile deviation of the signal-to-interference ratio (db)
    pub rsn: f64,   // Probability that the required SNR is achieved
    pub rt: f64,    // Probability that the required time spread T0 is not exceeded
    pub rf: f64,    // Probability that the required frequency spread f0 is not exceeded

    // Reliability
    pub bcr: f64,     // Basic circuit reliability
    pub ocr: f64,     // Overall circuit reliability without scattering
    pub ocrs: f64,    // Overall circuit reliability with scattering
    pub mir: f64,     // Multimode interference
    pub probocc: f64, // Probability of scattering occurring

    // There are a maximum of 5 CP from P.533-12 Table 1d)
    // See #define above for "Control point index names for readability"
    pub cp: [ControlPt; 5],

    // Antenna related parameters

    // Grw
    //  path->distance <= 7000 km
    //      Grw is the "lossless receiving antenna of gain Grw
    //      (dB relative to an isotropic radiator) in the direction of signal incidence"
    //      Grw will be the dominant mode gain
    //  path->distance >= 9000 km
    //      Grw is the "largest value of receiving antenna gain at the required azimuth in the
    //      elevation range 0 to 8 degrees."
    pub grw: f64,

    // Transmitter EIRP
    pub eirp: f64,

    // ITU-R P.533-12 5.2.1 modes considered "Up to three E modes (for paths up to 4000 km) and
    // up to six F2 modes are selected"
    // In part three of P.533-12 it would have been easier to make all nine modes in one array for digitally
    // modulated systems. To increase the readability and because the method often treats layers differently
    // the modes are separated by layer.
    pub md_f2: [Mode; MAX_F2_MDS],
    pub md_e: [Mode; MAX_E_MDS],

    // The following are conveniences for examining data
    // The variables *DMptr and DMidx are set in MedianAvailableReceiverPower()
    pub dm_ptr: Option<Mode>, // Pointer to the dominant mode
    pub dm_idx: i32,          // Index to the dominant mode (0-2) E layer (3-8) F2 layer

    // Noise Structure
    pub noise_p: NoiseParams,

    // End Calculated Parameters *****************************************************************************
}

impl Default for PathData {
    fn default() -> Self {
        let mut path = PathData {
            year: 0,
            month: 0,
            hour: 0,
            ssn: 0,
            modulation: 0,
            sor_l: 0,
            frequency: 0.0,
            bw: 0.0,
            txpower: 0.0,
            snrxxp: 0,
            snrr: 0.0,
            sirr: 0.0,
            f0: 0.0,
            t0: 0.0,
            a: 0.0,
            tw: 0.0,
            fw: 0.0,
            l_tx: Location { lat: 0.0, lng: 0.0 },
            l_rx: Location { lat: 0.0, lng: 0.0 },
            a_tx: Antenna::default(),
            a_rx: Antenna::default(),
            season: 0,
            distance: 0.0,
            ptick: 0.0,
            dmax: 0.0,
            b: 0.0,
            ele: 0.0,
            bmuf: 0.0,
            muf50: 0.0,
            muf90: 0.0,
            muf10: 0.0,
            opmuf: 0.0,
            opmuf90: 0.0,
            opmuf10: 0.0,
            n0_f2: 0,
            n0_e: 0,
            es: 0.0,
            el: 0.0,
            ei: 0.0,
            ep: 0.0,
            pr: 0.0,
            lz: 0.0,
            e0: 0.0,
            gap: 0.0,
            ly: 0.0,
            fm: 0.0,
            fl: 0.0,
            f: 0.0,
            fh: 0.0,
            gtl: 0.0,
            k: [0.0; 2],
            snr: 0.0,
            du_sn: 0.0,
            dl_sn: 0.0,
            snrxx: 0.0,
            sir: 0.0,
            du_si: 0.0,
            dl_si: 0.0,
            rsn: 0.0,
            rt: 0.0,
            rf: 0.0,
            bcr: 0.0,
            ocr: 0.0,
            ocrs: 0.0,
            mir: 0.0,
            probocc: 0.0,
            cp: std::array::from_fn(|_| ControlPt::default()),
            grw: 0.0,
            eirp: 0.0,
            md_f2: [Mode::default(); MAX_F2_MDS],
            md_e: [Mode::default(); MAX_E_MDS],
            dm_ptr: None,
            dm_idx: 0,
            noise_p: NoiseParams::default(),
            fof2: Vec::new(),
            m3kf2: Vec::new(),
            fof2var: Vec::new(),
        };

        // Initialize ionospheric parameter arrays
        let hrs = 24;   // 24 hours
        let lng = 241;  // 241 longitudes at 1.5 degree increments
        let lat = 121;  // 121 latitudes at 1.5 degree increments
        let ssn = 2;    // 2 SSN (12-month smoothed sun spot numbers) high and low

        path.fof2 = vec![vec![vec![vec![0.0f32; ssn]; lat]; lng]; hrs];
        path.m3kf2 = vec![vec![vec![vec![0.0f32; ssn]; lat]; lng]; hrs];

        // Initialize foF2 variability arrays
        let season = 3; // 3 seasons
        //      1) WINTER 2) EQUINOX 3) SUMMER
        let hrs = 24;   // 24 hours  
        let lat = 19;   // 19 latitude by 5
        //      0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90
        let ssn = 3;    // 3 SSN ranges
        //      1) R12 < 50 2) 50 <= R12 <= 100 3) R12 > 100
        let decile = 2; // 2 deciles 
        //  1) lower 2) upper

        path.fof2var = vec![vec![vec![vec![vec![0.0f64; decile]; ssn]; lat]; hrs]; season];

        path
    }
}
