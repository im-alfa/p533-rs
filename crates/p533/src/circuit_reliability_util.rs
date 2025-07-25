use crate::constants::*;
use crate::geometry;
use crate::path_data::{Location, PathData};

/// Digital modulation signal and interferers calculation
fn digital_modulation_signal_and_interferers(
    path: &mut PathData,
    i_s: &mut [i32; MAX_MDS],
    i_i: &mut [i32; MAX_MDS],
) -> f64 {
    /*
     DigitalModulationSignalandInterferers() Returns the signal that satisfies the amplitude ratio, A,
         and the time window, Tw. It also returns the indicies of the interfering modes in the array, iI[9]
         The latter is used in the calculation of overall circuit reliability, OCR

         INPUT
             struct PathData *path

         OUTPUT
             int iS[MAXMDS] - Index array of the modes that meet the signal criteria
             int iI[MAXMDS] - Index array of the modes that meet the interference criteria

         SUBROUTINES
            ElevationAngle()
            NumberofModes()
            ModeSort()
    */

    // This array is so that all modes can be examined together independant of E or F2 layer
    let mut m = vec![crate::path_data::Mode::default(); MAX_MDS];

    // The following 2 arrays, iPrw and itau, are modes indicies arrays for the digital BCR, SIR and OCR calculation
    let mut i_ew = [NOTINDEX; 9];
    let mut itau = [NOTINDEX; 9];

    // Initialize the order arrays
    for n in 0..MAX_MDS {
        i_ew[n % 9] = NOTINDEX;
        itau[n % 9] = NOTINDEX;
        i_i[n] = NOTINDEX;
        i_s[n] = NOTINDEX;
    }

    //*****************************************************************************************
    // Although the slant range. ptick, was calculated in MedianSkywaveFieldStrengthShort() is was
    // calculated under the E layer screening condition which is not relevant here so ptick must be
    // calculated here.
    if path.distance <= 9000.0 {
        // Determine the time delay for all modes that exist
        // E layer modes
        // E layer reflection height
        let hr = 110.0; // Reflection height

        for n in path.n0_e as usize..MAX_E_MDS {
            if path.md_e[n].bmuf != 0.0 {
                // Mode exists
                let dh = path.distance / (n as f64 + 1.0); // Hop distance
                let delta = crate::e_layer_screening_frequency::elevation_angle(dh, hr);
                let psi = dh / (2.0 * R0);
                let ptick = 2.0 * R0 * (f64::sin(psi) / f64::cos(delta - psi));
                path.md_e[n].tau = (n as f64 + 1.0) * (ptick / VOF_L) * 1000.0;
            }
        }

        // F2 layer modes
        // The reflection height for each F2 mode was calculated in ELayerScreeningFrequency()
        for n in path.n0_f2 as usize..MAX_F2_MDS {
            if path.md_f2[n].bmuf != 0.0 {
                // Mode exists
                let hr = path.md_f2[n].hr;
                let dh = path.distance / (n as f64 + 1.0); // Hop distance
                let delta = crate::e_layer_screening_frequency::elevation_angle(dh, hr);
                let psi = dh / (2.0 * R0);
                let ptick = 2.0 * R0 * (f64::sin(psi) / f64::cos(delta - psi));
                path.md_f2[n].tau = (n as f64 + 1.0) * (ptick / VOF_L) * 1000.0;
            }
        }

        // Do the following if there are 2 or more modes
        if number_of_modes(path) >= 2 {
            // For this calculation the layers don't matter so set up an array of all the modes
            // so that a single loop can be used
            // Point the M[] array at all of the modes in path
            for n in 0..MAX_E_MDS {
                m[n] = path.md_e[n].clone();
            }

            for n in 0..MAX_F2_MDS {
                m[n + 3] = path.md_f2[n].clone();
            }

            // P.533-12 Section 10.2.3 Reliability prediction procedure
            // Step 1: Determination of the dominant mode, Ew
            mode_sort(&mut m, &mut i_ew, DOMINANT); // iEw[0] is the dominant mode

            // Order the modes by time also
            mode_sort(&mut m, &mut itau, SOONEST); // itau[0] is the earlest mode

            // Step 2: All other active modes with strengths exceeding (Ew � A (dB)) are identified.
            // Step 3: The first arriving mode is identified, and all modes within the time window, Tw,
            // measured from the first arriving mode, are identified.
            // Step 4: For path lengths up to 7000 km, a power summation of the modes arriving within the
            // window is made, or for path lengths between 7000 and 9000 km the interpolation procedure is used.
            // the basic circuit reliability, BCR, is determined using the same method as in the analog modulation.
            // Step 2 gives an amplitude criteria and Step 3 gives a delay criteria and determine the modes which fulfill
            // both of these criteria. Step 4 details what to do with the modes that meet the criteria.
            // The zeroth element in the itau[] array is the soonest arriving mode.
            let deltat = m[itau[0] as usize].tau + path.tw / 1000.0; // The tau for each mode is in seconds where the time window, TW, is in mS
                                                                     // Time window criteria
                                                                     // The zeroth element in the iEw[] array is the dominant mode.
            let deltaa = m[i_ew[0] as usize].prw - path.a; // A ratio criteria

            // Initialize signal sum
            let mut ssum = 0.0; // Sum of the field strengths
            for n in 0..MAX_MDS {
                // Sum the mode as signals that satisfy the following criteria:
                //		i) The mode exists
                //		ii) The mode within the A ratio of the dominant mode median received power
                //		iii) The mode arrives within TW of the earlest arriving mode
                if m[n].bmuf != 0.0 {
                    if m[n].prw >= deltaa && m[n].tau <= deltat {
                        ssum += f64::powi(f64::powf(10.0, m[n].ew / 10.0), 2);
                        i_s[n] = n as i32;
                    } else {
                        // If the mode is not determined to be a signal then it is interference.
                        // Store the index of the interfering modes for later use in the calculation of
                        // the signal-to-interference ratio (See Table 3 P.842-4).
                        i_i[n] = n as i32;
                    }
                }
            }

            // Only determine Etw if there is a mode that satisfies the criteria above
            // otherwise Etw is set to the something small
            let etw = if ssum > 0.0 {
                10.0 * f64::log10(f64::sqrt(ssum))
            } else {
                TINY_DB as f64
            };

            // Use the combined mode power, Ssum, and the antenna gain, Grw
            let s = etw + path.grw - 20.0 * f64::log10(path.frequency) - 107.2;

            // The iI[] and iS[] arrays need to be ordered.
            for n in 0..MAX_MDS {
                for m in 0..MAX_MDS {
                    if i_i[n] < i_i[m] && i_i[n] != NOTINDEX && i_i[m] != NOTINDEX {
                        let j = i_i[n];
                        i_i[n] = i_i[m];
                        i_i[m] = j;
                    }

                    if i_s[n] < i_s[m] && i_s[n] != NOTINDEX && i_s[m] != NOTINDEX {
                        let j = i_s[n];
                        i_s[n] = i_s[m];
                        i_s[m] = j;
                    }
                }
            }

            s
        } else {
            // (NumberofModes(*path) < 2)
            // There is only one mode.
            path.pr
        }
    } else {
        // (path.distance > 9000.0)
        path.pr
    }
}

fn equatorial_scattering(path: &mut PathData, i_s: &[i32; MAX_MDS]) {
    /*
     EquatorialScattering() - Determines the equatorial scattering by the method described
            in P.533-12 Section 10.3 "Equatorial Scattering"

            INPUT
                struct PathData *path
                int iS[MAXMDS] - Index array of the modes identified as interference

            OUTPUT
                path.OCRs - Overall circuit reliability with scattering

           SUBROUTINES
               FindFlambdad()
               FindFTl()
    */

    let mut ptspread = [DBL_MIN; MAX_MDS]; // Time spread of the mode which satisfies the amplitude and time criteria
    const TSPREAD: f64 = 1.0; // Standard deviation of the time spread, taken as 1 mS
    let mut pfspread = [0.0f64; 2]; // Frequency spreading of the dominant mode
    const FSPREAD: f64 = 3.0; // Standard deviation of the frequency spread, taken as 3 Hz
    let mut probocc = [0.0f64; MAX_MDS]; // Probability of the occurrence of scattering

    // P.533-12 Section 10.3 "Equatorial Scattering"
    // Step 7: Find the time spread for the modes that are within the amplitude ratio, A, and time
    // window, Tw. The index for these modes can be found in the array iS[9] which was determined
    // by the routine DigitalModulationSignalandInterferers().

    if path.distance <= 9000.0 && i_s[0] != NOTINDEX {
        let tau = path.tw; // Time delay being considered

        for n in 0..MAX_MDS {
            // Initialize PTspread to something tiny, DBL_MIN
            ptspread[n] = DBL_MIN;
            if i_s[n] != NOTINDEX {
                let (pm, taum) = if i_s[n] < MAX_E_MDS as i32 {
                    // The mode is an E mode
                    (
                        path.md_e[i_s[n] as usize].ew,
                        path.md_e[i_s[n] as usize].tau,
                    )
                } else {
                    // The mode is a F2 mode
                    (
                        path.md_f2[i_s[n] as usize].ew,
                        path.md_f2[i_s[n] as usize].tau,
                    )
                };

                ptspread[n] = 0.056
                    * pm
                    * f64::exp(-f64::powi(tau - taum, 2) / (2.0 * f64::powi(TSPREAD, 2)));
            }
        }

        // Step 8: Determine the frequency spreading of the dominant mode
        // Find the dominant mode
        let pm = if i_s[0] < MAX_E_MDS as i32 {
            // The mode is an E mode
            path.md_e[i_s[0] as usize].ew
        } else {
            // The mode is a F2 mode
            path.md_f2[i_s[0] as usize].ew
        };

        // Center frequency
        let fm = path.frequency * 1e6; // (Hz)
                                       // Transmitted center frequency
                                       // Frequency window
        let f = path.fw; // (Hz)
                         // Frequency being considered
        pfspread[0] = 0.056 * pm * f64::exp(-f64::powi(f - fm, 2) / (2.0 * f64::powi(FSPREAD, 2)));
        pfspread[1] = 0.056 * pm * f64::exp(-f64::powi(-f - fm, 2) / (2.0 * f64::powi(FSPREAD, 2)));

        // Step 9: Determine the probability of scattering occuring
        let mut use_cp = FALSE; // Flag
                                // Determine if the time spread component is within the amplitude ratio, A, of the
                                // dominant mode power level
        for n in 0..MAX_MDS {
            if ptspread[n] != DBL_MIN {
                // Note: At this point pm is the dominant mode
                if pm - ptspread[n] >= path.a {
                    use_cp = TRUE;
                }
            }
        }

        // Determine is the frequency spread components are within the amplitude ratio, A, of the
        // dominant mode power level
        for n in 0..2 {
            if pm - pfspread[n] >= path.a {
                use_cp = TRUE;
            }
        }

        // Determine the coefficients for the probocc calculation that are independant of the control point
        let fr = 0.1 + 0.008 * f64::max(path.ssn as f64, 160.0); // Temp
        let fs = 0.55 + 0.45 * f64::sin(60.0 * D2R * (path.month as f64 + 1.0 - 1.5)); // Temp

        if use_cp == TRUE {
            // Use the control points
            for n in 0..MAX_MDS {
                // Examine all modes

                // Initialize variables for each mode
                probocc[n] = 0.0;
                let mut flambdad = 0.0; // Temp
                let mut ftl = 0.0; // Temp

                if i_s[n] != NOTINDEX && i_s[n] >= MAX_E_MDS as i32 {
                    // Does the mode exist and is it an F2 layer mode?
                    if i_s[n] == path.n0_f2 {
                        // Lowest order F2 mode
                        if path.distance <= path.dmax {
                            flambdad = find_flambdad(&path.cp[MP]);
                            ftl = find_ftl(&path.cp[MP]);
                        } else {
                            if ptspread[TD02] >= ptspread[RD02] {
                                flambdad = find_flambdad(&path.cp[TD02]);
                                ftl = find_ftl(&path.cp[TD02]);
                            } else {
                                flambdad = find_flambdad(&path.cp[RD02]);
                                ftl = find_ftl(&path.cp[RD02]);
                            }
                        }
                    } else {
                        // Higher order F2 modes
                        if path.distance <= path.dmax {
                            // Find the largest time scattering by brute force
                            if ptspread[T1K] >= ptspread[R1K] && ptspread[T1K] >= ptspread[MP] {
                                flambdad = find_flambdad(&path.cp[T1K]);
                                ftl = find_ftl(&path.cp[T1K]);
                            } else if ptspread[R1K] >= ptspread[T1K]
                                && ptspread[R1K] >= ptspread[MP]
                            {
                                flambdad = find_flambdad(&path.cp[R1K]);
                                ftl = find_ftl(&path.cp[R1K]);
                            } else if ptspread[MP] >= ptspread[R1K] && ptspread[MP] >= ptspread[T1K]
                            {
                                flambdad = find_flambdad(&path.cp[MP]);
                                ftl = find_ftl(&path.cp[MP]);
                            }
                        } else {
                            if ptspread[T1K] >= ptspread[R1K]
                                && ptspread[T1K] >= ptspread[TD02]
                                && ptspread[T1K] >= ptspread[MP]
                                && ptspread[T1K] >= ptspread[RD02]
                            {
                                flambdad = find_flambdad(&path.cp[T1K]);
                                ftl = find_ftl(&path.cp[T1K]);
                            } else if ptspread[R1K] >= ptspread[T1K]
                                && ptspread[R1K] >= ptspread[MP]
                                && ptspread[R1K] >= ptspread[TD02]
                                && ptspread[R1K] >= ptspread[RD02]
                            {
                                flambdad = find_flambdad(&path.cp[R1K]);
                                ftl = find_ftl(&path.cp[R1K]);
                            } else if ptspread[MP] >= ptspread[R1K]
                                && ptspread[MP] >= ptspread[T1K]
                                && ptspread[MP] >= ptspread[TD02]
                                && ptspread[MP] >= ptspread[RD02]
                            {
                                flambdad = find_flambdad(&path.cp[MP]);
                                ftl = find_ftl(&path.cp[MP]);
                            } else if ptspread[TD02] >= ptspread[RD02]
                                && ptspread[TD02] >= ptspread[MP]
                                && ptspread[TD02] >= ptspread[T1K]
                                && ptspread[TD02] >= ptspread[R1K]
                            {
                                flambdad = find_flambdad(&path.cp[TD02]);
                                ftl = find_ftl(&path.cp[TD02]);
                            } else if ptspread[RD02] >= ptspread[T1K]
                                && ptspread[RD02] >= ptspread[MP]
                                && ptspread[RD02] >= ptspread[R1K]
                                && ptspread[RD02] >= ptspread[TD02]
                            {
                                flambdad = find_flambdad(&path.cp[RD02]);
                                ftl = find_ftl(&path.cp[RD02]);
                            }
                        }
                    }
                }

                probocc[n] = flambdad * ftl * fr * fs;
            }

            // Find the biggest probocc
            let mut biggest = 0.0; // Temp
            for n in 0..MAX_MDS {
                if probocc[n] > biggest {
                    biggest = probocc[n];
                }
            }

            path.ocrs = path.ocr * (1.0 - biggest);
        } else {
            // No equatorial scattering
            path.ocrs = path.ocr;
        }
    } else {
        // For paths longer than 9000 km or no modes, no equatorial scattering
        path.ocrs = path.ocr;
    }
}

// Local defines
// Flags for digital reliability calculation
const DOMINANT: i32 = 0;
const SOONEST: i32 = 1;
const NOTINDEX: i32 = 99;

/// Calculate circuit reliability
/// Based on ITU-R P.533-12 and ITU-R P.842-4
pub fn circuit_reliability(path: &mut PathData) {
    /*
      CircuitReliability() finds the basic circuit reliability, BCR, and the overall circuit reliability, OCR,
            for analog and digital systems.

            The BCR is calculated as described in Table 1 ITU-R P.842-4.
            The calculation is broken down into 11 steps. Many of the parameters necessary for the computation have already been
            calculated elsewhere in this project. The median available receiver power of the wanted signal, which is
            Step 1 in Table 1 P.842-4, is calculated in the subroutine MedianAvailableReceiverPower(). The subroutine
            MedianAvailableReceiverPower() implements the calculation in P.533-12 Section 6 Median available receiver.
            Table 1 P.842-4 Steps 2, 5 and 8 are calculated in subroutine Noise() and summarized below.

                     P.842-4 Table 1 Step
                                2			Median noise factors for atmospheric, galactic and man-made noise
                                5			Lower decile deviation of atmospheric, galactic and Man-made noise
                                8			Upper decile deviation of atmospheric, galactic and Man-made noise

            This subroutine, CircuitReliability(), completes the remaining steps in P.842-4 Table 1.

            The simplified approximate BCR method from section 9 "BCR for digital modulation systems" P.842-4 is
            determined from the three probabilities:
                    i) Probability that the required signal-to-noise ratio, SN0, is achieved
                    ii) Probability that the required time spread, T0, at a level of –10 dB relative to the
                        peak signal amplitude is not exceeded
                    iii) Probability that the required frequency dispersion f0 at a level of –10 dB relative to the
                        peak signal amplitude is not exceeded

            This subroutine also determines the overall circuit reliability, OCR, via the method given in ITU-R P.P.842-4 Table 3
            for digital systems.

            This subroutine also determines the equatorial scattering occurrence probability and then calculates the OCR in the
            presence of scattering
    */

    // NORM is an array of independent variables that give
    // the normal cdf from 0.5 to 0.99 in 0.01 increments
    // with a standard deviataion of one and mean zero
    let norm: [f64; 50] = [
        0.0000000000,
        0.0250689082,
        0.0501535835,
        0.075269862,
        0.1004337206,
        0.1256613469,
        0.1509692155,
        0.1763741647,
        0.201893479,
        0.2275449764,
        0.2533471029,
        0.2793190341,
        0.3054807878,
        0.331853346,
        0.3584587930,
        0.3853204663,
        0.4124631294,
        0.4399131658,
        0.467698799,
        0.4958503478,
        0.5244005133,
        0.5533847202,
        0.5828415079,
        0.612812991,
        0.6433454057,
        0.6744897502,
        0.7063025626,
        0.7388468486,
        0.772193213,
        0.8064212461,
        0.8416212327,
        0.8778962945,
        0.9153650877,
        0.954165253,
        0.9944578841,
        1.036433391,
        1.080319342,
        1.12639113,
        1.174986792,
        1.226528119,
        1.281551564,
        1.340755033,
        1.405071561,
        1.47579103,
        1.554773595,
        1.644853625,
        1.750686073,
        1.880793606,
        2.053748909,
        2.326347874,
    ];

    let mut i_i: [i32; MAX_MDS] = [NOTINDEX; MAX_MDS]; // Interference mode indices
    let mut i_s: [i32; MAX_MDS] = [NOTINDEX; MAX_MDS]; // Signal mode indices

    // Table 2 P.842-4
    let table2_ld: [[f64; 10]; 2] = [
        [8.0, 12.0, 13.0, 10.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0],
        [11.0, 16.0, 17.0, 13.0, 11.0, 11.0, 11.0, 9.0, 8.0, 7.0],
    ];
    let table2_ud: [[f64; 10]; 2] = [
        [6.0, 8.0, 12.0, 13.0, 12.0, 9.0, 9.0, 8.0, 7.0, 7.0],
        [9.0, 11.0, 12.0, 13.0, 12.0, 9.0, 9.0, 8.0, 7.0, 7.0],
    ];

    let mut geomag = Location { lat: 0.0, lng: 0.0 };

    // For readability - access noise parameters
    let noise_p = path.noise_p;

    // Begin BCR Calculation *********************************************************************

    // Step 1: "Median available receiver power of wanted signal (dBW)"
    //		See MedianAvailableReceiverPower()

    // Step 2: Median noise factor for atmospheric noise, galactic and man-made noise
    //		See Noise() [P372.dll]

    // Step 3: "Median resultant signal-to-noise ratio (dB) for bandwidth b (Hz)"
    // There is a difference between digital and analog modulation in how the signal is found
    // Once the signal is found for each modulation the BCR calculation is identical
    let s = if path.modulation == ANALOG {
        path.pr
    } else {
        // The signal for digital modulation requires the calculation given in
        // ITU-R P.533-12 Section 10.2.3 Reliability prediction procedure.
        digital_modulation_signal_and_interferers(path, &mut i_s, &mut i_i)
    };

    let atmospheric_noise = noise_p.get_atmospheric_noise();
    let man_made_noise = noise_p.get_man_made_noise();
    let galactic_noise = noise_p.get_galactic_noise();

    path.snr = s
        - 10.0
            * (
                10.0_f64.powf(atmospheric_noise / 10.0) + // Atmospheric noise (FaA)
        10.0_f64.powf(man_made_noise / 10.0) +    // Man-made noise (FaM) 
        10.0_f64.powf(galactic_noise / 10.0)
                // Galactic noise (FaG)
            )
            .log10()
        - 10.0 * path.bw.log10()
        + 204.0; // kT₀ reference: -204 dBW/Hz

    // Step 4 & 7: "Signal upper decile deviation (day-to-day) (dB)" & "Signal lower decile deviation (day-to-day) (dB)"
    // Determine if the path crosses 60 degrees
    // From note 1 Table 2 P.842.4
    // If any point which "lies between control points located 1000 km from each end of the path" crosses
    // geomagnetic latitude 60 degrees then identify it.
    let mut gt60deg = 0; // 0 means less than 60 degrees
                         // Index for greater than 60 degrees geomagnetic latitude
    if path.distance > 2000.0 {
        // Note: There may be a degenerate case that looking at only three locations may not determine correctly
        // Check the mid-path
        geometry::geomagnetic_coords(path.cp[MP].l, &mut geomag);
        if geomag.lat >= 60.0 * D2R {
            gt60deg = 1; // 1 means a location is greater than 60 degrees
        }
        // Check the control point 1000 km from the transmitter
        geometry::geomagnetic_coords(path.cp[T1K].l, &mut geomag);
        if geomag.lat >= 60.0 * D2R {
            gt60deg = 1; // 1 means a location is greater than 60 degrees
        }
        // Check the control point 1000 km from the receiver
        geometry::geomagnetic_coords(path.cp[R1K].l, &mut geomag);
        if geomag.lat >= 60.0 * D2R {
            gt60deg = 1; // 1 means a location is greater than 60 degrees
        }
    }

    // Find the basic MUF index to retrieve the values from P.842-4 Table 2
    let f_bmuf_r = path.frequency / path.bmuf; // Frequency to basic MUF ratio

    let bmuf_idx = if f_bmuf_r <= 0.8 {
        0
    } else if f_bmuf_r <= 1.0 {
        1
    } else if f_bmuf_r <= 1.2 {
        2
    } else if f_bmuf_r <= 1.4 {
        3
    } else if f_bmuf_r <= 1.6 {
        4
    } else if f_bmuf_r <= 1.8 {
        5
    } else if f_bmuf_r <= 2.0 {
        6
    } else if f_bmuf_r <= 3.0 {
        7
    } else if f_bmuf_r <= 4.0 {
        8
    } else {
        9
    };

    // Deciles day-to-day
    let dl_sd = table2_ld[gt60deg][bmuf_idx]; // Signal lower decile deviation (day-to-day) (dB)
    let du_sd = table2_ud[gt60deg][bmuf_idx]; // Signal upper decile deviation (day-to-day) (dB)

    // Deciles hour-to-hour
    let du_sh = 5.0; // Signal upper decile deviation (hour-to-hour) (dB)
    let dl_sh = 8.0; // Signal lower decile deviation (hour-to-hour) (dB)

    // Step 6: "Upper decile deviation of resultant signal-to-noise ratio (dB)"
    let x = 10.0_f64.powf(noise_p.get_atmospheric_noise() / 10.0)
        + 10.0_f64.powf(noise_p.get_man_made_noise() / 10.0)
        + 10.0_f64.powf(noise_p.get_galactic_noise() / 10.0);

    let (du_a, dl_a) = noise_p.get_atmospheric_deciles();
    let (du_m, dl_m) = noise_p.get_man_made_deciles();
    let (du_g, dl_g) = noise_p.get_galactic_deciles();

    let y = 10.0_f64.powf((noise_p.get_atmospheric_noise() - dl_a) / 10.0)
        + 10.0_f64.powf((noise_p.get_man_made_noise() - dl_m) / 10.0)
        + 10.0_f64.powf((noise_p.get_galactic_noise() - dl_g) / 10.0);

    path.du_sn =
        ((10.0_f64 * (x / y).log10()).powi(2) + (du_sd as f64).powi(2) + (du_sh as f64).powi(2))
            .sqrt();

    // Step 9: "Lower decile deviation of resultant signal-to-noise ratio (dB)"
    // The value in variable x can be reused from Step 6 above.
    let y = 10.0_f64.powf((noise_p.get_atmospheric_noise() + du_a) / 10.0)
        + 10.0_f64.powf((noise_p.get_man_made_noise() + du_m) / 10.0)
        + 10.0_f64.powf((noise_p.get_galactic_noise() + du_g) / 10.0);

    path.dl_sn =
        ((10.0_f64 * (y / x).log10()).powi(2) + (dl_sd as f64).powi(2) + (dl_sh as f64).powi(2))
            .sqrt();

    // Step 11: "Basic circuit reliability for S/N >= or < S/Nr (%)"

    path.bcr = if path.snr >= path.snrr {
        let bcr_calc = f64::min(
            130.0 - 80.0 / (1.0 + (path.snr - path.snrr) / path.dl_sn),
            100.0,
        );
        bcr_calc
    } else {
        let bcr_calc = f64::max(
            80.0 / (1.0 + (path.snrr - path.snr) / path.du_sn) - 30.0,
            0.0,
        );
        bcr_calc
    };

    // End of BCR Calculation *********************************************************************

    // Begin simplified approximate BCR calculation for digital modulation systems ****************

    // The following calculation is based on a method found in ITU-R P.842-4
    // Section 9: "BCR for digital modulation systems".

    if path.modulation == DIGITAL && path.t0 != 0.0 && path.f0 != 0.0 {
        // Note: The following equivalences
        //		SNR0 = path.SNRr
        //		SNRm = path.SNR
        //		Dl = path.DlSN
        //		Du = path.DuSN
        // Time spread
        let d = path.distance; // for readability
                               // path distance
        let tm = if d <= 2000.0 {
            // Time spread
            f64::min(
                2.5e7 * (1.0 - (path.frequency / path.bmuf).powi(2)) * d.powf(-2.0),
                7.0 - 0.00175 * d,
            )
        } else {
            f64::min(
                4.27e-2 * (1.0 - (path.frequency / path.bmuf).powi(2)) * d.powf(0.65),
                3.5,
            )
        };

        // Frequency spread
        let fm = 0.02 * path.frequency * tm;

        // Time deciles
        let dt_u = 0.15 * tm; // Upper decile deviation of time spread (dB)
        let dt_l = 0.15 * tm; // Lower decile deviation of time spread (dB)

        // Frequency deciles
        let df_u = 0.1 * fm; // Upper decile deviation of frequency spread (dB)
        let df_l = 0.1 * fm; // Lower decile deviation of frequency spread (dB)

        // Probability that the required signal-to-noise ratio is achieved
        path.rsn = if path.snr >= path.snrr {
            f64::min(
                130.0 - 80.0 / (1.0 + (path.snr - path.snrr) / path.dl_sn),
                100.0,
            )
        } else {
            f64::max(
                80.0 / (1.0 + (path.snrr - path.snr) / path.du_sn) - 30.0,
                0.0,
            )
        };

        // Probability that the required time spread, T0, at a level of -10 dB relative to the peak
        // signal amplitude
        path.rt = if tm >= path.t0 {
            f64::min(130.0 - 80.0 / (1.0 + (path.t0 - tm) / dt_l), 100.0)
        } else {
            f64::max(80.0 / (1.0 + (tm - path.t0) / dt_u) - 30.0, 0.0)
        };

        // Probability that the required frequency spread f0 at a level of -10 dB relative to the peak
        // signal amplitude
        path.rf = if fm >= path.f0 {
            f64::min(130.0 - 80.0 / (1.0 + (path.f0 - fm) / df_l), 100.0)
        } else {
            f64::max(80.0 / (1.0 + (fm - path.f0) / df_u) - 30.0, 0.0)
        };
    }

    // End simplified approximate BCR calculation for digital modulation systems ******************

    // Begin OCR Calculation **********************************************************************

    if path.modulation == DIGITAL && path.distance <= 9000.0 {
        // Table 3 P.842-4
        // There are three sums for this calculation that can be accomplished in one loop:
        //		i) The interference sum with the protection ratio
        //		ii) The interference sum with the protection ratio and upper decile deviation
        //		iii) The interference sum with the protection ratio and lower decile deviation
        // First set all the deviations
        // Step 5: Set all the day-to-day deciles to 0 dB
        // Step 6: Upper decile deviation of the wanted signal, DuSh, was set to 5 dB above
        // the lower decile deviations of the interfering signal is set to 8 dB
        const DL_IH: f64 = 8.0; // Lower decile deviation of interference (dB)

        // Step 9: Lower decile deviation of the wanted signal, DlSh, was set to 8 dB above
        // the lower decile deviations of the interfering signal is set to 8 dB
        const DU_IH: f64 = 5.0; // Upper decile deviation of interference (dB)

        // Prepare the sums that will be used in the following steps:
        //		Step 4: Median resultant signal-to-interference signal (dB)
        //		Step 7: Upper decile deviation of resultant signal-to-interference ratio (dB)
        //		Step 10: Lower decile deviation of resultant signal-to-interference ratio (dB)
        // Initialize the sums
        let mut isum = 0.0; // The interference sum with the protection ratio
        let mut isumu = 0.0; // The interference sum with the protection ratio and upper decile deviation
        let mut isuml = 0.0; // The interference sum with the protection ratio and lower decile deviation

        // iI[] came from the routine DigitalModulationSignalandInterferers() which grouped
        // all E and F2 layers together. Consequently the E and F2 layers will have to be determined
        // separately here. Note: In the F2 layer loop the index is offset by 3 for the 3 E layer modes.
        // E layer loop
        for n in 0..MAX_MDS {
            if i_i[n] != NOTINDEX {
                if i_i[n] < MAX_E_MDS as i32 {
                    // E mode interference
                    isum += 10.0_f64.powf((path.md_e[i_i[n] as usize].prw - path.a) / 10.0);
                    isumu +=
                        10.0_f64.powf((path.md_e[i_i[n] as usize].prw - path.a + DU_IH) / 10.0);
                    isuml +=
                        10.0_f64.powf((path.md_e[i_i[n] as usize].prw - path.a - DL_IH) / 10.0);
                } else {
                    // F2 mode interference
                    isum += 10.0_f64.powf((path.md_f2[(i_i[n] - 3) as usize].prw - path.a) / 10.0);
                    isumu += 10.0_f64
                        .powf((path.md_f2[(i_i[n] - 3) as usize].prw - path.a + DU_IH) / 10.0);
                    isuml += 10.0_f64
                        .powf((path.md_f2[(i_i[n] - 3) as usize].prw - path.a - DL_IH) / 10.0);
                }
            }
        }

        // In the situation where all the modes are signal and there is no interference
        // there is nothing more to do
        // Return the values set by this program with the interference set to
        if isum == 0.0 {
            path.du_si = du_sh;
            path.dl_si = dl_sh;
            path.mir = 0.0;
            path.ocr = path.bcr * path.mir / 100.0;
            return;
        }

        // Step 4: Determine the signal-to-interference ratio
        path.sir = s - 10.0 * isum.log10();

        // Step 7: Determine the upper decile deviation of the signal-to-interference ratio
        path.du_si = (du_sh.powi(2) + (10.0 * (isum / isuml).log10()).powi(2)).sqrt();

        // Step 10: Determine the lower decile deviation of the signal-to-interference ratio
        path.dl_si = (dl_sh.powi(2) + (10.0 * (isumu / isum).log10()).powi(2)).sqrt();

        // Step 12: Circuit reliability in the presence of interference only for S/I >= or < S/Ir
        path.mir = if path.sir >= path.sirr {
            f64::min(
                130.0 - 80.0 / (1.0 + (path.sir - path.sirr) / path.dl_si),
                100.0,
            )
        } else {
            f64::max(
                80.0 / (1.0 + (path.sirr - path.sir) / path.du_si) - 30.0,
                0.0,
            )
        };

        // Overall Circuit reliability in the absence of scattering
        path.ocr = path.bcr * path.mir / 100.0;

        // Find the Overall Circuit reliability with scattering
        equatorial_scattering(path, &i_s);
    } else {
        // For analog modulation or long paths, set OCR = BCR
        path.ocr = path.bcr;
        path.ocrs = path.bcr;
    }

    // End OCR Calculation ************************************************************************

    // SNR for the required reliability
    //
    // For details please see
    // "CCIR Report 322 Noise Variation Parameters"
    // Technical Document 2813, June 1995
    // D. C. Lawrence
    // Naval Command, Control and Ocean Surveillance Center
    // RDT&E Division
    // http://www.dtic.mil/dtic/tr/fulltext/u2/a298722.pdf
    //
    // Note: The NORM[40] = 1.28
    //		 SNRXX = SNR50 +- t(XX%)*(D_u,l/1.28)

    if path.snrxxp < 50 {
        path.snrxx = path.snr + path.du_sn * norm[(50 - path.snrxxp) as usize] / norm[40];
    } else {
        // path.SNRXXp >= 50
        path.snrxx = path.snr - path.dl_sn * norm[(path.snrxxp - 50) as usize] / norm[40];
    }
}

/// Helper function to count the number of modes
fn number_of_modes(path: &PathData) -> i32 {
    let mut count = 0;

    // Count E layer modes
    for n in 0..MAX_E_MDS {
        if path.md_e[n].bmuf != 0.0 {
            count += 1;
        }
    }

    // Count F2 layer modes
    for n in 0..MAX_F2_MDS {
        if path.md_f2[n].bmuf != 0.0 {
            count += 1;
        }
    }

    count
}

/// Mode sorting function
fn mode_sort(m: &mut [crate::path_data::Mode], order: &mut [i32], criteria: i32) {
    /*
     ModeSort() is used to return the index to the mode of interest.

            INPUT
                struct Mode *M[MAXMDS] - 3 E modes and 6 F2 modes in one array
                int criteria - Flag that is either DOMINANT or SOONEST

            OUTPUT
                int order[MAXMDS] - Index array of the modes in the order desired

           SUBROUTINES
               None
    */

    // Initialize the order array
    let mut m_idx = 0; // This is so the array gets loaded from 0
    for n in 0..MAX_MDS {
        // Find the modes that exist
        if m[n].bmuf != 0.0 {
            order[m_idx] = n as i32;
            m_idx += 1;
        }
    }

    match criteria {
        DOMINANT => {
            for n in 0..MAX_MDS {
                if order[n] != NOTINDEX {
                    // The mode exists
                    for m_idx in 0..MAX_MDS {
                        if order[m_idx] != NOTINDEX {
                            // The mode exists
                            if m[order[n] as usize].prw < m[order[m_idx] as usize].prw {
                                let j = order[n];
                                order[n] = order[m_idx];
                                order[m_idx] = j;
                            }
                        }
                    }
                }
            }
        }
        SOONEST => {
            for n in 0..MAX_MDS {
                if order[n] != NOTINDEX {
                    // The mode exists
                    for m_idx in 0..MAX_MDS {
                        if order[m_idx] != NOTINDEX {
                            // The mode exists
                            if m[order[n] as usize].tau > m[order[m_idx] as usize].tau {
                                let j = order[n];
                                order[n] = order[m_idx];
                                order[m_idx] = j;
                            }
                        }
                    }
                }
            }
        }
        _ => {}
    }
}

/// Helper function to find F lambda d parameter
fn find_flambdad(cp: &crate::path_data::ControlPt) -> f64 {
    /*
     FindFlambdad() - Determines F sub lambda sub d in P.533-12 Appendix 1
            to Annex 1 "A model for scattering of HF signals"

        INPUT
            Control point to determine F sub lambda sub d

        OUTPUT
            returns F sub lambda sub d

        SUBROUTINES
           None
    */

    let lambdad = f64::abs(cp.dip[HR_100_KM]); // Magnetic dip parameter

    if lambdad >= 0.0 && lambdad < 15.0 * D2R {
        1.0
    } else if lambdad >= 15.0 * D2R && lambdad < 25.0 * D2R {
        f64::powi((25.0 - lambdad) / 10.0, 2) * ((lambdad - 10.0) / 5.0)
    } else if lambdad >= 25.0 * D2R && lambdad <= 90.0 * D2R {
        0.0
    } else {
        0.0
    }
}

/// Helper function to find F T l parameter
fn find_ftl(cp: &crate::path_data::ControlPt) -> f64 {
    /*
     FindTl() - Determines the F sub T sub l parameter in in P.533-12 Appendix 1
            to Annex 1 "A model for scattering of HF signals"

        INPUT
            Control point of interest

        OUTPUT
            returns time parameter F sub T sub l for the calculation of Prob sub occ

        SUBROUTINES
           None
    */

    let tl = cp.ltime; // Time parameter

    if tl > 0.0 && tl <= 3.0 {
        1.0
    } else if tl > 3.0 && tl <= 7.0 {
        f64::powi((7.0 - tl) / 4.0, 2) * ((tl - 1.0) / 2.0)
    } else if tl > 7.0 && tl <= 19.0 {
        0.0
    } else if tl > 19.0 && tl <= 20.0 {
        f64::powi(tl - 19.0, 2) * (41.0 - 2.0 * tl)
    } else if tl > 20.0 && tl <= 24.0 {
        1.0
    } else {
        0.0
    }
}
