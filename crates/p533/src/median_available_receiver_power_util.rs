use crate::constants::*;
use crate::median_skywave_field_strength_short_util::antenna_gain;
use crate::path_data::PathData;

/// Calculate median available receiver power
/// Based on ITU-R P.533-12 Section 6 "Median available receiver power"
pub fn median_available_receiver_power(path: &mut PathData) {
    /*
     MedianAvailableReceiverPower() - Calculates the median available receiver power by the method in
            ITU-R P.533-12 Section 6 "Median available receiver power".

            INPUT
                struct PathData *path

            OUTPUT
                path.Md_F2[].Prw - F2 mode received power
                path.Md_E[].Prw - E mode received power
                path.Pr - Median available received power
                path.Ep - The path field strength at the given path.distance.

           SUBROUTINE
               AntennaGain()
               DominantMode()
               AntennaGain08()
    */

    let grw: f64; // Mode gain paths < 9000 km and the overall receiver gain paths > 9000 km

    // Initialize
    let mut prw = TINY_DB as f64; // Greatest receive mode power
    let mut elevation = 2.0 * PI; // Antenna elevation

    match path.distance {
        // In each case the receiver gain, Grw, must be determined. Grw is calculated at the receiver elevation angles for
        // paths less than 9000 km. While for paths >= 9000 km the largest gain between 0 and 8 degrees is used.
        d if d <= 7000.0 => {
            // For each mode power to the total received power sum as given in Eqn (38) P.533-12
            // This can be done as each mode power is calculated
            let mut sum_pr = 0.0; // Summation of individual mode powers

            // Calculate the available signal power Prw (dBW) for each mode from
            // sky-wave field strength Ew (dB(1 µV/m)), frequency f (MHz) and Grw
            // lossless receiving antenna of gain.
            // Do any E-layer modes exist if so proceed
            // See "Modes considered" Section 5.2.1 P.533-12
            if path.n0_e != NO_LOWEST_MODE {
                for i in path.n0_e as usize..MAX_E_MDS {
                    if (i == path.n0_e as usize
                        && path.distance / (path.n0_e as f64 + 1.0) <= 2000.0)
                        || (i != path.n0_e as usize && path.md_e[i].bmuf != 0.0)
                    {
                        // Find the receiver gain for this mode.
                        path.md_e[i].grw = antenna_gain(path, &path.a_rx, path.md_e[i].ele, RXTOTX);

                        path.md_e[i].prw = path.md_e[i].ew + path.md_e[i].grw
                            - 20.0 * f64::log10(path.frequency)
                            - 107.2;

                        // Determine if this is the greatest received power
                        // If this is first time in the loop i == path.n0_E then initialize Prw
                        if prw < path.md_e[i].prw {
                            // Prw is the greatest power
                            prw = path.md_e[i].prw;

                            // Point to the dominant mode and set the dominant mode index.
                            path.dm_ptr = Some(path.md_e[i].clone());
                            path.dm_idx = i as i32;
                        }

                        // Add this mode to the sum.
                        sum_pr += f64::powf(10.0, path.md_e[i].prw / 10.0);
                    }
                }
            }

            // F2 modes
            // Do any F2-layer modes exist if so proceed
            if path.n0_f2 != NO_LOWEST_MODE {
                for i in path.n0_f2 as usize..MAX_F2_MDS {
                    if (i == path.n0_f2 as usize
                        && path.distance / (path.n0_f2 as f64 + 1.0) <= path.dmax
                        && path.md_f2[i].fs < path.frequency)
                        || (i != path.n0_f2 as usize
                            && path.md_f2[i].bmuf != 0.0
                            && path.md_f2[i].fs < path.frequency)
                    {
                        // Find the receiver gain for this mode.
                        path.md_f2[i].grw =
                            antenna_gain(path, &path.a_rx, path.md_f2[i].ele, RXTOTX);

                        path.md_f2[i].prw = path.md_f2[i].ew + path.md_f2[i].grw
                            - 20.0 * f64::log10(path.frequency)
                            - 107.2;

                        // Determine if this is the greatest received power.
                        // If there was an E mode then Prw is already set to that power
                        if prw < path.md_f2[i].prw {
                            // Prw is the greatest power.
                            prw = path.md_f2[i].prw;

                            // Point to the dominant mode and set the dominant mode index.
                            path.dm_ptr = Some(path.md_f2[i].clone());
                            path.dm_idx = (i + 3) as i32;
                        }

                        // Add this mode to the sum.
                        sum_pr += f64::powf(10.0, path.md_f2[i].prw / 10.0);
                    }
                }
            }

            // Now that the modes are calculated, set the path parameters.
            // Find the total received power.
            // If the SumPr is 0 then set the path.Pr to something small
            if sum_pr != 0.0 && path.dm_ptr.is_some() {
                path.pr = 10.0 * f64::log10(sum_pr);
                // The dominant mode is known so set any values in the path structure that are relevant.
                dominant_mode(path);
            } else {
                // There are no E modes and all the F2 modes are screened
                path.pr = TINY_DB as f64;
            }

            // Save the path field strength. For this distance select Es, which includes E layer screening.
            path.ep = path.es;
        }

        d if d > 7000.0 && d < 9000.0 => {
            // Determine the receiver gain.
            grw = antenna_gain_08(path, &path.a_rx, RXTOTX, &mut elevation);

            // Use the interpolated power, Ei.
            path.pr = path.ei + grw - 20.0 * f64::log10(path.frequency) - 107.2;

            // The path receiver gain Grw.
            path.grw = grw;

            // Save the path field strength, which at this distance is Ei.
            path.ep = path.ei;

            // Set the rx antenna elevation to the  long path rx elevation
            path.ele = elevation;
        }

        _ => {
            // path.distance >= 9000.0
            // Determine the receiver gain.
            grw = antenna_gain_08(path, &path.a_rx, RXTOTX, &mut elevation);

            // Use the combined mode power, El, and the antenna gain between 0 and 8 degrees, Grw.
            path.pr = path.el + grw - 20.0 * f64::log10(path.frequency) - 107.2;

            // The path receiver gain Grw.
            path.grw = grw;

            // Save the path field strength, which at this distance is El.
            path.ep = path.el;

            // Set the rx antenna elevation to the  long path rx elevation
            path.ele = elevation;
        }
    }
}

/// Stores values that are associated with the dominant mode to the path structure
fn dominant_mode(path: &mut PathData) {
    /*
     DominentMode() - Stores values that are associated with the dominant
            to the path structure. These are strictly not part of the standard
            P.533-12 but are provided for continuity in the analysis.

            INPUT
                struct PathData *path

            OUTPUT
                path->Grw
                path->ele

           SUBROUTINES
               None
    */

    // Select the dominant mode to represent the median path behavior.
    // The path receiver gain is the dominant mode receiver gain.
    if let Some(dm_ptr) = &path.dm_ptr {
        path.grw = dm_ptr.grw;

        // The path elevation angle is the dominant mode elevation angle.
        path.ele = dm_ptr.ele;
    }
}

/// Antenna gain 08 function
/// Determines the largest antenna gain in the range 0 to 8 degrees elevation
/// This is an EXACT translation of MedianSkywaveFieldStrengthLongUtil.AntennaGain08
fn antenna_gain_08(
    path: &PathData,
    ant: &crate::path_data::Antenna,
    direction: i32,
    elevation: &mut f64,
) -> f64 {
    /*
        AntennaGain08() - Determines the largest antenna gain in the range 0 to 8 degrees elevation.
            At present June 2013. There is a minimum elevation angle in the long model that is fixed
            at 3 degrees. Per Dambolt and Suessmann this routine which finds the maximum gain between
            0 and 8 degrees should not be altered under the assumption that it is improbable that the
            antenna gain determined by the proceedure would be less than 3 degrees.

            INPUT
                struct PathData path
                struct Antenna Ant
                int direction

            OUTPUT
                largest antenna gain in the range 0 to 8 degrees elevation

            SUBROUTINES
                AntennaGain()
    */

    // The structure Antenna Ant is used to tell the subroutine which antenna to calculate.

    // Find the largest value of transmitting gain at the required azimuth in the elevation
    // range 0 to 8 degrees
    let mut gmax = TINY_DB as f64;
    for i in 0..9 {
        let delta = i as f64 * D2R; // Elevation angle
        let g = crate::median_skywave_field_strength_short_util::antenna_gain(
            path, ant, delta, direction,
        );
        if g > gmax {
            gmax = g;
            *elevation = i as f64 * D2R;
        }
    }

    gmax
}
