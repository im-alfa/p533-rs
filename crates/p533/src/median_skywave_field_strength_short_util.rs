use crate::constants::*;
use crate::path_data::{PathData, ControlPt, Antenna, Location};
use crate::e_layer_screening_frequency::{elevation_angle, incidence_angle};

const TXEND: usize = 0;
const RXEND: usize = 1;
const NOIL: f64 = 9.14;

/// Calculate the median skywave field strength as described in ITU-R P.533-12 
/// Section 5.2 paths up to 7000 km
pub fn median_skywave_field_strength_short(path: &mut PathData) {
    // Only do this subroutine if the path is less than or equal to 9000 km if not exit
    if path.distance > 9000.0 {
        return;
    }
    
    // 5.2 Paths up to 7000 km
    // 5.2.1 Modes considered
    //		The modes considered will be determined by if statements within the 
    //		E and F2 layer calculation loops
    // 5.2.2 Field strength determination
    // For this calculation the SSN is restricted to 160 
    let ssn = f64::min(path.ssn as f64, MAX_SSN as f64); // Sun spot number
    
    // This procedure applies only to paths less than 7000 km where it is the only method used and
    // paths greater than 7000 km but less than 9000 km where the field strength is interpolated with
    // the long path method
    // Initialize the mirror reflection height for the E-Layer 
    const HR_E: f64 = 110.0;
    
    // Determine the mirror reflection height for the F2 layers
    // Use hr_F2 in this routine for readability
    // Note the path.distance is less than 9000 and path.distance is greater than dmax
    let hr_f2 = if path.distance > path.dmax {
        f64::min(1490.0 / path.cp[smallest_cp_fof2(path)].m3kf2 - 176.0, 500.0) // Find the smallest foF2 amongst the control points
    } else {
        path.cp[MP].hr // Use the midpoint to determine the mirror reflection height
    };
    
    /*******************************************************************************************************/
    // Although the calculations for E and F2 layers are for the most part the 
    // same in order to not obscure the calculation relative to the standard
    // each layer will be dealt within its own loop. 
    /*******************************************************************************************************/
    // Determine the local time at the midpath point
    let tz = (path.cp[MP].l.lng / (15.0 * D2R)) as i32; // Time zone at midpath
    let mpltime = (path.cp[MP].ltime + (tz % 24) as f64) as i32; // Midpath local time
    
    // Begin E modes median sky-wave field strength calculation
    // Does a low order E mode exist?
    if path.n0_e != NO_LOWEST_MODE {
        // Determine the E modes that satisfy the criteria.
        for n in path.n0_e as usize..MAX_E_MDS {
            if (n == path.n0_e as usize && path.distance / (path.n0_e as f64 + 1.0) <= 2000.0) // The lowest order mode is less than 2000 km.
                || (n > path.n0_e as usize && path.md_e[n].bmuf != 0.0) // higher order modes
            {
                // Find the elevation angle from equation 13 Section 5.1 Elevation angle
                // ITU-R P.533-12
                let delta = elevation_angle(path.distance / (n as f64 + 1.0), HR_E);
                
                // Store the elevation angle, but use delta elsewhere here for readability.
                path.md_e[n].ele = delta;
                
                // angle of incidence at height hr = 110 km
                let aoi110 = incidence_angle(delta, HR_E);
                
                // Vertical-incidence wave frequency
                let fv = path.frequency * f64::cos(aoi110);
                
                // Hop distance 
                let dh = path.distance / (n as f64 + 1.0);
                
                // Find the virtual slant range (19).
                let psi = dh / (2.0 * R0);
                
                path.ptick = f64::abs(2.0 * R0 * (f64::sin(psi) / f64::cos(delta + psi))) * (n as f64 + 1.0);
                
                // Determine the number of control points.
                let (at, fl, lh) = if path.distance <= 2000.0 {
                    // Find the loss due to all the absorption terms in Li
                    // The absorption term includes loss from solar zenith angles, ATnoon and phin(fv/foE)
                    let at = absorption_term(&path.cp[MP], path.month, fv);
                    
                    // Determine the longitudinal gyrofrequency
                    let fl = f64::abs(path.cp[MP].fh[HR_100_KM] * f64::sin(path.cp[MP].dip[HR_100_KM]));
                    
                    // Determine auroral and other signal losses
                    let lh = find_lh(&path.cp[MP], dh, mpltime, path.month);
                    
                    (at, fl, lh)
                } else {
                    // (path.distance > 2000.0) There are three control points
                    
                    // Find the loss due to all the absorption terms in Li
                    // The absorption term includes loss from solar zenith angles, ATnoon and phin(fv/foE)
                    let at = (absorption_term(&path.cp[MP], path.month, fv) +
                              absorption_term(&path.cp[T1K], path.month, fv) +
                              absorption_term(&path.cp[R1K], path.month, fv)) / 3.0;
                    
                    // Determine the average longitudinal gyrofrequency
                    let fl = (f64::abs(path.cp[MP].fh[HR_100_KM] * f64::sin(path.cp[MP].dip[HR_100_KM])) +
                              f64::abs(path.cp[T1K].fh[HR_100_KM] * f64::sin(path.cp[T1K].dip[HR_100_KM])) +
                              f64::abs(path.cp[R1K].fh[HR_100_KM] * f64::sin(path.cp[R1K].dip[HR_100_KM]))) / 3.0;
                    
                    // Determine auroral and other signal losses
                    let lh = (find_lh(&path.cp[MP], dh, mpltime, path.month) +
                              find_lh(&path.cp[T1K], dh, mpltime, path.month) +
                              find_lh(&path.cp[R1K], dh, mpltime, path.month)) / 3.0;
                    
                    (at, fl, lh)
                };
                
                // All the variable have been calculated to determine
                // Absorption loss (dB) for an n-hop mode, Li
                let li = (n as f64 + 1.0) * (1.0 + 0.0067 * ssn) * at / (f64::powf(path.frequency + fl, 2.0) * f64::cos(aoi110));
                
                // "Above-the-MUF" loss
                let lm = if path.frequency <= path.md_e[n].bmuf {
                    0.0
                } else {
                    f64::min(46.0 * f64::powf(path.frequency / path.md_e[n].bmuf - 1.0, 0.5) + 5.0, 58.0)
                };
                
                // Ground reflection loss
                let lg = 2.0 * (n as f64 + 1.0 - 1.0); // n is a C index starting at 0 instead of 1 
                
                // "Not otherwise included" loss
                path.lz = NOIL;
                
                // The ray path basic transmission loss for the mode under consideration
                path.md_e[n].lb = 32.45 + 20.0 * f64::log10(path.frequency) + 20.0 * f64::log10(path.ptick) + li + lm + lg + lh + path.lz;
                
                // Tx antenna gain in the desired direction (dB)
                let gt = antenna_gain(path, &path.a_tx, delta, TXTORX);
                
                // Transmit power
                let pt = path.txpower;
                
                path.md_e[n].ew = 136.6 + pt + gt + 20.0 * f64::log10(path.frequency) - path.md_e[n].lb;
            } else {
                break; // There is no E modes that satisfiy the criteria
            }
        }
    }
    
    /******************************************************************************************************************/
    // Begin F2 modes median sky-wave field strength calculation
    // Does a low order F2 mode exist?
    if path.n0_f2 != NO_LOWEST_MODE {
        // Determine the F2 mode that satisfy the criteria.
        for n in path.n0_f2 as usize..MAX_F2_MDS {
            // Modes Considered:
            // Any mode that has an E-layer maximum screening frequency that is less than the operating frequency. 
            // The lowest order mode must also have hops that are less than dmax km 
            if (n == path.n0_f2 as usize && path.distance / (path.n0_f2 as f64 + 1.0) <= path.dmax && path.md_f2[n].fs < path.frequency)
                || (n > path.n0_f2 as usize && path.md_f2[n].bmuf != 0.0 && path.md_f2[n].fs < path.frequency)
            {
                // higher order modes

                // Find the elevation angle from equation 13 Section 5.1 Elevation angle.
                // ITU-R P.533-12
                let delta = elevation_angle(path.distance / (n as f64 + 1.0), hr_f2);

                // Store the elevation angle
                path.md_f2[n].ele = delta;

                // angle of incidence at height hr = 110 km
                let aoi110 = incidence_angle(delta, 110.0);

                // Vertical-incidence wave frequency
                let fv = path.frequency * f64::cos(aoi110);

                // Hop distance 
                let dh = path.distance / (n as f64 + 1.0);

                // Find the virtual slant range (19)
                let psi = dh / (2.0 * R0);

                // Calculate the slant range
                path.ptick = f64::abs(2.0 * R0 * (f64::sin(psi) / f64::cos(delta + psi))) * (n as f64 + 1.0);

                // Determine the number of control points
                let (at, fl, lh) = if path.distance <= 2000.0 {
                    // Use the mid-path control point

                    // Find the loss due to all the absorption terms in Li
                    // The absorption term includes loss from solar zenith angles, ATnoon and phin(fv/foE)
                    let at = absorption_term(&path.cp[MP], path.month, fv);

                    // Determine the longitudinal gyrofrequency
                    let fl = f64::abs(path.cp[MP].fh[HR_100_KM] * f64::sin(path.cp[MP].dip[HR_100_KM]));

                    // Determine auroral and other signal losses
                    let lh = find_lh(&path.cp[MP], dh, mpltime, path.month);

                    (at, fl, lh)
                } else if 2000.0 < path.distance && path.distance <= path.dmax {
                    // There are three control points

                    // Find the loss due to all the absorption terms in Li
                    // The absorption term includes loss from solar zenith angles, ATnoon and phin(fv/foE)
                    let at = (absorption_term(&path.cp[MP], path.month, fv) +
                              absorption_term(&path.cp[T1K], path.month, fv) +
                              absorption_term(&path.cp[R1K], path.month, fv)) / 3.0;

                    // Determine the average longitudinal gyrofrequency
                    let fl = (f64::abs(path.cp[MP].fh[HR_100_KM] * f64::sin(path.cp[MP].dip[HR_100_KM])) +
                              f64::abs(path.cp[T1K].fh[HR_100_KM] * f64::sin(path.cp[T1K].dip[HR_100_KM])) +
                              f64::abs(path.cp[R1K].fh[HR_100_KM] * f64::sin(path.cp[R1K].dip[HR_100_KM]))) / 3.0;

                    // Determine auroral and other signal losses
                    let lh = (find_lh(&path.cp[MP], dh, mpltime, path.month) +
                              find_lh(&path.cp[T1K], dh, mpltime, path.month) +
                              find_lh(&path.cp[R1K], dh, mpltime, path.month)) / 3.0;

                    (at, fl, lh)
                } else {
                    // There are 5 control points.

                    // Find the loss due to all the absorption terms in Li.
                    // The absorption term includes loss from solar zenith angles, ATnoon and phin(fv/foE)
                    let at = (absorption_term(&path.cp[MP], path.month, fv) +
                              absorption_term(&path.cp[T1K], path.month, fv) +
                              absorption_term(&path.cp[R1K], path.month, fv) +
                              absorption_term(&path.cp[TD02], path.month, fv) +
                              absorption_term(&path.cp[RD02], path.month, fv)) / 5.0;

                    // Find the average longitudinal gyrofrequency
                    let fl = (f64::abs(path.cp[MP].fh[HR_100_KM] * f64::sin(path.cp[MP].dip[HR_100_KM])) +
                              f64::abs(path.cp[T1K].fh[HR_100_KM] * f64::sin(path.cp[T1K].dip[HR_100_KM])) +
                              f64::abs(path.cp[R1K].fh[HR_100_KM] * f64::sin(path.cp[R1K].dip[HR_100_KM])) +
                              f64::abs(path.cp[TD02].fh[HR_100_KM] * f64::sin(path.cp[TD02].dip[HR_100_KM])) +
                              f64::abs(path.cp[RD02].fh[HR_100_KM] * f64::sin(path.cp[RD02].dip[HR_100_KM]))) / 5.0;

                    // Determine auroral and other signal losses
                    let lh = (find_lh(&path.cp[MP], dh, mpltime, path.month) +
                              find_lh(&path.cp[T1K], dh, mpltime, path.month) +
                              find_lh(&path.cp[R1K], dh, mpltime, path.month) +
                              find_lh(&path.cp[TD02], dh, mpltime, path.month) +
                              find_lh(&path.cp[RD02], dh, mpltime, path.month)) / 5.0;

                    (at, fl, lh)
                };

                // All the variable have been calculated to determine
                // Absorption loss (dB) for an n-hop mode, Li.
                let li = (n as f64 + 1.0) * (1.0 + 0.0067 * ssn) * at / (f64::powf(path.frequency + fl, 2.0) * f64::cos(aoi110));

                // "Above-the-MUF" loss
                let lm = if path.frequency <= path.md_f2[n].bmuf {
                    0.0
                } else {
                    // (path.frequency > path.md_f2[n].bmuf)
                    if path.distance <= 3000.0 {
                        f64::min(36.0 * f64::powf(path.frequency / path.md_f2[n].bmuf - 1.0, 0.5) + 5.0, 60.0)
                    } else {
                        f64::min(70.0 * (path.frequency / path.md_f2[n].bmuf - 1.0) + 8.0, 80.0)
                    }
                };

                // Ground reflection loss
                let lg = 2.0 * (n as f64 + 1.0 - 1.0); // n is a C index starting at 0 instead of 1  

                // "Not otherwise included" loss
                path.lz = NOIL;

                // Calculate the total loss.
                path.md_f2[n].lb = 32.45 + 20.0 * f64::log10(path.frequency) + 20.0 * f64::log10(path.ptick) + li + lm + lg + lh + path.lz;

                // Tx antenna gain in the desired direction (dB)
                let gt = antenna_gain(path, &path.a_tx, delta, TXTORX);

                let pt = path.txpower;

                path.md_f2[n].ew = 136.6 + pt + gt + 20.0 * f64::log10(path.frequency) - path.md_f2[n].lb;
            }
        }
    }

    // Determine the overall resultant equivalent median sky-wave field strength, Es
    // See "Modes considered" Section 5.2.1 P.533-12
    // Es should be very small, reinitialize it for clarity.
    path.es = TINY_DB as f64;
    // Initialize the other locals needed here.
    let mut etw = 0.0; // Median field strength

    // Do any E-layer modes exist if so proceed
    if path.n0_e != NO_LOWEST_MODE {
        // Sum the lowest-order E-layer mode with hop length up to 2000 km
        // and the next two higher-order E-layer modes that exist;
        for n in path.n0_e as usize..MAX_E_MDS {
            if (n == path.n0_e as usize && path.distance / (path.n0_e as f64 + 1.0) <= 2000.0)
                || (n != path.n0_e as usize && path.md_e[n].bmuf != 0.0)
            {
                let mode_power = f64::powf(10.0, path.md_e[n].ew / 10.0);
                etw += mode_power;
                path.md_e[n].mc = TRUE;
            }
        }
    }

    // Do any F2-layer modes exist if so proceed
    if path.n0_f2 != NO_LOWEST_MODE {
        // Mode Considered: The lowest-order F2-layer mode with a hop length up to dmax (km) and 
        // maximally the next five F2-layer higher-order modes.
        // If that mode's E layer screening frequency is less than the operating frequency
        for n in path.n0_f2 as usize..MAX_F2_MDS {
            if (n == path.n0_f2 as usize && path.distance / (path.n0_f2 as f64 + 1.0) <= path.dmax && path.md_f2[n].fs < path.frequency)
                || (n != path.n0_f2 as usize && path.md_f2[n].bmuf != 0.0 && path.md_f2[n].fs < path.frequency)
            {
                let mode_power = f64::powf(10.0, path.md_f2[n].ew / 10.0);
                etw += mode_power;
                path.md_f2[n].mc = TRUE;
            }
        }
    }

    // Find the field strength if there are any modes to consider
    // If there are no modes than path.es will remain equal to TINY_DB 
    if etw != 0.0 {
        path.es = 10.0 * f64::log10(etw); // Field strength with E layer screening
    }
}

/// Determines the smallest Control point foF2
fn smallest_cp_fof2(path: &PathData) -> usize {
    let mut temp = 0;
    let mut idx = [0, 1, 2, 3, 4]; // This is what will change.
    
    // Sort by brute force
    for i in 0..5 {
        for j in 0..5 {
            if path.cp[idx[i]].fof2 > path.cp[idx[j]].fof2 {
                let temp_val = idx[i];
                idx[i] = idx[j];
                idx[j] = temp_val;
            }
        }
    }
    
    // return the last non zero index
    for i in 0..5 {
        if idx[i] != 0 {
            temp = idx[i];
        }
    }
    
    temp
}

/// Determines the absorption term from the three factors
fn absorption_term(cp: &ControlPt, month: i32, fv: f64) -> f64 {
    /*
         AbsorptionTerm() - Determines the absorption term from the three factors:
             i) Absorption factor at local noon
             ii) Absorption layer penetration factor
             iii) Diurnal absorption exponent
             These three factors can be seen in Figures 1, 2, and 3 of P.533-12 
             Section 5.2.2 "Field strength determination". The absorption factor is the product of these
             three factors
     
             INPUT
                 struct ControlPt CP - the Control point of interest
                 int month - The month index
                 double fv - Vertical-incidence wave frequency
     
             OUTPUT
                 returns the absorption term used to calculate Li in Eqn (20)
     
             SUBROUTINES
                ZeroCP()
                DiurnalAbsorptionExponent()
                SolarParameters()
                AbsorptionFactor()
                AbsorptionLayerPenetrationFactor()
     */

    let mut cp_0 = cp.clone(); // Temp control point

    // Find the diurnal absorption exponent, p.
    let p = diurnal_absorption_exponent(cp, month); // Diurnal absorption exponent

    // The solar zenith angle for the control point
    // Make sure that it doesn't exceed 102 degrees
    let chij = f64::min(cp.sun.sza, 102.0 * D2R); // Solar zenith angle at the j-th control point

    let fchij = f64::max(f64::powf(f64::cos(0.881 * chij), p), 0.02); // Eqn (21) with chij argument

    // Determine when noon is local time, then find the solar zenith angle.
    // To make this calculation, a temporary control point must be used 
    // where the only variable necessary to initialize is the longitude.
    cp_0 = cp.clone();
    // The hour that gets passed to SolarParameters() is a UTC fractional hour 
    // The local noon in UTC has already been calculated and stored by execution 
    // CalculateCPParameters()
    // The month is in the path structure thus: 
    // The hour for this calculation is CP[].Sun.lsn and
    // the month for this calculation is path->month.
    crate::calculate_cp_parameters::solar_parameters(&mut cp_0, month, cp.sun.lsn);

    // The solar zenith angle for the control point at noon local time.
    let chijnoon = cp_0.sun.sza; // Solar zenith angle at the j-th control point at noon

    let fchijnoon = f64::max(f64::powf(f64::cos(0.881 * chijnoon), p), 0.02); // Eqn (21) with chijnoon argument

    // Find the remaining absorption parameters.
    let atnoon = absorption_factor(cp, month); // Absorption factor at local noon and R12 = 0

    let phin = absorption_layer_penetration_factor(fv / cp.foe); // Absorption layer penetration factor

    atnoon * phin * fchij / fchijnoon
}

/// Determine the diurnal absorption exponent, p
fn diurnal_absorption_exponent(cp: &ControlPt, mut month: i32) -> f64 {
    /*
          DiurnalAbsorptionExponent() Determine the diurnal absorption exponent, p. 
             The exponent, p, is a function of the month and magnetic dip angle.
             The p vs magnetic dip angle graph is shown as Figure 3 ITU-R P.533-12.
     
             INPUT
                 struct ControlPt CP - control point of interest
                 int month - month index	
     
             OUTPUT
                 return p the diurnal absorption exponent
     
            SUBROUTINES
                None

         This is based on CONTP() in REC533
     
        For latitude = 3.4464, magnetic dip = 2.89 and month = 4 the diurnal 
        absorption exponent will be	p = 0.689             
     */

    let ppt = [30.0, 30.0, 30.0, 27.5, 32.5, 35.0, 37.5, 35.0, 32.5, 30.0, 30.0, 30.0];

    let pval1: [[[f64; 7]; 2]; 6] = [
        [
            [1.510, -0.353, -0.090, 0.191, 0.133, -0.067, -0.053],
            [1.400, -0.365, -1.212, -0.049, 1.187, 0.119, -0.400]
        ],
        [
            [1.490, -0.348, -0.055, 0.164, 0.160, -0.041, -0.080],
            [1.450, -0.119, -0.913, -0.640, 0.347, 0.458, 0.107]
        ],
        [
            [1.520, -0.410, -0.138, 0.308, 0.267, -0.113, -0.133],
            [1.500, -0.492, -0.958, 0.216, 0.267, -0.029, 0.187]
        ],
        [
            [1.580, -0.129, -0.228, -0.192, 0.200, 0.116, -0.027],
            [1.530, -0.468, -1.312, 0.096, 0.973, 0.057, -0.187]
        ],
        [
            [1.590, 0.002, -0.102, -0.579, -0.467, 0.522, 0.613],
            [1.490, -0.937, -1.622, 1.365, 1.720, -0.873, -0.453]
        ],
        [
            [1.600, -0.060, -0.175, -0.037, 0.147, -0.008, -0.027],
            [1.460, -0.881, -1.595, 0.901, 2.133, -0.395, -0.933]
        ]
    ];

    let pval2: [[[f64; 7]; 2]; 6] = [
        [
            [1.60, -0.030, -0.135, -0.137, 0.053, 0.072, 0.027],
            [1.43, -0.902, -1.667, 0.905, 2.480, -0.383, -1.173]
        ],
        [
            [1.59, -0.032, -0.083, -0.119, 0.000, 0.031, 0.053],
            [1.46, -0.831, -1.653, 0.708, 2.320, -0.257, -1.067]
        ],
        [
            [1.59, -0.060, -0.180, -0.181, 0.267, 0.081, -0.107],
            [1.51, -0.809, -1.740, 0.750, 2.240, -0.301, -0.960]
        ],
        [
            [1.57, -0.189, -0.207, -0.005, 0.293, 0.004, -0.107],
            [1.52, -0.433, -1.015, -0.017, 0.440, 0.115, 0.080]
        ],
        [
            [1.55, -0.292, -0.275, 0.093, 0.427, -0.026, -0.187],
            [1.44, -0.279, -0.770, -0.266, 0.053, 0.245, 0.267]
        ],
        [
            [1.51, -0.347, -0.082, 0.160, 0.093, -0.048, -0.027],
            [1.40, -0.355, -1.212, -0.102, 1.187, 0.172, -0.400]
        ]
    ];

    // Initialize the exponent p
    let mut p = 0.0;

    // Initialize the modified magnetic dip angle (degrees)
    let mut moddip = f64::abs(f64::atan2(cp.dip[HR_100_KM], f64::sqrt(f64::cos(cp.l.lat)))); // Modified magnetic dip (or latitude) (degrees)

    if moddip > 70.0 * D2R {
        moddip = 70.0 * D2R;
    }

    if cp.l.lat < 0.0 {
        month += 6;
        if month > 11 {
            month -= 12;
        }
    }

    let pp = ppt[month as usize] * D2R;

    let i = if moddip > pp {
        moddip = -1.0 + 2.0 * (moddip - pp) / (70.0 * D2R - pp);
        1
    } else {
        moddip = -1.0 + 2.0 * moddip / pp;
        0
    };

    let mut sx = 1.0;

    for j in 0..7 {
        let a = if month <= 5 {
            pval1[month as usize][i][j]
        } else {
            pval2[(month - 6) as usize][i][j]
        };

        p += a * sx;
        sx *= moddip;
    }

    p
}

/// Calculates the absorption factor ATnoon as shown Figure 1 ITU-R P.533-12
fn absorption_factor(cp: &ControlPt, month: i32) -> f64 {
    /*
          AbsorptionFactor() Calculates the absorption factor ATnoon as shown Figure 1 ITU-R P.533-12 
         
                 INPUT
                     struct ControlPt CP - Control point of interest
                     int month - month index	
         
                 OUTPUT
                     returns ATnoon - The absorption factor at local noon and R12 = 0.
         
                 SUBROUTINES
                    None
         
         This routine is based on CONTAT() in REC533.             
     */

    let atno: [[f64; 29]; 9] = [
        [
            323.9, 297.5, 274.5, 256.4, 244.2, 235.0, 229.5, 226.1, 226.8, 229.0, 232.5, 237.0, 243.4, 249.9, 258.1, 267.5, 277.5, 283.3, 283.2, 273.1, 257.0, 232.1, 201.4, 171.5, 146.0, 123.0, 103.1, 83.0, 66.6
        ],
        [
            312.1, 285.1, 263.1, 251.8, 249.5, 250.9, 254.5, 260.3, 266.7, 272.3, 277.8, 280.3, 283.9, 284.5, 284.4, 283.0, 278.6, 273.0, 265.7, 256.3, 244.8, 232.0, 218.1, 204.5, 189.9, 172.3, 155.3, 135.5, 116.2
        ],
        [
            347.7, 321.9, 302.5, 293.8, 291.4, 289.3, 292.1, 296.6, 304.3, 313.0, 321.7, 333.8, 342.6, 349.6, 355.2, 355.6, 352.2, 341.7, 327.3, 308.4, 286.0, 265.0, 244.1, 223.8, 202.8, 181.8, 160.8, 141.6, 123.4
        ],
        [
            338.0, 313.2, 297.0, 290.2, 292.1, 299.4, 308.0, 320.4, 331.6, 340.7, 347.8, 353.8, 357.0, 360.0, 359.8, 358.3, 355.8, 350.8, 344.5, 332.7, 316.4, 292.5, 266.1, 236.4, 214.0, 193.8, 177.5, 165.0, 155.9
        ],
        [
            328.1, 303.8, 287.7, 282.5, 284.4, 289.4, 294.8, 303.6, 312.9, 322.7, 332.3, 343.8, 350.6, 358.7, 364.3, 365.8, 362.4, 356.0, 346.7, 333.0, 318.8, 299.7, 282.1, 260.5, 240.5, 220.6, 203.9, 186.3, 173.0
        ],
        [
            305.1, 288.5, 275.2, 273.7, 278.6, 288.9, 302.5, 319.3, 333.6, 346.3, 356.3, 364.7, 371.7, 373.6, 374.2, 373.1, 370.5, 365.1, 358.5, 347.7, 335.0, 320.3, 299.1, 276.6, 253.2, 230.7, 214.0, 196.6, 185.3
        ],
        [
            345.4, 319.4, 298.7, 290.1, 290.0, 291.8, 296.3, 302.9, 312.1, 320.1, 327.8, 334.1, 340.2, 343.3, 345.7, 346.5, 345.3, 341.1, 334.5, 321.7, 304.2, 286.8, 265.9, 244.8, 224.1, 204.5, 183.6, 164.1, 145.2
        ],
        [
            341.9, 314.8, 295.3, 277.9, 265.0, 258.2, 254.4, 255.8, 257.3, 262.9, 268.5, 279.0, 287.5, 295.2, 299.6, 300.2, 298.9, 291.5, 279.0, 262.6, 245.7, 227.0, 203.6, 182.3, 163.2, 147.1, 133.9, 119.9, 110.8
        ],
        [
            318.8, 293.3, 268.3, 251.7, 240.4, 233.1, 229.4, 228.8, 230.5, 235.5, 239.7, 242.6, 245.4, 247.5, 248.9, 249.9, 248.5, 244.4, 237.3, 225.6, 213.5, 195.2, 172.7, 151.3, 131.1, 113.1, 100.1, 89.0, 80.0
        ]
    ];

    /*
        Note: the month index into the array is as follows:
                Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec
            i =  0   1   2   3   4   5   5   4   6   7   8   0
     */

    let i = match month {
        JUL => 5,
        AUG => 4,
        SEP => 6,
        OCT => 7,
        NOV => 8,
        DEC => 0,
        _ => month
    };

    let mut x = f64::abs(cp.l.lat * R2D); // ?
    if x >= 70.0 {
        x = 69.99; // This is for the (int) casting of X so that j is not >= 28.
    }
    x /= 2.5;
    let j = x as usize; // Index
    x -= j as f64;
    let atnoon = atno[i as usize][j + 1] * x + atno[i as usize][j] * (1.0 - x); // Absorption factor at local noon and R12 = 0 

    atnoon
}

/// Determines the absorption layer penetration factor shown in Figure 2 P.533-12
fn absorption_layer_penetration_factor(t: f64) -> f64 {
    /*
         AbsorptionLayerPenetrationFactor() - Determines the absorption layer penetration factor
             shown in Figure 2 P.533-12.
     
             INPUT
                 T - the ratio of the Vertical-incidence wave frequency and foE
     
             OUTPUT
                 returns phi the absorption layer penetration factor
     
             SUBROUTINES
                None

        This routine is based on PHIFUN() in REC533.
     */

    let phi = if t < 0.0 {
        0.0
    } else if t <= 1.0 {
        let x = (t - 0.475) / 0.475;
        let mut phi = (((((-0.093 * x + 0.04) * x + 0.127) * x - 0.027) * x + 0.044) * x + 0.159) * x + 0.225;
        phi = f64::min(phi, 0.53);
        phi
    } else if t <= 2.2 {
        let x = (t - 1.65) / 0.55;
        let mut phi = (((((0.043 * x - 0.07) * x - 0.027) * x + 0.034) * x + 0.054) * x - 0.049) * x + 0.375;
        phi = f64::min(phi, 0.53);
        phi
    } else if t <= 10.0 {
        0.34 + (10.0 - t) * 0.02 / 7.8
    } else {
        0.34
    };

    // Multiply by the scaling factor.
    phi / 0.34
}

/// Finds the value of Lh from Table 2 ITU-R P.533-12 "Values of Lh giving auroral and other signal losses"
fn find_lh(cp: &ControlPt, dh: f64, hour: i32, month: i32) -> f64 {
    /*	
     *	FindLh() - Finds the value of Lh from Table 2 ITU-R P.533-12 "Values of Lh giving auroral and other signal losses".
     *
     *		INPUT
     *			struct ControlPt CP
     *			double dh - hop distance
     *			int hour - hour index
     *			int month - month index
     *
     *		OUTPUT
     *			returns the value of Lh
     *	
     */

    /*	
     * The upper three blocks of the array 
     *	Lh[Transmission range][Season][Mid-path local time][Geomagnetic Latitude] is for 
     *	transmission ranges less than or equal to 2500 km.
     *	The lower three blocks are for transmission ranges greater than 2500 km
     *	Each block of [8][8] is for one of three seasons: Winter, Equinox and Summer
     *	in each block eight columns are for the mid-path local time, t, in 3-hour increments
     *				  Eight rows are for the geomagnetic latitude, Gn
     */

    // a) Transmission ranges less than or equal to 2500 km
    let lh: [[[[f64; 8]; 8]; 3]; 3] = [
        [
            [
                // a) Transmission ranges less than or equal to 2500 km
                // Winter
                [2.0, 6.6, 6.2, 1.5, 0.5, 1.4, 1.5, 1.0],
                [3.4, 8.3, 8.6, 0.9, 0.5, 2.5, 3.0, 3.0],
                [6.2, 15.6, 12.8, 2.3, 1.5, 4.6, 7.0, 5.0],
                [7.0, 16.0, 14.0, 3.6, 2.0, 6.8, 9.8, 6.6],
                [2.0, 4.5, 6.6, 1.4, 0.8, 2.7, 3.0, 2.0],
                [1.3, 1.0, 3.2, 0.3, 0.4, 1.8, 2.3, 0.9],
                [0.9, 0.6, 2.2, 0.2, 0.2, 1.2, 1.5, 0.6],
                [0.4, 0.3, 1.1, 0.1, 0.1, 0.6, 0.7, 0.3]
            ],
            [
                // Equinox
                [1.4, 2.5, 7.4, 3.8, 1.0, 2.4, 2.4, 3.3],
                [3.3, 11.0, 11.6, 5.1, 2.6, 4.0, 6.0, 7.0],
                [6.5, 12.0, 21.4, 8.5, 4.8, 6.0, 10.0, 13.7],
                [6.7, 11.2, 17.0, 9.0, 7.2, 9.0, 10.9, 15.0],
                [2.4, 4.4, 7.5, 5.0, 2.6, 4.8, 5.5, 6.1],
                [1.7, 2.0, 5.0, 3.0, 2.2, 4.0, 3.0, 4.0],
                [1.1, 1.3, 3.3, 2.0, 1.4, 2.6, 2.0, 2.6],
                [0.5, 0.6, 1.6, 1.0, 0.7, 1.3, 1.0, 1.3]
            ],
            [
                //Summer
                [2.2, 2.7, 1.2, 2.3, 2.2, 3.8, 4.2, 3.8],
                [2.4, 3.0, 2.8, 3.0, 2.7, 4.2, 4.8, 4.5],
                [4.9, 4.2, 6.2, 4.5, 3.8, 5.4, 7.7, 7.2],
                [6.5, 4.8, 9.0, 6.0, 4.8, 9.1, 9.5, 8.9],
                [3.2, 2.7, 4.0, 3.0, 3.0, 6.5, 6.7, 5.0],
                [2.5, 1.8, 2.4, 2.3, 2.6, 5.0, 4.6, 4.0],
                [1.6, 1.2, 1.6, 1.5, 1.7, 3.3, 3.1, 2.6],
                [0.8, 0.6, 0.8, 0.7, 0.8, 1.6, 1.5, 1.3]
            ]
        ],
        [
            [
                // b) Transmission ranges greater than 2500 km
                // Winter
                [1.5, 2.7, 2.5, 0.8, 0.0, 0.9, 0.8, 1.6],
                [2.5, 4.5, 4.3, 0.8, 0.3, 1.6, 2.0, 4.8],
                [5.5, 5.0, 7.0, 1.9, 0.5, 3.0, 4.5, 9.6],
                [5.3, 7.0, 5.9, 2.0, 0.7, 4.0, 4.5, 10.0],
                [1.6, 2.4, 2.7, 0.6, 0.4, 1.7, 1.8, 3.5],
                [0.9, 1.0, 1.3, 0.1, 0.1, 1.0, 1.5, 1.4],
                [0.6, 0.6, 0.8, 0.1, 0.1, 0.6, 1.0, 0.5],
                [0.3, 0.3, 0.4, 0.0, 0.0, 0.3, 0.5, 0.4]
            ],
            [
                // Equinox
                [1.0, 1.2, 2.7, 3.0, 0.6, 2.0, 2.3, 1.6],
                [1.8, 2.9, 4.1, 5.7, 1.5, 3.2, 5.6, 3.6],
                [3.7, 5.6, 7.7, 8.1, 3.5, 5.0, 9.5, 7.3],
                [3.9, 5.2, 7.6, 9.0, 5.0, 7.5, 10.0, 7.9],
                [1.4, 2.0, 3.2, 3.8, 1.8, 4.0, 5.4, 3.4],
                [0.9, 0.9, 1.8, 2.0, 1.3, 3.1, 2.7, 2.0],
                [0.6, 0.6, 1.2, 1.3, 0.8, 2.0, 1.8, 1.3],
                [0.3, 0.3, 0.6, 0.6, 0.4, 1.0, 0.9, 0.6]
            ],
            [
                //Summer
                [1.9, 3.8, 2.2, 1.1, 2.1, 1.2, 2.3, 2.4],
                [1.9, 4.6, 2.9, 1.3, 2.2, 1.3, 2.8, 2.7],
                [4.4, 6.3, 5.9, 1.9, 3.3, 1.7, 4.4, 4.5],
                [5.5, 8.5, 7.6, 2.6, 4.2, 3.2, 5.5, 5.7],
                [2.8, 3.8, 3.7, 1.4, 2.7, 1.6, 4.5, 3.2],
                [2.2, 2.4, 2.2, 1.0, 2.2, 1.2, 4.4, 2.5],
                [1.4, 1.6, 1.4, 0.6, 1.4, 0.8, 2.9, 1.6],
                [0.7, 0.8, 0.7, 0.3, 0.7, 0.4, 1.4, 0.8]
            ]
        ],
        [
            [
                // b) Transmission ranges greater than 2500 km
                // Winter
                [1.5, 2.7, 2.5, 0.8, 0.0, 0.9, 0.8, 1.6],
                [2.5, 4.5, 4.3, 0.8, 0.3, 1.6, 2.0, 4.8],
                [5.5, 5.0, 7.0, 1.9, 0.5, 3.0, 4.5, 9.6],
                [5.3, 7.0, 5.9, 2.0, 0.7, 4.0, 4.5, 10.0],
                [1.6, 2.4, 2.7, 0.6, 0.4, 1.7, 1.8, 3.5],
                [0.9, 1.0, 1.3, 0.1, 0.1, 1.0, 1.5, 1.4],
                [0.6, 0.6, 0.8, 0.1, 0.1, 0.6, 1.0, 0.5],
                [0.3, 0.3, 0.4, 0.0, 0.0, 0.3, 0.5, 0.4]
            ],
            [
                // Equinox
                [1.0, 1.2, 2.7, 3.0, 0.6, 2.0, 2.3, 1.6],
                [1.8, 2.9, 4.1, 5.7, 1.5, 3.2, 5.6, 3.6],
                [3.7, 5.6, 7.7, 8.1, 3.5, 5.0, 9.5, 7.3],
                [3.9, 5.2, 7.6, 9.0, 5.0, 7.5, 10.0, 7.9],
                [1.4, 2.0, 3.2, 3.8, 1.8, 4.0, 5.4, 3.4],
                [0.9, 0.9, 1.8, 2.0, 1.3, 3.1, 2.7, 2.0],
                [0.6, 0.6, 1.2, 1.3, 0.8, 2.0, 1.8, 1.3],
                [0.3, 0.3, 0.6, 0.6, 0.4, 1.0, 0.9, 0.6]
            ],
            [
                //Summer
                [1.9, 3.8, 2.2, 1.1, 2.1, 1.2, 2.3, 2.4],
                [1.9, 4.6, 2.9, 1.3, 2.2, 1.3, 2.8, 2.7],
                [4.4, 6.3, 5.9, 1.9, 3.3, 1.7, 4.4, 4.5],
                [5.5, 8.5, 7.6, 2.6, 4.2, 3.2, 5.5, 5.7],
                [2.8, 3.8, 3.7, 1.4, 2.7, 1.6, 4.5, 3.2],
                [2.2, 2.4, 2.2, 1.0, 2.2, 1.2, 4.4, 2.5],
                [1.4, 1.6, 1.4, 0.6, 1.4, 0.8, 2.9, 1.6],
                [0.7, 0.8, 0.7, 0.3, 0.7, 0.4, 1.4, 0.8]
            ]
        ]
    ];

    // Initialize the geomagnetic coordinates.
    let mut gn = Location { lat: 0.0, lng: 0.0 };

    // Find the geomagnetic coordinates for location of the control point.
    crate::geometry::geomagnetic_coords(cp.l, &mut gn);

    // Determine the season index for the Lh array.
    let season = what_season_for_lh(cp.l, month); // season index

    // Lh[Transmission range][season][geomagnetic latitude][mid-path local time]
    // Determine the indices
    // Transmit range index
    let txrange = if dh <= 2500.0 { 0 } else { 1 };

    // Midpoint local time index
    let mplt = match hour {
        1..=3 => 0,
        4..=6 => 1,
        7..=9 => 2,
        10..=12 => 3,
        13..=15 => 4,
        16..=18 => 5,
        19..=21 => 6,
        _ => 7 // 22..=23 or 0
    };

    let gn_lat = f64::abs(gn.lat);

    let gmlat = if gn_lat >= 77.5 * D2R {
        0
    } else if gn_lat >= 72.5 * D2R && gn_lat < 77.5 * D2R {
        1
    } else if gn_lat >= 67.5 * D2R && gn_lat < 72.5 * D2R {
        2
    } else if gn_lat >= 62.5 * D2R && gn_lat < 67.5 * D2R {
        3
    } else if gn_lat >= 57.5 * D2R && gn_lat < 62.5 * D2R {
        4
    } else if gn_lat >= 52.5 * D2R && gn_lat < 57.5 * D2R {
        5
    } else if gn_lat >= 47.5 * D2R && gn_lat < 52.5 * D2R {
        6
    } else if gn_lat >= 42.5 * D2R && gn_lat < 47.5 * D2R {
        7
    } else {
        // gn_lat < 42.5 * D2R
        return 0.0; // Nothing more to do return 0.0 as Lh
    };

    // Lh[transmission range][season][geomagnetic latitude][mid-path local time]
    lh[txrange][season][gmlat][mplt]
}

/// Determines the index into the Lh array dependent on the month and latitude
fn what_season_for_lh(l: Location, month: i32) -> usize {
    /*
     *	WhatSeasonforLh() - Determines the index into the Lh array dependent on the month and latitude.
     *
     *		INPUT
     *			struct Location L - location of interest
     *			int month - month index
     *
     *		OUTPUT
     *			returns the season
     */

    let season = if l.lat >= 0.0 {
        // Northern hemisphere and the equator
        match month {
            DEC | JAN | FEB => WINTER,
            MAR | APR | MAY | SEP | OCT | NOV => EQUINOX,
            JUN | JUL | AUG => SUMMER,
            _ => EQUINOX
        }
    } else {
        // Southern hemisphere 
        match month {
            JUN | JUL | AUG => WINTER,
            MAR | APR | MAY | SEP | OCT | NOV => EQUINOX,
            DEC | JAN | FEB => SUMMER,
            _ => EQUINOX
        }
    };

    season
}

/// Finds the antenna gain at the desired elevation, delta
pub fn antenna_gain(path: &PathData, ant: &Antenna, delta: f64, direction: i32) -> f64 {
    /*
        AntennaGain() - Finds the antenna gain at the desired elevation, delta

            INPUT
                struct PathData path
                struct Antenna Ant
                double delta

            OUTPUT
                returns the interpolated antenna gain at the desired elevation, delta

            SUBROUTINES
                Bearing()
                BilinearInterpolation()
    */

    // The structure Antenna Ant is used to tell the subroutine which antenna to calculate.

    let mut freq_index = 0usize;
    /* If we have pattern data for multiple frequencies, find the index of the
     * frequency closest to the path.frequency.
     */
    if ant.freqn > 1 {
        let mut min_freq_delta = f64::MAX;
        for i in 0..ant.freqn {
            let freq_delta = f64::abs(ant.freqs[i as usize] - path.frequency);
            if freq_delta < min_freq_delta {
                min_freq_delta = freq_delta;
                freq_index = i as usize;
            }
        }
    }

    // The elevations and azimuths have to be in degrees because the antenna pattern is indexed in degrees.

    // delta is in radians convert to degrees
    let delta = delta * R2D;

    // Bearing from transmitter to receiver
    let b = match direction {
        // Determine the bearing
        // From the tx to rx.
        TX_TO_RX => crate::geometry::bearing(path.l_tx, path.l_rx) * R2D,
        RX_TO_TX => crate::geometry::bearing(path.l_rx, path.l_tx) * R2D,
        _ => 0.0
    };

    // Now determine the gain at the elevation, delta
    // Find the indices to determine the neighbors for the gain interpolation.
    let delta_u = f64::ceil(delta) as usize; // Upper elevation indices
    let delta_l = f64::floor(delta) as usize; // Lower elevation indices

    // The bearing might wrap around.
    let br = (f64::ceil(b) as i32 % 360) as usize; // Left bearing indices
    let bl = (f64::floor(b) as i32 % 360) as usize; // Right bearing indices

    // Identify the neighbors.
    let ll = ant.pattern[freq_index][bl][delta_l];
    let lr = ant.pattern[freq_index][br][delta_l];
    let ul = ant.pattern[freq_index][bl][delta_u];
    let ur = ant.pattern[freq_index][br][delta_u];

    // Determine the fractional column and row.
    // The distance between indices is fixed at 1 degree.
    let r = delta - f64::floor(delta); // The fractional part of the row  (elevation)
    let c = b - f64::floor(b); // The fractional part of the column (azimuth)
    let g = crate::calculate_cp_parameters::bilinear_interpolation(ll, lr, ul, ur, r, c); // Interpolated gain

    g
}
