use crate::constants::*;
use crate::path_data::{PathData, ControlPt, Location};

#[derive(Debug, Clone)]
pub struct Neighbor {
    pub l: Location,                // the physical location of the point
    pub fof2: [f64; 2],            // the values for the ionospheric parameters at the point (SSN=0, SSN=100)
    pub m3kf2: [f64; 2],           // the values for the ionospheric parameters at the point (SSN=0, SSN=100)
    pub j: i32,                     // the j (lng) index into the gridpoint map
    pub k: i32,                     // the k (lat) index into the gridpoint map
}

impl Default for Neighbor {
    fn default() -> Self {
        Neighbor {
            l: Location { lat: 0.0, lng: 0.0 },
            fof2: [0.0; 2],
            m3kf2: [0.0; 2],
            j: 0,
            k: 0,
        }
    }
}

/// Calculate control point parameters
/// Finds the ionospheric parameters foF2 and M3kF2 for high and low SSN, 
/// gyrofrequency and magnetic dip at 100 and 300 km, and solar parameters
pub fn calculate_cp_parameters(path: &PathData, cp: &mut ControlPt) {
    calculate_cp_parameters_impl(cp, &path.fof2, &path.m3kf2, path.hour, path.ssn, path.month);
}

/// Calculate control point parameters implementation with individual parameters
/// This allows calling without borrowing issues from the main path struct
pub fn calculate_cp_parameters_impl(
    cp: &mut ControlPt, 
    fof2_data: &Vec<Vec<Vec<Vec<f32>>>>, 
    m3kf2_data: &Vec<Vec<Vec<Vec<f32>>>>, 
    hour: i32, 
    ssn: i32, 
    month: i32
) {
    /*
        CalculateCPParameters() finds the ionosheric parameters foF2 and M3kF2 for the high (SSN = 100) and low (SSN = 0), the gyrofrequency and magnetic dip at 100 and 300 km 
            and the solar parameters for the control point. The ionospheric parameters are determined by the bi-linear interpolation method in P.1144-3 (2001) 
            and linear interpolation by the desired SSN (Found in path.SSN - the user selected SSN) the routine determines the foF2 and M3F2 at the point of interest. 
            This routine further determines foE at the point of interest. 

            INPUT
                struct PathData *path
                struct ControlPt *here - This is a pointer to the control point of interest.

            OUTPUT
                There are four subroutines in this program that calculate the parameters for 
                the control point. Which subroutine calulates which parameter is summarized below

                via IonosphericParameters()
                    here.foF2 - Critical frequency of the F2 layer
                    here.M3kF2 - Critical frequency of the F2 layer for 3000 km
                via FindFoE()
                    here.foE - Critical frequency fo the E layer
                    here.ltime - Local time
                via magfit()
                    here.dip[2] - Magnetic dip calculated at 100 and 300 km
                    here.fH[2] - Gyrofreqency calculated at 100 and 300 km
                via SolarParameters() 
                    here.Sun.ha - Hour angle (radians)
                    here.Sun.sha -  Sunrise/Sunset hour angle (radians)
                    here.Sun.sza - Solar zenith angle (radians)
                    here.Sun.decl - Solar declination (radians)
                    here.Sun.eot- Equation of time (minutes)
                    here.Sun.lsr - local sunrise (hours)
                    here.Sun.lsn - local solar noon (hours)
                    here.Sun.lss - local sunset (hours)

            SUBROUTINES
                IonosphericParameters()
                SolarParameters()
                FindfoE()
                magfit()
    */

    /*
     * Find the ionospheric parameters foF2 and M3kF2 at the control point here.
     * If here is not on a grid point then use bilinear interpolation.
     */
    ionospheric_parameters(cp, fof2_data, m3kf2_data, hour, ssn);

    /*
     * Calculate the magnetic dip and gyrofrequency for this location
     */
    magfit(cp, 100.0);
    magfit(cp, 300.0);

    /*
     * Calculate the solar parameters. 
     * These parametres are used in the MedianSkywaveFieldStrengthLong() calculation. Because this 
     * routine determines the control point parameters for all methods, find the solar parameters now
     * before entering the conditional loop for the foE calculation.
     */
    // Find the solar parameters for the control point.
    solar_parameters(cp, month, hour as f64);

    /*
     * Calculate foE by the method outlined in P.1239-2. 
     */
    find_foe(cp, month, hour as f64, ssn);

    /* 
     * At each control point the gyrofrequency and magnetic dip must also be calculated.
     * The calculation is done at two heights: 
     *		height = 300 km is for the determination of MUF and the long path (> 9000 km).
     *      height = 100 km is used in the determination of absorption on the int paths (< 9000 km).
     */

    magfit(cp, 100.0);
    magfit(cp, 300.0);
}

/// Finds foF2 and M(3000)F2 by Bilinear Interpolation
fn ionospheric_parameters(
    cp: &mut ControlPt,
    fof2: &Vec<Vec<Vec<Vec<f32>>>>,
    m3kf2: &Vec<Vec<Vec<Vec<f32>>>>,
    hour: i32,
    mut ssn: i32,
) {
    /*     
      IonosphericParameters() - Finds foF2 and M(3000)F2 by Bilinear Interpolation.
     
          INPUTS
             struct ControlPt *here,
              double ****foF2
             double ****M3kF2
             int hour
             int SSN
     
         OUTPUT
             here->foF2
             here->M3KF2

        SUBROUTINES
            BilinearInterpolation()

         I am indebted to Peter Suessman for the use of his program, iongrid ver 1.80, which was used extensively
         to verify the method in this routine. 
    */

    let mut ul = Neighbor::default(); // upper left
    let mut ur = Neighbor::default(); // upper right
    let mut lr = Neighbor::default(); // lower right
    let mut ll = Neighbor::default(); // lower left

    let mut ned = Neighbor::default();

    // For the gridpoint maps at an increment of 1.5 degrees in lat and long
    const ZERO_LAT: i32 = 60; // The max latitude is 121 
    const ZERO_LNG: i32 = 120; // The max longitude is 241

    // This routine is dependent on the 1.5 degree increment
    const INC: f64 = 1.5 * D2R; // Increment for the gridpoint maps (1.5 * pi) / 180 = 0.0261799388

    const LNG: i32 = 241; // 241 Longitudes at 1.5 degree increments
    const LAT: i32 = 121; // 121 latitudes at 1.5 degree increments

    /*
     * Find the neighborhood around the point of interest
     * The neighbors are determined differently for each quadrant so that the neighbors are in the correct order relative to the 
     * array indices. This is done to simplify the interpolation code. There are no quadrants in the foF2 and M3kF2 data array. In the 
     * foF2 and M3kF2 array the 0,0 point is the southwest corner and the indices increase north and east. 
     * First determine the quadrant
     */
    if cp.l.lat >= 0.0 {
        if cp.l.lng >= 0.0 {
            // NE quad
            ll.k = ZERO_LAT + (cp.l.lat / INC) as i32; // North is positive
            ll.j = ZERO_LNG + (cp.l.lng / INC) as i32; // East is positive
            lr.k = ll.k;
            lr.j = ll.j + 1;
            ur.k = ll.k + 1;
            ur.j = ll.j + 1;
            ul.k = ll.k + 1;
            ul.j = ll.j;
            // Check the rollover conditions - corners and edges
            if ll.j != LNG - 1 {
                // not on E edge
                if ll.k != LAT - 1 {
                    // not on N edge
                    // don't do anything
                } else {
                    // N edge
                    // fix the N rollover - no rollover N to S
                    ur.k = ll.k;
                    ul.k = ll.k;
                }
            } else {
                // E edge
                if ll.k != LAT - 1 {
                    // not on N edge
                    // fix the E rollover - rollover E to W
                    lr.j = 0;
                    ur.j = 0;
                } else {
                    // NE corner
                    // Collapse to 1D
                    lr.k = ll.k;
                    lr.j = 0;
                    ur.k = ll.k;
                    ur.j = 0;
                    ul.k = ll.k;
                    ul.j = ll.j;
                }
            }
        } else {
            // NW quad
            lr.k = ZERO_LAT + (cp.l.lat / INC) as i32; // North is positive
            lr.j = ZERO_LNG + (cp.l.lng / INC) as i32; // West is negative	
            ll.k = lr.k;
            ll.j = lr.j - 1;
            ul.k = lr.k + 1;
            ul.j = lr.j - 1;
            ur.k = lr.k + 1;
            ur.j = lr.j;
            // Check the rollover conditions - corners and edges
            if lr.j != 0 {
                // not on W edge
                if lr.k != LAT - 1 {
                    // not on N edge
                    // don't do anything
                } else {
                    // N edge
                    // fix the N rollover - no rollover N to S
                    ur.k = lr.k;
                    ul.k = lr.k;
                }
            } else {
                // W edge
                if lr.k != LAT - 1 {
                    // not on N edge
                    // fix the W rollover - rollover W to E
                    ll.j = LNG - 1;
                    ul.j = LNG - 1;
                } else {
                    // NW corner
                    // Collapse to 1D
                    ll.k = lr.k;
                    ll.j = LNG - 1;
                    ur.k = lr.k;
                    ur.j = lr.j;
                    ul.k = lr.k;
                    ul.j = LNG - 1;
                }
            }
        }
    } else {
        if cp.l.lng >= 0.0 {
            // SE quad
            ul.k = ZERO_LAT + (cp.l.lat / INC) as i32; // South is negative
            ul.j = ZERO_LNG + (cp.l.lng / INC) as i32; // East is positive
            ur.k = ul.k;
            ur.j = ul.j + 1;
            ll.k = ul.k - 1;
            ll.j = ul.j;
            lr.k = ul.k - 1;
            lr.j = ul.j + 1;
            // Check the rollover conditions - corners and edges
            if ul.j != LNG - 1 {
                // not on E edge
                if ul.k != 0 {
                    // not on SE corner
                    // don't do anything
                } else {
                    // S edge
                    // fix the S rollover - no rollover S to N
                    ll.k = ul.k;
                    lr.k = ul.k;
                }
            } else {
                // E edge
                if ul.k != 0 {
                    // not on SE corner
                    // fix the E rollover - rollover E to W
                    lr.j = 0;
                    ur.j = 0;
                } else {
                    // SE corner
                    // Collapse to 1D
                    lr.k = ul.k;
                    lr.j = 0;
                    ur.k = ul.k;
                    ur.j = 0;
                    ll.k = ul.k;
                    ll.j = ul.j;
                }
            }
        } else {
            // SW quad
            ur.k = ZERO_LAT + (cp.l.lat / INC) as i32; // South is negative
            ur.j = ZERO_LNG + (cp.l.lng / INC) as i32; // West is negative
            ul.k = ur.k;
            ul.j = ur.j - 1;
            ll.k = ur.k - 1;
            ll.j = ur.j - 1;
            lr.k = ur.k - 1;
            lr.j = ur.j;
            // Check the rollover conditions - corners and edges
            if ur.j != 0 {
                // not on W edge
                if ur.k != 0 {
                    // not on SW corner
                    // don't do anything
                } else {
                    // S edge
                    // fix the N rollover - no rollover N to S
                    lr.k = ur.k;
                    ll.k = ur.k;
                }
            } else {
                // W edge
                if ur.k != 0 {
                    // not on SW corner
                    // fix the W rollover - rollover W to E
                    ll.j = LNG - 1;
                    ul.j = LNG - 1;
                } else {
                    // SW corner
                    // Collapse to 1D
                    lr.k = ur.k;
                    lr.j = ur.j;
                    ll.k = ur.k;
                    ll.j = LNG - 1;
                    ul.k = ur.k;
                    ul.j = LNG - 1;
                }
            }
        }
    }

    // Determine the lat and lng in degrees
    ul.l.lat = ul.k as f64 * INC - PI / 2.0;
    ul.l.lng = ul.j as f64 * INC - PI;
    ll.l.lat = ll.k as f64 * INC - PI / 2.0;
    ll.l.lng = ll.j as f64 * INC - PI;
    ur.l.lat = ur.k as f64 * INC - PI / 2.0;
    ur.l.lng = ur.j as f64 * INC - PI;
    lr.l.lat = lr.k as f64 * INC - PI / 2.0;
    lr.l.lng = lr.j as f64 * INC - PI;

    /*
     * At this point you have the neighborhood around the point now you can populate the foF2 and
     * M3kF2 for each of the adjacent points
     */
    let n = hour; // Temp hours index

    for m in 0..2 {
        // SSN
        // Upper Left
        ul.fof2[m] = fof2[n as usize][ul.j as usize][ul.k as usize][m] as f64;
        ul.m3kf2[m] = m3kf2[n as usize][ul.j as usize][ul.k as usize][m] as f64;
        // Upper Right
        ur.fof2[m] = fof2[n as usize][ur.j as usize][ur.k as usize][m] as f64;
        ur.m3kf2[m] = m3kf2[n as usize][ur.j as usize][ur.k as usize][m] as f64;
        // Lower Left
        ll.fof2[m] = fof2[n as usize][ll.j as usize][ll.k as usize][m] as f64;
        ll.m3kf2[m] = m3kf2[n as usize][ll.j as usize][ll.k as usize][m] as f64;
        // Lower Right
        lr.fof2[m] = fof2[n as usize][lr.j as usize][lr.k as usize][m] as f64;
        lr.m3kf2[m] = m3kf2[n as usize][lr.j as usize][lr.k as usize][m] as f64;
    }

    // Now you are ready to interpolate the value at the point of interest
    // determine the fractional "column" j and fractional "row" k for the bilinear interpolation calculation
    let frack = f64::abs(cp.l.lat / INC) - (f64::abs(cp.l.lat / INC) as i32) as f64; // Fractional row distance
    // Fractional "column" j and fractional "row" k
    let fracj = f64::abs(cp.l.lng / INC) - (f64::abs(cp.l.lng / INC) as i32) as f64; // Fractional column distance
    for m in 0..2 {
        ned.fof2[m] = bilinear_interpolation(ll.fof2[m], lr.fof2[m], ul.fof2[m], ur.fof2[m], frack, fracj);
        ned.m3kf2[m] = bilinear_interpolation(ll.m3kf2[m], lr.m3kf2[m], ul.m3kf2[m], ur.m3kf2[m], frack, fracj);
    }

    // End of calculation for foF2 and M3kF2

    /*
     * Now interpolate by the SSN. Note the SSN maximum has been restricted to a maximum of 160 ITU-R P.533-12.
     * "For most purposes it is adequate to assume a linear relationship with R12 for both foF2 and M(3000)F2." 
     * ITU-R P.1239-2 (10-2009)
     * Note the index on foF2 and M3kF2 in the neighbor structure is for the SSN = 0 (index = 0) and SSN = 100 (index = 1)
     */
    ssn = i32::min(ssn, MAX_SSN);
    cp.fof2 = (ned.fof2[1] * ssn as f64 + ned.fof2[0] * (100.0 - ssn as f64)) / 100.0;
    cp.m3kf2 = (ned.m3kf2[1] * ssn as f64 + ned.m3kf2[0] * (100.0 - ssn as f64)) / 100.0;

    // End of calculation for foF2 and M3kF2
}

/// Calculate solar parameters for the control point
pub fn solar_parameters(cp: &mut ControlPt, month: i32, hour: f64) {
    /*     
         SolarParameters() - Calculate the solar parameters at the control point for the given 
             time and month.

             INPUT
                 struct ControlPt *here - The control point of interest
                 int month - Month index
                 double hour - Decimal hours

             OUTPUT
                 here->Sun.ha - Hour angle (radians)
                 here->Sun.sha -  Sunrise/Sunset hour angle (radians)
                 here->Sun.sza - Solar zenith angle (radians)
                  here->Sun.decl - Solar declination (radians)
                  here->Sun.eot- Equation of time (minutes)
                 here->Sun.lsr - local sunrise (hours)
                 here->Sun.lsn - local solar noon (hours)
                 here->Sun.lss - local sunset (hours)

            SUBROUTINES
                None

        Thanks to the following references
        See www.analemma.com/Pages/framesPage.html
        See holbert.faculty.asu.edu/eee463/SolarCalcs.pdf
        Although W is + and E is - and the time zones are also reversed 
        See www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
    */

    let a = 0.98565327; // Average angle per day
    let b = 3.98891967; // Minutes per degree of Earth's rotation
    let s = f64::sin(23.45 * D2R); // Earth's tilt sine
    let c = f64::cos(23.45 * D2R); // Earth's tilt cosine
    let v = 78.746118 * D2R; // Value of nu on March 21st

    // The day of the year (doty) array allows us to determine the day count of the day of interest
    let doty = [0, 31, 59, 90, 120, 152, 181, 212, 243, 273, 304, 334];

    // Determine the local time, hours, minutes, seconds and time zone
    let ltime = hour + (cp.l.lng / (15.0 * D2R)) as i32 as f64; // Local time 
    // Local time 
    let tzone = (cp.l.lng / (15.0 * D2R)) as i32 as f64; // hours
    // Time zone
    // At present this code only works for the 15th day of the month
    // If this changes a day field should be added to the path structure
    // and passed into this routine
    let day = 15;

    let d = doty[month as usize] as f64 + day as f64 + hour / 24.0;

    // Calculate the Equation-of-Time
    // First find the time due to the elliptic orbit of the Earth
    // The average day is 360 degrees / 365.25 days a year assuming a circular orbit equals 0.985653
    // Assume that the perihelion ( The Earth is closest to the sun ) is on January 2nd.
    // So in D days of the year the earth moves through lambda degrees
    let lambda = a * D2R * (d - 2.0);

    // Determine the arc length due to an elliptical orbit
    // ( 360 degrees / PI ) * 0.016713 the shape factor of the elliptic equals 1.915169
    let nu = lambda + 1.915169 * D2R * f64::sin(lambda);

    // Find the angles associated with the tilt of the Earth
    // epsilon is the mean sun angle of the Earth after N - 80 days
    let mut epsilon = a * D2R * (d - 80.0);

    // epsilon is +- PI/2
    if epsilon >= 270.0 * D2R {
        epsilon -= 2.0 * PI;
    } else if epsilon >= 90.0 * D2R {
        epsilon -= PI;
    }

    // The angle of the true sun is beta
    let beta = f64::atan(c * f64::tan(epsilon));

    // Equation of Time = tilt effect + elliptic effect
    // Where 0.398892 is the minutes per degree of Earth's rotation 
    // 1440 minutes per day /361 degrees per day 
    cp.sun.eot = b * (epsilon - beta + (lambda - nu)) * R2D;

    // Solar declination in radians
    cp.sun.decl = f64::asin(s * f64::sin(f64::sin(a * (d - 2.0) * D2R) * 0.016713 + a * (d - 2.0) * D2R - v));

    // Find the hour angle which can be found from the solar time corrected for the local longitude and the eot
    let toffset = (cp.l.lng / (15.0 * D2R) - tzone) * 60.0 + cp.sun.eot; // minutes

    let tst = ltime * 60.0 + toffset; // Apparent/Local/True solar time in minutes
    // True solar time
    cp.sun.ha = (tst / 4.0 - 180.0) * D2R; // radians

    // Hour angle at sunrise and sunset in radians
    cp.sun.sha = f64::acos(f64::cos(90.833 * D2R) / (f64::cos(cp.l.lat) * f64::cos(cp.sun.decl)) - f64::tan(cp.l.lat) * f64::tan(cp.sun.decl));

    // The cosine of the solar zenith angle can be found
    let mut cosphi = f64::sin(cp.l.lat) * f64::sin(cp.sun.decl) + f64::cos(cp.l.lat) * f64::cos(cp.sun.decl) * f64::cos(cp.sun.ha); // cosine of the solar zenith angle

    /* (watch out for the roundoff errors) */
    if f64::abs(cosphi) > 1.0 {
        if cosphi >= 0.0 {
            cosphi = 1.0;
        } else {
            cosphi = -1.0;
        }
    }

    cp.sun.sza = f64::acos(cosphi); // Solar zenith angle which will be positive even in the for southern latitudes

    // Switch the sign of the longitude for the time calculation
    // Local Sunrise relative to UTC in fractional hours
    cp.sun.lsr = (720.0 + (-cp.l.lng - cp.sun.sha) * R2D * 4.0 - cp.sun.eot) / 60.0;

    // Local Sunset relative to UTC in fractional hours
    cp.sun.lss = (720.0 + (-cp.l.lng + cp.sun.sha) * R2D * 4.0 - cp.sun.eot) / 60.0;

    // Local Solar noon relative to UTC in fractional hours
    cp.sun.lsn = (720.0 + -cp.l.lng * R2D * 4.0 - cp.sun.eot) / 60.0;

    // Roll over the times. Note: add 24 because for the fmod(x, 24) x might be negative 
    cp.sun.lsr = (cp.sun.lsr + 24.0) % 24.0;
    cp.sun.lss = (cp.sun.lss + 24.0) % 24.0;
    cp.sun.lsn = (cp.sun.lsn + 24.0) % 24.0;

    // Store the UTC time to here structure
    cp.ltime = hour;
}

/// Calculate foE for the control point
fn find_foe(cp: &mut ControlPt, month: i32, hour: f64, mut ssn: i32) {
    /* 
        FindFoE() - Determines the critical frequency of the E-layer, foE, by the method in 
        ITU-R P.1239. This routine assumes that the solar parameters have been calculated 
        for the control point before execution.

        INPUT
            struct ControlPt *here
            int month
            int hour
            int SSN

        OUTPUT
            here.foE critical frequency of the E-layer determined for the control point here

        SUBROUTINES
            None
    */

    // Temps for the calculation of foE
    let m: f64;
    let n: f64;
    let x: f64;
    let y: f64; // temps
    // End of Temporary Variables

    // Restrict the ssn to 160
    ssn = i32::min(ssn, MAX_SSN);

    /*
     * Now find the foE for the control point "here"
     * There are four parts to the calculation 
     *		A = solar activity factor
     *		B = seasonal factor
     *		C = main latitude factor
     *		D = time-of-day factor
     */
    /* 
     * Calculation for A : solar activity factor
     * First find phi sub 12 (phi) by eqn (2) in P.1239-2 (2009)
     */
    let phi = 63.7 + 0.728 * ssn as f64 + 0.00089 * f64::powf(ssn as f64, 2.0); // monthly mean 10.7 cm solar radio flux
    let a = 1.0 + 0.0094 * (phi - 66.0); // solar activity factor

    /*
     * Calculation for B : seasonal factor
     */
    if f64::abs(cp.l.lat) < 32.0 * D2R {
        m = -1.93 + 1.92 * f64::cos(cp.l.lat);
    } else {
        // fabs(here.L.lat) >= 32*D2R
        m = 0.11 - 0.49 * f64::cos(cp.l.lat);
    }

    if f64::abs(cp.l.lat - cp.sun.decl) < 80.0 * D2R {
        n = cp.l.lat - cp.sun.decl;
    } else {
        n = 80.0 * D2R;
    }

    let b = f64::powf(f64::cos(n), m); // seasonal factor

    /* 
     * Calculation for C : main latitude factor
     */
    if f64::abs(cp.l.lat) < 32.0 * D2R {
        x = 23.0;
        y = 116.0;
    } else {
        // fabs(here.L.lat) >= 32*D2R
        x = 92.0;
        y = 35.0;
    }

    let c = x + y * f64::cos(cp.l.lat); // main latitude factor
    /*
     * Calculation of D : time-of-day factor
     */
    // In each case in determining D from the solar zenith angle (here.sza) the p exponent is determined the same
    let p = if f64::abs(cp.l.lat) <= 12.0 * D2R { 1.31 } else { 1.2 }; // coefficient

    // Now calculate D conditional on the solar zenith angle (here.sza)
    let d = if cp.sun.sza <= 73.0 * D2R {
        f64::powf(f64::cos(cp.sun.sza), p)
    } else if cp.sun.sza > 73.0 * D2R && cp.sun.sza < PI / 2.0 {
        // Twilight is 90 degrees
        let dsza = 6.27e-13 * f64::powf(cp.sun.sza * R2D - 50.0, 8.0) * D2R; // delta solar zenith angle
        f64::powf(f64::cos(cp.sun.sza - dsza), p)
    } else {
        // (here.sza >= 90.0*D2R )
        // In this case local sunset and sunrise must be known.
        // Find h the number of hours after sunset
        let hour_adj = ((hour as i32 + 1) % 24) as f64; // Adjust time and roll over

        let h; // coefficient

        if cp.sun.lss >= cp.sun.lsr && hour_adj >= cp.sun.lss && hour_adj >= cp.sun.lsr {
            h = hour_adj - cp.sun.lss;
        } else if cp.sun.lss < cp.sun.lsr && hour_adj >= cp.sun.lss && hour_adj < cp.sun.lsr {
            h = hour_adj - cp.sun.lss;
        } else if cp.sun.lss >= cp.sun.lsr && hour_adj < cp.sun.lss && hour_adj < cp.sun.lsr {
            h = 24.0 - cp.sun.lss + hour_adj;
        } else {
            h = 0.0;
        }

        // If it is night determine if here is in a polar region and during a period of polar winter
        // The Norwegian territory of Svalbad is known to experiences a civil polar night lasting 
        // from about 11 November until 30 January. 
        // Civil Polar night is when the sun is 6 degrees below the horizon. 
        // Because the arctic circle is at 66.5622 degrees, civil twilight would be 72.5622 degrees. Although civil twilight is chosen here,
        // the precise angle of ionospheric twilight is open to debate. Half the month of November to the 1st of February experiences polar 
        // night. This program works on median months, so polar night will be defined as ...
        if (cp.l.lat > 72.5622 * D2R &&
            (month == NOV || month == DEC || month == JAN))
            ||
            (cp.l.lat < -72.5622 * D2R && (month == MAY || month == JUN ||
                                           month == JUL)) {
            // Northern hemisphere polar winter || Southern hemisphere polar winter
            f64::powf(0.072, p) * f64::exp(25.2 - 0.28 * cp.sun.sza * R2D)
        } else {
            // Choose the larger of the two calculations
            f64::max(f64::powf(0.072, p) * f64::exp(-1.4 * h),
                f64::powf(0.072, p) * f64::exp(25.2 - 0.28 * cp.sun.sza * R2D))
        }
    };

    // Choose the larger of the foE calculations
    cp.foe = f64::max(f64::powf(a * b * c * d, 0.25),
        f64::powf(0.004 * f64::powf(1.0 + 0.021 * phi, 2.0), 0.25));
}

/// Calculate magnetic parameters at specified height
fn magfit(cp: &mut ControlPt, height: f64) {
    /*
     magfit() calculates the magnetic dip and the gyrofrequency
        This calculation is described in Section 2 of Rec P.1239 (1997) equations (5) thru (11).
        This subroutine is patterned after the Fortran subroutine by the same name in the ITS
        propagation package, 2006.
        The input lat and long are in radians.

         Initialized the arrays
         P = The Associated Legendre function
         DP = The derivative of P
         CT = Appears to be the Associated Legendre function coefficients as a function of m and n
         G & H = Numerical coefficients for the field model (gauss)

         INPUT
             struct ControlPt *here - Control point of interest
             double height - Height at which the calculation is made

          OUTPUT
             here.dip[hr] - Magnetic dip
             here.fH[hr] - Gyrofrequency

         SUBROUTINES
            None

     For the following input parameters produce the following outputs for the Fh, gyrofrequency, and I, magnetic dip angle.
     
     lat_d =   0.732665 radians or (41 + 58/60 + 43/3600) degrees 
     long_d =  -1.534227 radians or -(87 + 54/60 + 17/3600) degrees
     height = 1800;

     Fh = 0.733686 
     I = 1.241669

     This routine is based on MAGFIT.FOR in REC533.

    */

    let mut p: [[f64; 7]; 7] = [
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ];

    let mut dp: [[f64; 7]; 7] = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ];

    let g: [[f64; 7]; 7] = [
        [0.000000, 0.304112, 0.024035, -0.031518, -0.041794, 0.016256, -0.019523],
        [0.000000, 0.021474, -0.051253, 0.062130, -0.045298, -0.034407, -0.004853],
        [0.000000, 0.000000, -0.013381, -0.024898, -0.021795, -0.019447, 0.003212],
        [0.000000, 0.000000, 0.000000, -0.0064960, 0.007008, -0.000608, 0.021413],
        [0.000000, 0.000000, 0.000000, 0.000000, -0.002044, 0.002775, 0.001051],
        [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000697, 0.000227],
        [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.001115]
    ];

    let h: [[f64; 7]; 7] = [
        [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
        [0.000000, -0.057989, 0.033124, 0.014870, -0.011825, -0.000796, -0.005758],
        [0.000000, 0.000000, -0.001579, -0.004075, 0.010006, -0.002000, -0.008735],
        [0.000000, 0.000000, 0.000000, 0.000210, 0.000430, 0.004597, -0.003406],
        [0.000000, 0.000000, 0.000000, 0.000000, 0.001385, 0.002421, -0.000118],
        [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.001218, -0.001116],
        [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.000325]
    ];

    let ct: [[f64; 7]; 7] = [
        [0.0000000, 0.0000000, 0.33333333, 0.266666666, 0.25714286, 0.25396825, 0.25252525],
        [0.0000000, 0.0000000, 0.00000000, 0.200000000, 0.22857142, 0.23809523, 0.24242424],
        [0.0000000, 0.0000000, 0.00000000, 0.000000000, 0.14285714, 0.19047619, 0.21212121],
        [0.0000000, 0.0000000, 0.00000000, 0.000000000, 0.00000000, 0.11111111, 0.16161616],
        [0.0000000, 0.0000000, 0.00000000, 0.000000000, 0.00000000, 0.00000000, 0.09090909],
        [0.0000000, 0.0000000, 0.00000000, 0.000000000, 0.00000000, 0.00000000, 0.00000000],
        [0.0000000, 0.0000000, 0.00000000, 0.000000000, 0.00000000, 0.00000000, 0.00000000]
    ];

    let hr = match height {
        // dip and fH can only be calculated for 2 heights in this project
        100.0 => HR_100_KM,
        300.0 => HR_300_KM,
        _ => HR_100_KM
    };

    let ar = R0 / (R0 + height);

    let mut fz = 0.0;
    let mut fx = 0.0;
    let mut fy = 0.0;

    for n in 1..=6 {
        let mut sumz = 0.0;
        let mut sumx = 0.0;
        let mut sumy = 0.0;

        for m in 0..=n {
            if n == m {
                p[m][n] = f64::cos(cp.l.lat) * p[m - 1][n - 1];
                dp[m][n] = f64::cos(cp.l.lat) * dp[m - 1][n - 1] + f64::sin(cp.l.lat) * p[m - 1][n - 1];
            } else if n != 1 {
                p[m][n] = f64::sin(cp.l.lat) * p[m][n - 1] - ct[m][n] * p[m][n - 2];
                dp[m][n] = f64::sin(cp.l.lat) * dp[m][n - 1] - f64::cos(cp.l.lat) * p[m][n - 1] -
                           ct[m][n] * dp[m][n - 2];
            } else {
                p[m][n] = f64::sin(cp.l.lat) * p[m][n - 1];
                dp[m][n] = f64::sin(cp.l.lat) * dp[m][n - 1] - f64::cos(cp.l.lat) * p[m][n - 1];
            }

            sumz += p[m][n] * (g[m][n] * f64::cos(m as f64 * cp.l.lng) + h[m][n] * f64::sin(m as f64 * cp.l.lng));
            sumx += dp[m][n] * (g[m][n] * f64::cos(m as f64 * cp.l.lng) + h[m][n] * f64::sin(m as f64 * cp.l.lng));
            sumy += m as f64 * p[m][n] *
                (g[m][n] * f64::sin(m as f64 * cp.l.lng) - h[m][n] * f64::cos(m as f64 * cp.l.lng));
        }

        fz += f64::powf(ar, n as f64 + 2.0) * (n as f64 + 1.0) * sumz;
        fx -= f64::powf(ar, n as f64 + 2.0) * sumx;
        fy += f64::powf(ar, n as f64 + 2.0) * sumy;
    }

    cp.dip[hr] = f64::atan(fz / f64::sqrt(f64::powf(fx, 2.0) + f64::powf(fy / f64::cos(cp.l.lat), 2.0)));
    cp.fh[hr] = 2.8 * f64::sqrt(f64::powf(fx, 2.0) + f64::powf(fy / f64::cos(cp.l.lat), 2.0) + f64::powf(fz, 2.0));
}

/// Bilinear interpolation function implementing ITU-R P.1144-5 method
pub fn bilinear_interpolation(ll: f64, lr: f64, ul: f64, ur: f64, r: f64, c: f64) -> f64 {
    ll * ((1.0 - r) * (1.0 - c)) +
    ul * (r * (1.0 - c)) +
    lr * ((1.0 - r) * c) +
    ur * (r * c)
}
