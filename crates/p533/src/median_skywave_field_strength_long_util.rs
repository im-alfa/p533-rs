use crate::path_data::PathData;
use crate::constants::*;
use crate::geometry;
use crate::e_layer_screening_frequency;
use crate::calculate_cp_parameters;
use std::f64::consts::PI;

// Local Define
const NOIL: f64 = -0.17;
const MINELEANGLEL: f64 = 3.0; // Minimum elevation angle for long paths

// Control point array, CP, defines to enhance readability.
const MAXCP: usize = 28; // There are a potential 26 90 km penetration points and 2 control points from Table 1a.
const TD_M2: usize = 26; // For this routine this will be the index to the Control point at T + d0/2.
const RD_M2: usize = 27; // For this routine this will be the index to the Control point at R - d0/2.

// Define for FindfL()
const NOTIME: i32 = 99;

// Define for WinterAnomaly()
const SOUTH: usize = 1;
const NORTH: usize = 0;

/// Median sky-wave field strength calculation for long paths (> 9000 km)
/// Based on ITU-R P.533-12 Section 5.3 "Paths longer than 7000 km"
pub fn median_skywave_field_strength_long(path: &mut PathData) {
    /*     
      MedianSkywaveFieldStrengthLong() - Determines the signal strength for paths greater than 9000 km in accordance with 
             P.533-12 Section 5.3 "Paths longer than 7000 km". This routine is patterned after the FTZ() method found in REC533(). 
             The basis of this routine here and in the algorithm in P.533-12 is from work by	Thomas Dambolt and Peter Suessman. 
     
             INPUT
                 struct PathData *path
     
             OUTPUT
                 path.El - Median field strength (dB(1uV/m))
     
             SUBROUTINES
                ElevationAngle()
                ZeroCP()
                GreatCirclePoint()
                CalculateCPParameters()
                AntennaGain08()
                findfM()
                findfL()
     */

    // fL Calculation

    // fM Calculation

    let mut cp = vec![vec![crate::path_data::ControlPt::default(); 24]; MAXCP]; // Temp

    // Initialize variables
    let mut elevation = 360.0; // Antenna elevation

    // This procedure applies only to paths greater than 7000 km. 
    //   i)  For paths greater than 7000 km the field strength is interpolated 
    //	     with the short (< 7000 km) path method.
    //   ii) For paths greater than 9000 km this is the only method used
    if path.distance >= 7000.0 {
        // For both the fM and fL calculation the reflection height is 300 km
        const HR: f64 = 300.0; // Mirror reflection height

        // fL Calculation Geometry
        // The following determines the hops for the penetration points
        // The number of hops is determined by finding equal hops
        // that are 3000 km or less
        let mut n = 0;    // Temp
        while path.distance / (n as f64 + 1.0) > 3000.0 {
            n += 1;
        }

        // Number of hops
        let n_l = n; // Number of hops 

        // Hop distance
        let d_l = path.distance / (n_l as f64 + 1.0); // Hop distance

        // Determine the elevation angle for fL hops
        let delta_l = e_layer_screening_frequency::elevation_angle(d_l, HR); // Elevation angle

        // fM Calculation Geometry
        // The following determines the hops for the control points
        // The number of hops is determined by finding equal hops
        // that are 4000 km or less. For this calculation the minimum
        // elevation angle is taken into account. If the the hops are 
        // determineded to be 4000 km or less and the elevation angle is 
        // found to be below the minimum another hop is added
        n = 0;
        while path.distance / (n as f64 + 1.0) > 4000.0 {
            n += 1;
        }

        // Number of hops
        let mut n_m = n; // Number of hops

        // Hop distance
        let mut d_m = path.distance / (n_m as f64 + 1.0); // Hop distance

        // Determine the elevation angle
        let mut delta_m = e_layer_screening_frequency::elevation_angle(d_m, HR); // Elevation angle

        // Is the calculated elevation angle more than the minimum
        if delta_m < MINELEANGLEL * D2R {
            n_m += 1; // Add a hop

            // Hop distance
            d_m = path.distance / (n_m as f64 + 1.0);
            // Elevation angle
            delta_m = e_layer_screening_frequency::elevation_angle(d_m, HR);
        }

        /**********************************************************************************************************
           Control and Penetration point initialization for the reference frequencies

           For the fL calculation, there will be required 24 hours of data at the locations of interest.
           Both the calculation of the upper and lower reference frequencies, fM and fL respectively,
           require 24 hours of data. 
           For the calculation of fM, the 24 hours of data are required at the control points described in 
           Table 1a) P.533-12. 
           For the calculation of fL, the 24 hours of data are required at the 90 km penetration points. There are 
           potentially 13 hops tx to rx since the circumference of the earth is 40,075.16 and the hop length 
           for this routine is 3000. For every hop there are two penetrations of the 90 km D layer.
           Consequently, 24 x 2 (fM) and 24 x 26 (fL) control points are needed. 
           Therefore, for this calculation 24 x 15 control points are required. 
         *********************************************************************************************************/
        // Determine the 90 km penetration points 
        let i90 = e_layer_screening_frequency::incidence_angle(delta_l, 90.0); // Angle of incidence at a height of 90 km

        // Determine where the rays penetrate the 90 km height to calculate the dips.
        // Find the 90-km height half-hop distance.
        let phi = PI / 2.0 - delta_l - i90; // 90 km penetration angle

        let dh90 = R0 * phi; // 90 km penetration distance

        // The path structure is used to determine the data at the control points 
        // Store the path.hour
        let hour = path.hour; // Temp

        for j in 0..24 { // hours		
            path.hour = j as i32;

            for i in 0..=n_l { // 90-km penetration points 
                // Zero the elements of the two control points
                zero_cp(&mut cp[2 * i][j]);
                zero_cp(&mut cp[2 * i + 1][j]);

                // There are two control points per hop.
                // First the end nearest the tx for this hop.
                let fracd = (i as f64 * d_l + dh90) / path.distance;
                geometry::great_circle_point(path.l_tx, path.l_rx, &mut cp[2 * i][j], path.distance, fracd);
                calculate_cp_parameters::calculate_cp_parameters_impl(&mut cp[2 * i][j], &path.fof2, &path.m3kf2, path.hour, path.ssn, path.month);

                cp[2 * i][j].hr = 90.0;

                // Next the end nearest to the receiver for this hop
                let fracd = ((i as f64 + 1.0) * d_l - dh90) / path.distance;
                geometry::great_circle_point(path.l_tx, path.l_rx, &mut cp[2 * i + 1][j], path.distance, fracd);
                calculate_cp_parameters::calculate_cp_parameters_impl(&mut cp[2 * i + 1][j], &path.fof2, &path.m3kf2, path.hour, path.ssn, path.month);

                cp[2 * i + 1][j].hr = 90.0;
            }

            // Initialize control points (T + d0/2 & R - d0/2) from Table 1a) as the last two control points in the array.
            // First determine the fractional distances and then find the point on the great circle between tx and rx.
            let fracd = 1.0 / (2.0 * (n_m as f64 + 1.0)); // T + d0/2 as a fraction of the total path length
            geometry::great_circle_point(path.l_tx, path.l_rx, &mut cp[TD_M2][j], path.distance, fracd);
            let fracd = 1.0 - 1.0 / (2.0 * (n_m as f64 + 1.0)); // R - d0/2 as a fraction of the total path length
            geometry::great_circle_point(path.l_tx, path.l_rx, &mut cp[RD_M2][j], path.distance, fracd);
            // All distances for the control points are relative to the tx.

            // Find foF2, M(3000)F2 and foE these control points.
            calculate_cp_parameters::calculate_cp_parameters_impl(&mut cp[TD_M2][j], &path.fof2, &path.m3kf2, path.hour, path.ssn, path.month);
            calculate_cp_parameters::calculate_cp_parameters_impl(&mut cp[RD_M2][j], &path.fof2, &path.m3kf2, path.hour, path.ssn, path.month);

            cp[TD_M2][j].x = 0.0;
            cp[TD_M2][j].foe = 0.0;
            cp[TD_M2][j].hr = 300.0; // For this calculation the reflection height is fixed at 300 km.

            cp[RD_M2][j].x = 0.0;
            cp[RD_M2][j].foe = 0.0;
            cp[RD_M2][j].hr = 300.0; // For this calculation the reflection height is fixed at 300 km.
        }

        // Restore the path.hour
        path.hour = hour;
        /**********************************************************************
           End control point initialization for the reference frequencies.
        ***********************************************************************/
        // The data is now available in the control point array, CP, for the calculation of the reference frequencies.
        // Find the virtual slant range (19)
        let psi = d_m / (2.0 * R0);
        path.ptick = ((2.0 * R0 * (psi.sin() / (delta_m + psi).cos())).abs()) * (n_m as f64 + 1.0);

        // Free space field strength
        path.e0 = 139.6 - 20.0 * path.ptick.log10();

        let a_tx_copy = path.a_tx.clone();
        path.gtl = antenna_gain_08(path, &a_tx_copy, TXTORX, &mut elevation);

        // Focusing on long distance gain limited to 15 dB
        let d = path.distance; // Path distance for focus gain term

        path.gap = f64::min(10.0 * (d / (R0 * (d / R0).sin().abs())).log10(), 15.0);

        path.ly = NOIL;

        // Transmitter power dB (1kW)
        let pt = path.txpower; // Transmitter power

        // Mean gyrofrequency
        path.fh = (cp[TD_M2][path.hour as usize].fh[HR_300_KM] + cp[RD_M2][path.hour as usize].fh[HR_300_KM]) / 2.0;

        // Find the MUF fM
        find_mufs_and_fm(path, &cp, n_m, d_m);

        // Lower frequency
        find_fl(path, &cp, n_l, d_l, path.ptick, path.fh, i90);

        let f = path.frequency; // path.frequency

        // Calculate Etl in sections
        let mut etl = (path.fl + path.fh).powi(2) / (f + path.fh).powi(2) + (f + path.fh).powi(2) / (path.fm + path.fh).powi(2); // Resultant median field strength
        etl *= (path.fm + path.fh).powi(2) / ((path.fm + path.fh).powi(2) + (path.fl + path.fh).powi(2));

        path.f = 1.0 - etl;

        etl = path.e0 * (1.0 - etl);
        etl = etl - 30.0 + pt + path.gtl + path.gap - path.ly;
        path.el = etl;

        /**************************************************************
           End of the calculation for paths greater than 7000 km 
        ***************************************************************/

        // In the case that the path is greater than 9000 km, this is the only method that is used so set the path parameters appropriately
        if path.distance > 9000.0 {
            // The elevation angle for the fM calculation is what reperesents the long path
            path.ele = delta_m;

            // Copy the control points to the path structure
            copy_cp(&cp[RD_M2][path.hour as usize], &mut path.cp[RD02]);
            copy_cp(&cp[TD_M2][path.hour as usize], &mut path.cp[TD02]);

            // Copy the two extreme penetration points to control points in the path structure
            copy_cp(&cp[0][path.hour as usize], &mut path.cp[T1K]);
            copy_cp(&cp[2 * n_l][path.hour as usize], &mut path.cp[R1K]);

            // Path dmax
            path.dmax = 4000.0;
        }
    }
}

fn zero_cp(cp: &mut crate::path_data::ControlPt) {
    cp.l.lat = 0.0;
    cp.l.lng = 0.0;
    cp.distance = 0.0;
    cp.foe = 0.0;
    cp.fof2 = 0.0;
    cp.m3kf2 = 0.0;
    cp.dip = [0.0; 2];
    cp.fh = [0.0; 2];
    cp.ltime = 0.0;
    cp.hr = 0.0;
    cp.x = 0.0;
    cp.sun.lsn = 0.0;
    cp.sun.lsr = 0.0;
    cp.sun.lss = 0.0;
    cp.sun.sza = 0.0;
    cp.sun.decl = 0.0;
    cp.sun.eot = 0.0;
    cp.sun.sha = 0.0;
    cp.sun.ha = 0.0;
}

/// Finds the MUFs and fM
fn find_mufs_and_fm(path: &mut PathData, cp: &[Vec<crate::path_data::ControlPt>], _hops: usize, d_m: f64) {
    /*     
      FindfM() Finds the upper reference frequency, fM, from 24 hours of calculated MUFs.
        Determines the the basic and operation MUFs for the long path. This routine also
        finds the 10% and 90% decile values for the MUF
     
             INPUT
                 struct PathData *path
                 struct ControlPt CP[MAXCP][24] - 24 hours of control point data
                 int hops - Number of hops
                double dh - hop length
     
             OUTPUT
                 path.BMUF
                path.MUF50
                path.MUF10
                path.MUF90
                path.OPMUF
                path.OPMUF10 
                path.OPMUF90
                path.fM

            SUBROUTINES
                Bearing()
                FindfoF2var()
     */

    let mut f4 = [0.0; 2];           // F2(4000)MUF at the control points  
    let mut fz = [0.0; 2];           // F2(Zero)MUF at the control points

    let mut f_bm_min = [0.0; 2];       // The lowest value of f4 in 24 hours at control points

    // Values used in the determination of K
    let w = [0.1, 0.2];
    let x = [1.2, 0.2];
    let y = [0.6, 0.4];
    let mut f_bm = [[0.0; 24]; 2];      // F2(D)MUF at the control points for 24 hours 

    let mut noon = [0; 2];             // temp local noon index 
    
    // The folowing equation for fD comes from 
    // "Predicting the Performance of BAND 7 Communications Systems Using Electronic Computers"
    // NBS Report 7619, D. Lucas, 1963. The coefficients have been changed for  
    // hop distances in km. See page 92 step 12 in section "IX Mathmatical Expressions"
    // Calculate fD "Distance reduction factor"
    let f_d = ((((((-2.40074637494790e-24 * d_m +
                  25.8520201885984e-21) * d_m +
                 -92.4986988833091e-18) * d_m +
                102.342990689362e-15) * d_m +
               22.0776941764705e-12) * d_m +
              87.4376851991085e-9) * d_m +
             29.1996868566837e-6) * d_m;

    // What is the index that will be local noon (UTC) at the T + d0/2 and R - d0/2?
    // The lat and lng of the control points are all the same so the hour is arbitrarily so use 1
    // Subtract one from the value noon[*] since it will be used as an index and not a countable hour 
    noon[0] = (12.0 - cp[TD_M2][1].l.lng / (15.0 * D2R)) as i32 - 1; // The longitude is subtracted from 12 because W is negative
    noon[1] = (12.0 - cp[RD_M2][1].l.lng / (15.0 * D2R)) as i32 - 1; //

    // Roll over the time
    noon[0] = (noon[0] + 24) % 24;
    noon[1] = (noon[1] + 24) % 24;

    // Initialize the variables to find the minimum MUF in 24 hours.
    f_bm_min[0] = 100.0;
    f_bm_min[1] = 100.0;
    for t in 0..24 { // Local time counter
        // Calculate F2(4000)MUF
        f4[0] = 1.1 * cp[TD_M2][t].fof2 * cp[TD_M2][t].m3kf2;
        // Calculate F2(ZERO)MUF
        fz[0] = cp[TD_M2][t].fof2 + 0.5 * cp[TD_M2][t].fh[HR_300_KM];

        // Calculate the Basic MUF
        f_bm[0][t] = fz[0] + (f4[0] - fz[0]) * f_d;
        f_bm_min[0] = f64::min(f_bm[0][t], f_bm_min[0]);

        // Calculate F2(4000)MUF	
        f4[1] = 1.1 * cp[RD_M2][t].fof2 * cp[RD_M2][t].m3kf2;
        // Calculate F2(ZERO)MUF
        fz[1] = cp[RD_M2][t].fof2 + 0.5 * cp[RD_M2][t].fh[HR_300_KM];

        f_bm[1][t] = fz[1] + (f4[1] - fz[1]) * f_d;
        f_bm_min[1] = f64::min(f_bm[1][t], f_bm_min[1]);
    }

    // Before proceeding, finding the forward azimuth at the midpoint is required.
    // The azimuth is used to interpolate the W, X and Y values to calculate K.
    let mut a = geometry::bearing(path.cp[MP].l, path.l_rx); // Azimuth at the midpoint

    // Now use A to interpolate the W, X and Y values.
    if a > PI {
        a -= PI;
    }

    if a >= PI / 2.0 {
        a -= PI / 2.0;
    } else {
        a = PI / 2.0 - a;
    }

    let ew = a / (PI / 2.0); // Interpolation value to fine W, X and Y
    let iw = w[0] * (1.0 - ew) + w[1] * ew; // Interpolated W
    let iy = y[0] * (1.0 - ew) + y[1] * ew; // Interpolated Y
    // Interpolated values of W, X and Y respectively
    let ix = x[0] * (1.0 - ew) + x[1] * ew; // Interpolated X

    // fBM,min has been determined for both control points now K  
    for n in 0..2 { // Local counters
        // determine K
        path.k[n] = 1.2 + iw * (f_bm[n][path.hour as usize] / f_bm[n][noon[n] as usize]) + ix * ((f_bm[n][noon[n] as usize] / f_bm[n][path.hour as usize]).powf(1.0 / 3.0) - 1.0) + iy * (f_bm_min[n] / f_bm[n][noon[n] as usize]).powi(2);
    }

    path.fm = f64::min(path.k[0] * f_bm[0][path.hour as usize], path.k[1] * f_bm[1][path.hour as usize]);

    // Before leaving this routine determine if any of the other MUF parameters should be set
    if path.distance > 9000.0 {
        let smaller_cp: usize; // Index that indicates the control point where the smaller basic MUF 
        if f_bm[0][path.hour as usize] < f_bm[1][path.hour as usize] {
            smaller_cp = TD_M2;
            // This is the Basic MUF for the path if the path is greater than 9000 km, otherwise it is calculated in MUFBasic()
            path.bmuf = f_bm[0][path.hour as usize];
        } else {
            smaller_cp = RD_M2;
            // This is the Basic MUF for the path if the path is greater than 9000 km, otherwise it is calculated in MUFBasic()
            path.bmuf = f_bm[1][path.hour as usize];
        }

        // Determine the MUF deciles
        let mut decile = DL;             // decile flag
        // Find the deltal in the foF2var array
        let deltal = crate::muf_variability::find_fof2_var(path, cp[smaller_cp][path.hour as usize].ltime, cp[smaller_cp][path.hour as usize].l.lat, decile); // lower decile of MUF

        decile = DU; // Upper MUF decile
        // Find the deltau in the foF2var array
        let deltau = crate::muf_variability::find_fof2_var(path, cp[smaller_cp][path.hour as usize].ltime, cp[smaller_cp][path.hour as usize].l.lat, decile); // Upper decile of MUF

        // Determine the decile MUFs
        path.muf50 = path.bmuf;
        path.muf10 = path.muf50 * deltau;
        path.muf90 = path.muf50 * deltal;

        // Operation MUF and decile Operation MUF
        path.opmuf = path.fm;
        path.opmuf10 = path.opmuf * deltau; // largest 10% OPMUF
        path.opmuf90 = path.opmuf * deltal; // largest 90% OPMUF
    }
}

/// Finds the lower reference frequency fL
fn find_fl(path: &mut PathData, cp: &[Vec<crate::path_data::ControlPt>], hops: usize, _dh: f64, ptick: f64, fh: f64, i90: f64) {
    /*	 
        FindfL() - Determines the lower reference frequency from 24 hours of solar zenith angles.
     
            INPUT
                struct PathData *path
                struct ControlPt CP[MAXCP][24] - 24 hours of control point data including solar zenith angles
                int hops - Number of hops
                double dh - Hop distance
                double ptick - Slant range
                double fH - Mean gyrofrequency
                double i90 - Incident angle at 90 km
     
            OUTPUT
                return fL the upper reference frequency
     
            SUBROUTINE
                WinterAnomaly()
     */

    let mut sum_cos_chi = [0.0; 24];
    let mut f_l = [0.0; 24];
    let mut dt: f64;      // Temp

    let mut prev: usize;
    let mut now: usize;

    for t in 0..24 {
        sum_cos_chi[t] = 0.0;
        for i in 0..(2 * (hops + 1)) {
            let chi = cp[i][t].sun.sza; // solar zenith angle
            if chi > 0.0 && chi < PI / 2.0 {
                // Sum the control points
                sum_cos_chi[t] += chi.cos().sqrt();
            }
        }
    }

    // Determine winter-anomaly factor, Aw at the path mid-point
    let aw = winter_anomaly(path.cp[MP].l.lat, path.month);

    let f_ln = (path.distance / 3000.0).sqrt();

    // Calculate all 24 values of the LUF first
    // This is necessary first so that we can find the transition of the LUF from
    // day-to-night and not the night-to-day transition.
    // Note: For this calculation the SSN can be greater than MAXSSN
    for i in 0..24 { // Calculate fL for 24 hours
        f_l[i] = f64::max((5.3 * ((1.0 + 0.009 * path.ssn as f64) * sum_cos_chi[i] / (i90.cos() * (9.5e6 / ptick).ln())).sqrt() - fh) * (aw + 1.0), f_ln);
    }

    let mut tr = NOTIME; // Initialize the reference time indicator.

    // Time that requires special treatment
    // Find the first local time that fL[i] <= 2.0*fLN decay from day-LUF to night-LUF
    for n in 0..24 { // Calculate fL for 24 hours
        now = n;
        prev = (n + 24 - 1) % 24; // Find the previous time and roll over

        if tr == NOTIME { // only find one reference time in 24 hours
            if f_l[prev] >= 2.0 * f_ln && f_l[now] <= 2.0 * f_ln { // Find a time where the LUF is decreasing into night
                tr = now as i32;
                dt = (2.0 * f_ln - f_l[tr as usize]) / (f_l[prev] - f_l[tr as usize]);
                f_l[tr as usize] = 0.7945 * f_l[prev] * (dt * (1.0 - 0.7945) + 0.7945);
                if f_l[now] < f_l[tr as usize] {
                    f_l[now] = f_l[tr as usize];
                }
            }
        }
    }

    if tr != NOTIME {
        for i in 1..4 {
            now = ((tr + i) + 24) as usize % 24; // Find now and roll over
            prev = (now + 24 - 1) % 24; // Find previous time and roll over
            f_l[now] = f64::max(f_l[prev] * 0.7945, f_l[now]);
        }
    }

    // Find fL[] for the "present hour" for further calculations 
    // The "present hour" is calculated at path.hour + 1. 
    // Elsewhere in the code it is assumed that UTC 1 uses the index path.hour = 0
    // to represent time from 0:00 to 0:59 
    // In this calculation UTC 1 implies 1:00 to 1:59 UTC so the "present hour" 
    // is path.hour + 1. 
    // This is a consequence of the routine being based on the Fortran program FTZ() 
    // which of course uses indexing beginning at 1
    // Find the "present hour" and roll it over if necessary
    now = ((path.hour + 1) + 24) as usize % 24;

    // Set the value of fL[]
    path.fl = f_l[now];
}

/// Calculates the winter anomaly factor
fn winter_anomaly(lat: f64, month: i32) -> f64 {
    /*
      WinterAnomaly() calculates the winter anomaly factor for long paths (> 9000 km).
            This function uses the Table 5 P.533-12 to find Aw for any latitude.
            Values other than 60 degrees latitude are determined by interpolation.
     
            INPUT
                double lat - Latitude of the point of interest
                double month - The month of interest
     
            OUTPUT
                return Aw the interpolated winter anomaly factor 
     
            SUBROUTINES
                None
     */

    // Winter-anomaly factor
    let aw = [
        [0.30, 0.00],  // January
        [0.15, 0.00],  // February
        [0.03, 0.00],  // March
        [0.00, 0.03],  // April
        [0.00, 0.15],  // May
        [0.00, 0.30],  // June
        [0.00, 0.30],  // July
        [0.00, 0.15],  // August
        [0.00, 0.03],  // September
        [0.03, 0.00],  // October
        [0.15, 0.00],  // November
        [0.30, 0.00]   // December
    ];

    let i_ns = if lat < 0.0 { SOUTH } else { NORTH }; // North-South index

    let lat_abs = lat.abs(); // All latitudes must be positive
    let lat_deg = lat_abs * R2D;

    if lat_deg <= 30.0 || lat_deg >= 90.0 {
        return 0.0;
    }

    // Interpolate with the peak winter anomaly factor at 60 degrees. 
    if lat_deg < 60.0 {
        aw[(month - 1) as usize][i_ns] * (lat_deg - 30.0) / 30.0
    } else {
        aw[(month - 1) as usize][i_ns] * (90.0 - lat_deg) / 30.0
    }
}

/// Antenna gain calculation for 0-8 degrees
fn antenna_gain_08(path: &mut PathData, ant: &crate::path_data::Antenna, direction: i32, elevation: &mut f64) -> f64 {
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
    let mut g_max = TINY_DB as f64;
    for i in 0..9 {
        let delta = i as f64 * D2R; // Elevation angle
        let g = crate::median_skywave_field_strength_short_util::antenna_gain(path, &ant, delta, direction);
        if g > g_max {
            g_max = g;
            *elevation = i as f64 * D2R;
        }
    }

    g_max
}

/// Copies one control point to another
fn copy_cp(this_cp: &crate::path_data::ControlPt, that_cp: &mut crate::path_data::ControlPt) {
    /*
      CopyCp() Copies one control point to another for the situation where they can not point
            to the same memeory location.
     
            INPUT
                struct ControlPt thisCP - The source control point that contains the information to copy 
                struct ControlPt thatCP - The target conntrol point that will be coppied to. 
     
            OUTPUT
                None 
     
            SUBROUTINES
                None
     */

    that_cp.dip[0] = this_cp.dip[0];
    that_cp.dip[1] = this_cp.dip[1];

    that_cp.distance = this_cp.distance;

    that_cp.fh[0] = this_cp.fh[0];
    that_cp.fh[1] = this_cp.fh[1];

    that_cp.foe = this_cp.foe;
    that_cp.fof2 = this_cp.fof2;
    that_cp.hr = this_cp.hr;
    that_cp.l.lat = this_cp.l.lat;
    that_cp.l.lng = this_cp.l.lng;
    that_cp.ltime = this_cp.ltime;
    that_cp.m3kf2 = this_cp.m3kf2;

    that_cp.sun.sza = this_cp.sun.sza;
    that_cp.sun.sha = this_cp.sun.sha;
    that_cp.sun.decl = this_cp.sun.decl;
    that_cp.sun.eot = this_cp.sun.eot;
    that_cp.sun.ha = this_cp.sun.ha;
    that_cp.sun.lsn = this_cp.sun.lsn;
    that_cp.sun.lsr = this_cp.sun.lsr;
    that_cp.sun.lss = this_cp.sun.lss;

    that_cp.x = this_cp.x;
}
