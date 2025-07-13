use crate::constants::*;
use crate::path_data::*;
use crate::geometry;

/// Sets the path structure output values to default values
pub fn initialize_path(path: &mut PathData) {
    // It is assumed that ValidatePath() has run by this point to assure that the
    // input data is correct.

    /********************************************************************************************/
    /* First, zero out all variables that can potentially have values in any part of this program. */
    /********************************************************************************************/

    // Initialize the path to be set elsewhere.

    // Clear the variables that will be determined elsewhere.
    path.eirp = TINY_DB as f64;
    path.opmuf = 99.9;
    path.opmuf10 = 99.9;
    path.opmuf90 = 99.9;
    path.distance = 999999.9;
    path.season = 99;

    // Initialize all the data in the structures.

    // Initializing modes
    initialize_modes(&mut path.md_f2);
    initialize_modes_e(&mut path.md_e);

    // End initializing modes

    // To initialize the control point, the distance between the tx and rx needs to be determined.
    // Find the great circle distance between the tx and rx.
    path.distance = geometry::great_circle_distance(path.l_tx, path.l_rx);

    // There is a degenerate case where path distance is zero. 
    // If the distance is zero set it to epsilon as an approximation
    if path.distance == 0.0 {
        path.distance = DBL_EPSILON;
    }

    // Determine if this is a long path. If so, adjust the distance.
    if path.sor_l == LONG_PATH {
        path.distance = R0 * PI * 2.0 - path.distance;
    }

    // Initialize the control points
    initialize_cps(path);
    // End initializing control points

    // Initialize the path variables.
    // For several calculations you need to know the season.
    path.season = what_season(path.cp[MP].l, path.month);
}

/// This routine zeros the five potential control points. The locations of the T - d0/2 and R - d0/2 
/// control points are determined elsewhere because they are dependent on n0 ( the lowest-order F2 mode)
fn initialize_cps(path: &mut PathData) {
    // Initialize five control points
    for i in 0..5 {
        path.cp[i].l.lat = 0.0;
        path.cp[i].l.lng = 0.0;
        path.cp[i].distance = 0.0;
        path.cp[i].foe = 0.0;
        path.cp[i].fof2 = 0.0;
        path.cp[i].m3kf2 = 0.0;
        path.cp[i].fh[HR_300_KM] = 0.0;
        path.cp[i].dip[HR_300_KM] = 0.0;
        path.cp[i].fh[HR_100_KM] = 0.0;
        path.cp[i].dip[HR_100_KM] = 0.0;
        path.cp[i].sun.lsn = 0.0;
        path.cp[i].sun.lsr = 0.0;
        path.cp[i].sun.lss = 0.0;
        path.cp[i].sun.sza = 0.0;
        path.cp[i].sun.decl = 0.0;
        path.cp[i].sun.eot = 0.0;
        path.cp[i].sun.sha = 0.0;
        path.cp[i].sun.ha = 0.0;
        path.cp[i].ltime = 0.0;
        path.cp[i].hr = 0.0;
        path.cp[i].x = 0.0;
    }

    /**************************************************************/
    /* Now calculate everything that can be before p533() begins. */
    /**************************************************************/
    // Now find MP, T1k and R1k control points.
    // First, determine the fractional distances then find the point on the great circle between tx and rx.
    // All distances for the control points are relative to the tx.
    // Mid-point control point M
    // MP control point which is always used in the calculation of p533()
    let fracd = 0.5; // fractional distance
    geometry::great_circle_point(path.l_tx, path.l_rx, &mut path.cp[MP], path.distance, fracd);

    // Find foF2, M(3000)F2 and foE the MP control point
    crate::calculate_cp_parameters::calculate_cp_parameters_impl(
        &mut path.cp[MP], 
        &path.fof2, 
        &path.m3kf2, 
        path.hour, 
        path.ssn, 
        path.month
    );

    // The next two control points depend on the total path length. If the path is not at least 2000 km then there is no
    // point in determining control points 1000 km from each end. 
    if path.distance >= 2000.0 {
        // R1k Control point - Fractional distance R - 1000
        let fracd = (path.distance - 1000.0) / path.distance;
        geometry::great_circle_point(path.l_tx, path.l_rx, &mut path.cp[R1K], path.distance, fracd);

        // T1k Control point - Fractional distance T + 1000
        let fracd = 1000.0 / path.distance;
        geometry::great_circle_point(path.l_tx, path.l_rx, &mut path.cp[T1K], path.distance, fracd);

        // Find foF2, M(3000)F2 and foE these control points
        crate::calculate_cp_parameters::calculate_cp_parameters_impl(
            &mut path.cp[T1K], 
            &path.fof2, 
            &path.m3kf2, 
            path.hour, 
            path.ssn, 
            path.month
        );
        crate::calculate_cp_parameters::calculate_cp_parameters_impl(
            &mut path.cp[R1K], 
            &path.fof2, 
            &path.m3kf2, 
            path.hour, 
            path.ssn, 
            path.month
        );
    }

    // Note the local sunrise, sunset and noon at the CPs are found when IonParameters() is run.

    // The other control points at T + d0/2 and R - d0/2 require that the lowest-order propagating mode
    // be determined. This calculation is performed in MUFBasic().
}

/// Initializes the F2 modes
fn initialize_modes(m: &mut [Mode; MAX_F2_MDS]) {
    for i in 0..MAX_F2_MDS {
        m[i].bmuf = 0.0;
        m[i].muf90 = 0.0;
        m[i].muf50 = 0.0;
        m[i].muf10 = 0.0;
        m[i].opmuf = 0.0;
        m[i].opmuf10 = 0.0;
        m[i].opmuf90 = 0.0;
        m[i].fprob = 0.0;
        m[i].deltal = 0.0;
        m[i].deltau = 0.0;
        m[i].fs = 0.0;
        m[i].lb = -TINY_DB as f64;
        m[i].ew = TINY_DB as f64;
        m[i].prw = TINY_DB as f64;
        m[i].grw = TINY_DB as f64;
        m[i].hr = 0.0;
        m[i].tau = 0.0;
        m[i].ele = 0.0;
        m[i].mc = FALSE;
    }
}

/// Initializes the E modes
fn initialize_modes_e(m: &mut [Mode; MAX_E_MDS]) {
    for i in 0..MAX_E_MDS {
        m[i].bmuf = 0.0;
        m[i].muf90 = 0.0;
        m[i].muf50 = 0.0;
        m[i].muf10 = 0.0;
        m[i].opmuf = 0.0;
        m[i].opmuf10 = 0.0;
        m[i].opmuf90 = 0.0;
        m[i].fprob = 0.0;
        m[i].deltal = 0.0;
        m[i].deltau = 0.0;
        m[i].fs = 0.0;
        m[i].lb = -TINY_DB as f64;
        m[i].ew = TINY_DB as f64;
        m[i].prw = TINY_DB as f64;
        m[i].grw = TINY_DB as f64;
        m[i].hr = 0.0;
        m[i].tau = 0.0;
        m[i].ele = 0.0;
        m[i].mc = FALSE;
    }
}

/// Determines the month and latitude dependent index which represents the season.
/// The index is used to locate appropriate variables in the foF2var array.
fn what_season(l: Location, month: i32) -> i32 {
    let mut season = EQUINOX as i32;

    if l.lat >= 0.0 {
        // Northern hemisphere and the equator
        match month {
            NOV | DEC | JAN | FEB => {
                season = WINTER as i32;
            },
            MAR | APR | SEP | OCT => {
                season = EQUINOX as i32;
            },
            MAY | JUN | JUL | AUG => {
                season = SUMMER as i32;
            },
            _ => {},
        }
    } else {
        // Southern Hemisphere 
        match month {
            MAY | JUN | JUL | AUG => {
                season = WINTER as i32;
            },
            MAR | APR | SEP | OCT => {
                season = EQUINOX as i32;
            },
            NOV | DEC | JAN | FEB => {
                season = SUMMER as i32;
            },
            _ => {},
        }
    }

    season
}
