use crate::constants::*;
use crate::e_layer_screening_frequency;
use crate::geometry;
use crate::path_data::*;

/// Calculate basic MUF as described in ITU-R P.533-12 Section 3.5 F2-layer basic MUF
/// and Section 3.3 E-layer basic MUF.
pub fn calculate_muf_basic(path: &mut PathData) {
    // Only do this subroutine if the path is less than or equal to 9000 km if not exit
    if path.distance > 9000.0 {
        return;
    }

    // Initialize lowest-order mode indices
    path.n0_f2 = NO_LOWEST_MODE;
    path.n0_e = NO_LOWEST_MODE;

    // First, determine the F2 layer Basic MUF
    // Determine the mirror reflection height ( hr ) at the midpoint.
    let hr = f64::min(1490.0 / path.cp[MP].m3kf2 - 176.0, 500.0); // Mirror reflection height

    // Store hr to the control point but use hr in this routine for readability.
    path.cp[MP].hr = hr;

    // The mode is just determined by geometry. I will assume that the earth is a sphere. So find the lowest mode by geometry.
    // The mirror reflection height will tell you if multiple hops are necessary.
    // Determine where the horizon is. If the horizon is less than half the distance you must have multiple hops.
    // Try the other modes with a max of 6.
    // There is a minimum elevation angle is MINELEANGLE degrees for the short model.
    let minele = MIN_ELE_ANGLE_S * D2R; // Min elevation in P533 was 3 degrees in radians

    // Find the hop distance to the horizon (dh)
    // The angle of incidence for the minimum elevation (minele) and mirror reflection height (hr)
    let aoi = e_layer_screening_frequency::incidence_angle(minele, hr); // Angle of incidence

    // Now solve for the arc length of the half-hop length and multiply by the earth's radius.
    // Then to find the total hop length (dh) multiply by 2.0
    let mut dh = (PI - aoi - (PI / 2.0 + minele)) * R0 * 2.0; // The hop distance

    // *************************** Not Used ******************************************
    // Determine the maximum hop length for the path if the reflection height is 500 km
    let _dhmax =
        (PI - e_layer_screening_frequency::incidence_angle(minele, 500.0) - (PI / 2.0 + minele))
            * R0
            * 2.0; // The maximum hop distance at a 500km reflection height
                   //*********************************************************************************

    // The hop distance cannot be longer than 4000 km.
    dh = f64::min(dh, 4000.0);

    // Find the lowest-order F2 mode
    for n0 in 0..MAX_F2_MDS {
        if dh > path.distance / (n0 as f64 + 1.0) {
            // Is the mirror reflection height horizon less than the n0 hop distance?
            // At this point lowest-order mode is known – store it.
            path.n0_f2 = n0 as i32;
            break; // You have found the lowest-order mode. Jump out of the loop.
        }
    }

    // Check to see if a low order mode exists F2 layers. If it does not skip Sections 3.5.1 and 3.5.2
    // since they are dependent on the lowest order F2 mode existing
    if path.n0_f2 != NO_LOWEST_MODE {
        let n0 = path.n0_f2 as usize;

        // There will be several calculations that follow that require dmax and B calculated at the midpoint.
        // Both will be calculated here for clarity.
        // The calculation of dmax is dependent on the calculation of B, Eqn (6) in
        // Section 3.5.1.1 Paths up to dmax (km) ITU-R P.533-12.
        // The dmax is determined at mid-path and then stored to the path structure. The path->dmax
        // will be used elsewhere in p533().

        // Limit dmax to 4000 km
        path.dmax = f64::min(calc_dmax(&mut path.cp[MP]), 4000.0);

        // For readability use a local variable.
        let dmax = path.dmax; // path->dmax for readability

        // 3.5.1 Lowest-order mode
        if path.distance <= path.dmax {
            // 3.5.1.1 Paths up to dmax (km)
            // For this calculation the hop distance is path->distance/(n0+1).
            // The path->dmax has already been calculated at the midpoint for the remainder of these calculations.
            let b_mp = calc_b_mut(&mut path.cp[MP]);
            let n0f2dmuf = calc_f2dmuf(
                &path.cp[MP],
                path.distance / (n0 as f64 + 1.0),
                path.dmax,
                b_mp,
            ); // Lowest-order mode basic MUF

            // The lowest-order F2 mode has been determined. Load it into the path structure
            path.md_f2[path.n0_f2 as usize].bmuf = n0f2dmuf; // lowest-order mode basic MUF

            // The n0F2DMUF is the basic F2 MUF of the path.
            path.bmuf = n0f2dmuf;

            // NOTE: There will be no T + d0/2 and R - d0/2 control points for this path
        } else {
            // 3.5.1.2 Paths longer than dmax (km)
            // In this case we have to load two new control points at T - d sub 0/2 and R - d sub 0/2.
            // To find these new locations,
            // first determine the fractional distances and then find the point on the great circle between tx and rx.
            let fracd = 1.0 / (2.0 * (n0 as f64 + 1.0)); // T + d sub 0/2 as a fraction of the total path length
                                                         // Fractional distance
            geometry::great_circle_point(
                path.l_tx,
                path.l_rx,
                &mut path.cp[TD02],
                path.distance,
                fracd,
            );
            let fracd = 1.0 - 1.0 / (2.0 * (n0 as f64 + 1.0)); // R - d sub 0/2 as a fraction of the total path length
            geometry::great_circle_point(
                path.l_tx,
                path.l_rx,
                &mut path.cp[RD02],
                path.distance,
                fracd,
            );
            // All distances for the control points are relative to the tx.

            // Find foF2, M(3000)F2 and foE these control points.
            crate::calculate_cp_parameters::calculate_cp_parameters_impl(
                &mut path.cp[TD02],
                &path.fof2,
                &path.m3kf2,
                path.hour,
                path.ssn,
                path.month,
            );
            crate::calculate_cp_parameters::calculate_cp_parameters_impl(
                &mut path.cp[RD02],
                &path.fof2,
                &path.m3kf2,
                path.hour,
                path.ssn,
                path.month,
            );

            // Determine the F2 basic MUF at each control point
            // For these control points calculate the basic MUF at F2(dmax)MUF
            // Note in this case for equation (3) in P.533-12 section 3.5.1.1

            // Calculate B values first (this stores x in each control point)
            let b_td02 = calc_b_mut(&mut path.cp[TD02]);
            let b_rd02 = calc_b_mut(&mut path.cp[RD02]);

            let f2dmuf = [
                calc_f2dmuf(
                    &path.cp[TD02],
                    path.distance / (n0 as f64 + 1.0),
                    dmax,
                    b_td02,
                ),
                calc_f2dmuf(
                    &path.cp[RD02],
                    path.distance / (n0 as f64 + 1.0),
                    dmax,
                    b_rd02,
                ),
            ];

            path.md_f2[n0].bmuf = f64::min(f2dmuf[0], f2dmuf[1]); // Basic MUF is the lower of the two control point MUFs.

            // The path MUF is small of the two MUFs calculated at the control points T + d0/2 and R - d0/2
            path.bmuf = path.md_f2[n0].bmuf;
        }

        // 3.5.2 Higher-order modes (paths up to 9 000 km)
        if path.distance <= 9000.0 {
            for n in (n0 + 1)..MAX_F2_MDS {
                if path.distance <= path.dmax {
                    // 3.5.2.1 Paths up to dmax (km)
                    let b_mp = calc_b_mut(&mut path.cp[MP]);
                    path.md_f2[n].bmuf =
                        calc_f2dmuf(&path.cp[MP], path.distance / (n as f64 + 1.0), dmax, b_mp);
                } else {
                    //3.5.2.2 Paths longer than dmax (km)
                    // The higher-order modes calculated here imply that the lowest-order mode has already set up the
                    // two new control points and the basic MUF has been determined. The path->BMUF is then n0F2(D)MUF.
                    // Now find the mode basic MUF as a function of F2(dmax[CP[MP]])MUF and a scaling factor.
                    // The scaling factor is
                    //     Mn/Mn0 = nF2(D)MUF/n0F2(D)MUF
                    //        where d is the hop length of the mode and D is the total distance.
                    //**************************************************************
                    // Note: For this calculation dmax is allowed to exceed 4000 km
                    //**************************************************************

                    let dmax_td02 = calc_dmax(&mut path.cp[TD02]);
                    let b_td02_n0 = calc_b_from_stored_x(&path.cp[TD02]);
                    let dmax_rd02 = calc_dmax(&mut path.cp[RD02]);
                    let b_rd02_n0 = calc_b_from_stored_x(&path.cp[RD02]);

                    let mn0 = [
                        calc_f2dmuf(
                            &path.cp[TD02],
                            path.distance / (n0 as f64 + 1.0),
                            dmax_td02,
                            b_td02_n0,
                        ),
                        calc_f2dmuf(
                            &path.cp[RD02],
                            path.distance / (n0 as f64 + 1.0),
                            dmax_rd02,
                            b_rd02_n0,
                        ),
                    ];

                    let dmax_td02_n = calc_dmax(&mut path.cp[TD02]);
                    let b_td02_n = calc_b_from_stored_x(&path.cp[TD02]);
                    let dmax_rd02_n = calc_dmax(&mut path.cp[RD02]);
                    let b_rd02_n = calc_b_from_stored_x(&path.cp[RD02]);

                    let mn = [
                        calc_f2dmuf(
                            &path.cp[TD02],
                            path.distance / (n as f64 + 1.0),
                            dmax_td02_n,
                            b_td02_n,
                        ),
                        calc_f2dmuf(
                            &path.cp[RD02],
                            path.distance / (n as f64 + 1.0),
                            dmax_rd02_n,
                            b_rd02_n,
                        ),
                    ];

                    path.md_f2[n].bmuf = path.bmuf * f64::min(mn[0] / mn0[0], mn[1] / mn0[1]);
                }
            }
        }
    }
    // End calculation for F2 Layer Basic MUF

    // E layer basic MUF calculation
    if path.distance < 4000.0 {
        // Determine the lowest order E mode
        let hr = 110.0; // The mirror reflection height is 110 km for the E layer basic MUF calculation

        // Determine the elevation angle
        let delta = minele; // Elevation angle

        // The angle of incidence for the minimum elevation (minele) and mirror reflection height (hr)
        let aoi = e_layer_screening_frequency::incidence_angle(delta, hr);

        // Now solve for the arc length of the half-hop length and multiply by the earth's radius. Then to find the total
        // hop length (dh) multiply by 2.0.
        let mut dh = (PI - aoi - (PI / 2.0 + minele)) * R0 * 2.0;

        // The hop distance can't be longer than 4000 km.
        dh = f64::min(dh, 4000.0);

        for n0 in 0..3 {
            if dh > path.distance / (n0 as f64 + 1.0) {
                // Is the mirror reflection height horizon less than the n0 hop distance?
                path.n0_e = n0 as i32;
                break; // You have found the lowest-order mode. Jump out of the loop.
            }
        }

        // Is there a lowest order E mode?
        if path.n0_e != NO_LOWEST_MODE {
            let n0 = path.n0_e as usize;

            // There are three E paths.
            for n in n0..MAX_E_MDS {
                // Save the reflection height, although for E layers it is always 110.0 km.
                path.md_e[n].hr = hr;
                // Find the hop length for this mode.
                let dh = f64::min(path.distance / (n as f64 + 1.0), 4000.0);
                // The angle of incidence for hr = 110 km.
                // First, find the angle associated with the half-hop distance, psi.
                let _psi = dh / (2.0 * R0); // The angle associated with the hop length d, d = R*psi.
                                            // Find the elevation angle, delta.
                let delta = e_layer_screening_frequency::elevation_angle(dh, hr);
                // The incident angle i110 is
                let i110 = e_layer_screening_frequency::incidence_angle(delta, hr); // The incident angle for a reflection height of 110 km.

                // The control point is the midpoint
                if path.distance < 2000.0 {
                    path.md_e[n].bmuf = path.cp[MP].foe / f64::cos(i110);
                } else if path.distance >= 2000.0 && path.distance <= 4000.0 {
                    path.md_e[n].bmuf =
                        f64::min(path.cp[R1K].foe, path.cp[T1K].foe) / f64::cos(i110);
                } else {
                    // path.distance > 4000.0
                    path.md_e[n].bmuf = 0.0; // The basic MUF is zero past 4000 km.
                                             // This is redundant but to make the point set the basic E layer MUF to 0.0.
                }

                // Determine the lowest-order E layer mode.
                if path.md_e[n].bmuf != 0.0 && path.n0_e == NO_LOWEST_MODE {
                    path.n0_e = n as i32;
                }
            }
        } else {
            // path.distance >= 4000.0
            // There are no E layer modes for path.distance >= 4000.0.
        }
    }
    // End calculation for the E layer Basic MUF.

    // Finally, determine the path basic MUF.
    if path.n0_e != NO_LOWEST_MODE {
        if path.n0_f2 != NO_LOWEST_MODE {
            // F2 and E modes exist
            path.bmuf = f64::max(
                path.md_e[path.n0_e as usize].bmuf,
                path.md_f2[path.n0_f2 as usize].bmuf,
            );
        } else {
            // E modes exist. F2 mode do not.
            path.bmuf = path.md_e[path.n0_e as usize].bmuf;
        }
    } else {
        // path.n0_e == NO_LOWEST_MODE
        if path.n0_f2 != NO_LOWEST_MODE {
            // F2 modes exist. E modes do not.
            path.bmuf = path.md_f2[path.n0_f2 as usize].bmuf;
        } else {
            // This is an error condition
            path.bmuf = TOO_BIG;
        }
    }
}

/// Finds the dmax for the path. The dmax can be calculated to be greater than 4000 km.
/// Presently, June 2013, dmax can be greater than 4000 km for the calculation of higher order MUFs.
/// all other locations it will typically be restricted to 4000 km.
pub fn calc_dmax(cp: &mut ControlPt) -> f64 {
    // Now you can find B
    let b = calc_b_mut(cp);

    // Determine d sub max
    let dmax = 4780.0
        + (12610.0 + 2140.0 / cp.x.powi(2) - 49720.0 / cp.x.powi(4) + 688900.0 / cp.x.powi(6))
            * (1.0 / b - 0.303);

    dmax
}

/// Readonly version of calc_dmax that doesn't modify the control point
pub fn calc_dmax_readonly(cp: &ControlPt) -> f64 {
    // Calculate x without storing it
    let x = if cp.foe != 0.0 {
        f64::max(cp.fof2 / cp.foe, 2.0)
    } else {
        2.0
    };

    // Now you can find B
    let b = cp.m3kf2 - 0.124
        + (cp.m3kf2.powi(2) - 4.0) * (0.0215 + 0.005 * f64::sin(7.854 / x - 1.9635));

    // Determine d sub max
    let dmax = 4780.0
        + (12610.0 + 2140.0 / x.powi(2) - 49720.0 / x.powi(4) + 688900.0 / x.powi(6))
            * (1.0 / b - 0.303);

    dmax
}

/// Determines the intermediate values B from Eqn (6) P.533-12
fn calc_b_mut(cp: &mut ControlPt) -> f64 {
    // Determine x which is the ratio of foF2 to foE at the midpoint.
    // Check to see if foE exists
    cp.x = if cp.foe != 0.0 {
        f64::max(cp.fof2 / cp.foe, 2.0)
    } else {
        2.0
    };

    // Now you can find B
    let b = cp.m3kf2 - 0.124
        + (cp.m3kf2.powi(2) - 4.0) * (0.0215 + 0.005 * f64::sin(7.854 / cp.x - 1.9635));

    b
}

/// Determines the intermediate values B from Eqn (6) P.533-12
pub fn calc_b(cp: &ControlPt) -> f64 {
    // Determine x which is the ratio of foF2 to foE at the midpoint.
    // Check to see if foE exists
    let x = if cp.foe != 0.0 {
        f64::max(cp.fof2 / cp.foe, 2.0)
    } else {
        2.0
    };

    // Now you can find B
    let b = cp.m3kf2 - 0.124
        + (cp.m3kf2.powi(2) - 4.0) * (0.0215 + 0.005 * f64::sin(7.854 / x - 1.9635));

    b
}

pub fn calc_b_from_stored_x(cp: &ControlPt) -> f64 {
    // Use the stored x value (should have been set by calc_dmax)
    let b = cp.m3kf2 - 0.124
        + (cp.m3kf2.powi(2) - 4.0) * (0.0215 + 0.005 * f64::sin(7.854 / cp.x - 1.9635));
    b
}

/// Performs Eqn (4) P.533-12
fn calc_cd(d: f64, dmax: f64) -> f64 {
    let z = 1.0 - 2.0 * d / dmax;

    let cd = 0.74 - 0.591 * z - 0.424 * z.powi(2) - 0.090 * z.powi(3)
        + 0.088 * z.powi(4)
        + 0.181 * z.powi(5)
        + 0.096 * z.powi(6);

    cd
}

/// Determines the F2 Layer MUF from Eqn (3) P.533-12
pub fn calc_f2dmuf(cp: &ControlPt, distance: f64, dmax: f64, b: f64) -> f64 {
    // If the distance is less than dmax use the distance
    let d = if distance <= dmax { distance } else { dmax };

    // ITU-R P.533-12 Eqn (4)
    let cd = calc_cd(d, dmax); // Cd at D

    let d = 3000.0; // From ITU-R P.533-12 "C sub 3000 : value of Cd for D = 3 000 km"

    let c3k = calc_cd(d, dmax); // Cd at 3000 km

    let f2dmuf =
        (1.0 + cd / c3k * (b - 1.0)) * cp.fof2 + cp.fh[HR_300_KM] / 2.0 * (1.0 - distance / dmax);

    f2dmuf
}
