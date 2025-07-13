use crate::constants::*;
use crate::path_data::*;

/// Calculate E layer screening frequency
pub fn calculate_e_layer_screening_frequency(path: &mut PathData) {
    // E layer screening is restricted to paths no longer than 4000 km.
    if path.distance > 4000.0 {
        return;
    }

    // Determine the foE for this calculation.
    let foe = if path.distance <= 2000.0 {
        path.cp[MP].foe
    } else {
        f64::max(path.cp[T1K].foe, path.cp[R1K].foe)
    };

    // Now find the E-Layer screening for the F2 modes that exist.
    // The hop of the F2 mode is related to the index i. Since i is an C index use the hop number n = k + 1; 
    if path.n0_f2 != NO_LOWEST_MODE {
        for k in (path.n0_f2 as usize)..MAX_F2_MDS {
            // Determine the hop distance
            let dh = path.distance / (k as f64 + 1.0); // Hop distance

            if path.distance <= path.dmax {
                path.md_f2[k].hr = mirror_reflection_height(path, &path.cp[MP], dh);
            } else if path.distance > path.dmax {
                // In this case you have to find the mirror reflection height at all the control points and take the mean.
                // Assume that the hop distance is path->dmax.
                path.md_f2[k].hr = (mirror_reflection_height(path, &path.cp[TD02], dh) +
                                   mirror_reflection_height(path, &path.cp[MP], dh) +
                                   mirror_reflection_height(path, &path.cp[RD02], dh)) / 3.0;
            }

            // Find the elevation angle from equation 13 Section 5.1 Elevation angle.
            // ITU-R P.533-12
            let deltaf = elevation_angle(dh, path.md_f2[k].hr); // Elevation angle calculated for F2 layer

            // angle of incidence at height hr = 110 km
            let i = incidence_angle(deltaf, 110.0); // Angle of incidence

            // Now calculate the E layer screening frequency.
            if path.distance <= 2000.0 {
                path.md_f2[k].fs = 1.05 * path.cp[MP].foe / f64::cos(i);
            } else {
                // (path.distance > 2000)
                // Use the larger of the foE at the control points 1000 km from either end.
                path.md_f2[k].fs = 1.05 * f64::max(path.cp[T1K].foe, path.cp[R1K].foe) / f64::cos(i);
            }
        }
    }
}

/// Calculates the mirror reflection height by the method in ITU-R P.533-12 Section 5.1 "Elevation angle".
fn mirror_reflection_height(path: &PathData, cp: &ControlPt, d: f64) -> f64 {
    let mut hr = 0.0; // Mirror reflection height

    // Determine the critical frequency ratio
    let x = cp.fof2 / cp.foe;

    let y = f64::max(x, 1.8);

    let delta_m = 0.18 / (y - 1.4) + 0.096 * (f64::min(path.ssn as f64, 160.0) - 25.0) / 150.0;

    let xr = path.frequency / cp.fof2;

    let h = 1490.0 / (cp.m3kf2 + delta_m) - 316.0;

    if x > 3.33 && xr >= 1.0 {
        // a)
        let e1 = -0.09707 * xr.powi(3) + 0.6870 * xr.powi(2) - 0.7506 * xr + 0.6;

        let f1 = if xr <= 1.71 {
            -1.862 * xr.powi(4) + 12.95 * xr.powi(3) - 32.03 * xr.powi(2) + 33.50 * xr - 10.91
        } else {
            // (xr > 1.71)
            1.21 + 0.2 * xr
        };

        let g = if xr <= 3.7 {
            -2.102 * xr.powi(4) + 19.50 * xr.powi(3) - 63.15 * xr.powi(2) - 44.73
        } else {
            19.25
        };

        let ds = 160.0 + (h + 43.0) * g;
        let a = (d - ds) / (h + 140.0);
        let a1 = 140.0 + (h - 47.0) * e1;
        let b1 = 150.0 + (h - 17.0) * f1 - a1;

        let h = if b1 >= 0.0 && a >= 0.0 {
            a1 + b1 * 2.4_f64.powf(-a)
        } else {
            a1 + b1
        };

        hr = f64::min(h, 800.0);
    } else if x > 3.33 && xr < 1.0 {
        // b)
        let z = f64::max(xr, 0.1);
        let e2 = 0.1906 * z.powi(2) + 0.00583 * z + 0.1936;
        let a2 = 151.0 + (h - 47.0) * e2;
        let f2 = 0.645 * z.powi(2) + 0.883 * z + 0.162;
        let b2 = 141.0 + (h - 24.0) * f2 - a2;
        let df = f64::min(0.115 * d / (z * (h + 140.0)), 0.65);
        let b = -7.535 * df.powi(4) + 15.75 * df.powi(3) - 8.834 * df.powi(2) - 0.378 * df + 1.0;

        let h = if b2 >= 0.0 {
            a2 + b2 * b
        } else {
            a2 + b2
        };

        hr = f64::min(h, 800.0);
    } else {
        // c) x <= 3.33
        let j = -0.7126 * y.powi(3) + 5.863 * y.powi(2) - 16.13 * y + 16.07;
        let u = 8.0e-5 * (h - 80.0) * (1.0 + 11.0 * y.powf(-2.2)) + 1.2e-3 * h * y.powf(-3.6);
        hr = f64::min(115.0 + h * j + u * d, 800.0);
    }

    hr
}

/// Determines the elevation angle from P.533-12 equation (13) Section 5.1 Elevation angle
/// given the hop distance (dh) and the mirror reflection height (hr)
pub fn elevation_angle(dh: f64, hr: f64) -> f64 {
    let ele = 1.0 / f64::tan(dh / (2.0 * R0)) - R0 / (R0 + hr) / f64::sin(dh / (2.0 * R0));
    f64::atan(ele)
}

/// Determine the angle of incidence from P.533-12 equation (12).
/// Section 4 E-layer maximum screening frequency (fs) given the 
/// mirror reflection height (hr) and the elevation angle (deltaf)
pub fn incidence_angle(deltaf: f64, hr: f64) -> f64 {
    f64::asin(R0 * f64::cos(deltaf) / (R0 + hr))
}
