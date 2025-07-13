use crate::path_data::PathData;
use crate::constants::*;

/// Interpolate field strength between 7000km and 9000km paths
/// Based on ITU-R P.533-12 Section 5.4 "Paths between 7000 and 9000 km"
pub fn between_7000km_and_9000km(path: &mut PathData) {
    /*     
      Between7000kmand9000km() Interpolates between the path-distances of 7000 and 9000 km. 
             The interpolation method is given in ITU-R P.533-12 Section 5.4 "Paths between 7000 and 9000 km"
     
         INPUT 
             struct PathData *path
     
         OUTPUT
             path-Ei - The interpolated path field strength (dB(1uV/m))
             path.BMUF - The path basic MUF (MHz)

        SUBROUTINES
            CalcB()
            Calcdmax()
            CalcF2DMUF()                     
     */

    let mut bmuf = [0.0; 2]; // Array of basic MUFs at two control points

    if path.distance > 7000.0 && path.distance < 9000.0 {
        let xl = 10.0f64.powf(path.el / 100.0); // Linear field strength of path.El
        let xs = 10.0f64.powf(path.es / 100.0); // Linear field strength of path.Es

        let xi = xs + (path.distance - 7000.0) / 2000.0 * (xl - xs); // Linear interpolated field strength

        path.ei = 100.0 * xi.log10(); // Eqn (42) P.533-12

        // Calculate the basic MUF according to P.533-12 Section 5.4 "Paths between 7000 and 9000 km"
        let dmax = f64::min(crate::muf_basic::calc_dmax_readonly(&path.cp[TD02]), 4000.0); // Maximum hop distance as calculated by Eqn (5) P.533-12
        bmuf[0] = crate::muf_basic::calc_f2dmuf(&path.cp[TD02], path.distance / (path.n0_f2 as f64 + 1.0), dmax, crate::muf_basic::calc_b(&path.cp[TD02]));

        let dmax = f64::min(crate::muf_basic::calc_dmax_readonly(&path.cp[RD02]), 4000.0);
        bmuf[1] = crate::muf_basic::calc_f2dmuf(&path.cp[RD02], path.distance / (path.n0_f2 as f64 + 1.0), dmax, crate::muf_basic::calc_b(&path.cp[RD02]));

        path.bmuf = f64::min(bmuf[0], bmuf[1]);
    }
}
