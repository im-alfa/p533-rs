use crate::constants::*;
use crate::path_data::PathData;

/// Bilinear interpolation function implementing ITU-R P.1144-5 method
/// 
/// # Parameters
/// - `ll`: Lower left neighbor
/// - `lr`: Lower right neighbor  
/// - `ul`: Upper left neighbor
/// - `ur`: Upper right neighbor
/// - `r`: Fraction row
/// - `c`: Fractional column
/// 
/// # Returns
/// Interpolated value
fn bilinear_interpolation(ll: f64, lr: f64, ul: f64, ur: f64, r: f64, c: f64) -> f64 {
    ll * ((1.0 - r) * (1.0 - c)) +
    ul * (r * (1.0 - c)) +
    lr * ((1.0 - r) * c) +
    ur * (r * c)
}

/// Calculates the MUF variability for F2 and E layers including probability calculations
/// This implements the procedures from ITU-R P.533-12 for determining MUF percentiles
/// and mode probability based on frequency and MUF relationships
pub fn calculate_muf_variability(path: &mut PathData) {
    // Only do this subroutine if the path is less than or equal to 9000 km if not exit
    if path.distance > 9000.0 {
        return;
    }

    // Section 3.65 seems to imply that the basic MUF for the path and the 50% MUF are the same so for this calculation set them equal
    path.muf50 = path.bmuf;

    // Calculate F2 MUF variability
    for i in 0..MAX_F2_MDS {
        if path.md_f2[i].bmuf != 0.0 {
            // If the basic MUF is set then the layer exists
            
            // Section 3.6 P.533-12 indicates that the basic MUF and MUF(50) are the same.
            path.md_f2[i].muf50 = path.md_f2[i].bmuf;
            
            // For F2 layer, the decile factors are determined from the foF2var array
            let hour = path.cp[MP].ltime;
            let lat = path.cp[MP].l.lat;
            
            // Get the upper and lower decile factors
            let deltal = find_fof2_var(&path, hour, lat, 0);  // Lower decile (10%)
            let deltau = find_fof2_var(&path, hour, lat, 1);  // Upper decile (90%)
            
            path.md_f2[i].deltal = deltal;
            path.md_f2[i].deltau = deltau;
            
            // Calculate the MUF10 and MUF90 from the basic MUF
            path.md_f2[i].muf10 = path.md_f2[i].deltau * path.md_f2[i].muf50;
            path.md_f2[i].muf90 = path.md_f2[i].deltal * path.md_f2[i].muf50;
            
            // Now determine the probability that the mode can be supported
            if path.frequency < path.md_f2[i].muf50 {
                // Equation 9 P.533-12
                path.md_f2[i].fprob = f64::min(
                    1.3 - 0.8 / (1.0 + (1.0 - path.frequency / path.md_f2[i].muf50) / (1.0 - path.md_f2[i].deltal)),
                    1.0
                );
            } else {
                // Equation 10 P.533-12 (path.frequency >= path.md_f2[i].muf50)
                path.md_f2[i].fprob = f64::max(
                    0.8 / (1.0 + (path.frequency / path.md_f2[i].muf50 - 1.0) / (path.md_f2[i].deltau - 1.0)) - 0.3,
                    0.0
                );
            }
        }
    }
    
    // Calculate E layer MUF variability
    for i in 0..MAX_E_MDS {
        if path.md_e[i].bmuf != 0.0 {
            // If the basic MUF is set then the layer exists
            
            // Section 3.6 P.533-12 indicates that the basic MUF and MUF(50) are the same.
            path.md_e[i].muf50 = path.md_e[i].bmuf;
            
            // For E layer, use fixed decile factors as per ITU-R P.533-12
            path.md_e[i].deltal = 0.95;  // Lower decile (10%)
            path.md_e[i].deltau = 1.05;  // Upper decile (90%)
            
            // Calculate the MUF10 and MUF90 from the basic MUF
            path.md_e[i].muf10 = path.md_e[i].deltau * path.md_e[i].muf50;
            path.md_e[i].muf90 = path.md_e[i].deltal * path.md_e[i].muf50;
            
            // Now determine the probability that the mode can be supported
            if path.frequency < path.md_e[i].muf50 {
                path.md_e[i].fprob = f64::min(
                    1.3 - 0.8 / ((1.0 - path.frequency / path.md_e[i].muf50) / (1.0 - path.md_e[i].deltal)),
                    1.0
                );
            } else {
                path.md_e[i].fprob = f64::max(
                    0.8 / ((path.frequency / path.md_e[i].muf50 - 1.0) / (path.md_e[i].deltau - 1.0)) - 3.0,
                    0.0
                );
            }
        }
    }
    
    // Now assume that the same procedure that applies to the basic MUF does for the variability MUFs
    // As was done in the MUFBasic() determine the largest of the MUF90 and MUF10
    // Pick the MUF extrema for the path for the E and F2 layer
    let mut emuf10 = 0.0; // Temp 10% E layer MUF
    let mut emuf90 = 0.0; // Temp 90% E layer MUF
    
    for i in 0..MAX_E_MDS {
        if path.md_e[i].bmuf != 0.0 {
            // If the Basic MUF is set the layer exists
            if path.md_e[i].muf90 > emuf90 {
                emuf90 = path.md_e[i].muf90;
            }
            if path.md_e[i].muf10 > emuf10 {
                emuf10 = path.md_e[i].muf10;
            }
        }
    }
    
    let mut f2muf10 = 0.0; // Temp 10% F2 layer MUF
    let mut f2muf90 = 0.0; // Temp 90% F2 layer MUF
    
    for i in 0..MAX_F2_MDS {
        if path.md_f2[i].bmuf != 0.0 {
            // If the Basic MUF is set the layer exists
            if path.md_f2[i].muf90 > f2muf90 {
                f2muf90 = path.md_f2[i].muf90;
            }
            if path.md_f2[i].muf10 > f2muf10 {
                f2muf10 = path.md_f2[i].muf10;
            }
        }
    }
    
    // Determine the 90% and 10% MUF amongst all existant modes
    path.muf90 = f64::max(emuf90, f2muf90); // largest 90% MUF
    path.muf10 = f64::max(emuf10, f2muf10); // largest 10% MUF
}

/// Determines the variation in foF2 given the season, hour, latitude and decile
/// This routine uses the bilinear interpolation method in ITU-R P.1144-5
/// 
/// # Parameters
/// - `path`: Path data structure containing foF2var array
/// - `hour`: Hour of interest
/// - `lat`: Latitude of interest  
/// - `decile`: Upper or lower decile index (0 = lower/10%, 1 = upper/90%)
/// 
/// # Returns
/// Interpolated foF2 variability value
pub fn find_fof2_var(path: &PathData, hour: f64, lat: f64, decile: usize) -> f64 {
    let lat = f64::abs(lat / (5.0 * D2R)); // 5 degree increments
    
    // Determine the fractional column and row
    // The distance between indices is 1.0
    let r = lat - f64::floor(lat); // The fractional part of the row
    let c = hour - f64::floor(hour); // The fractional part of the column
    let mut lat_l = f64::floor(lat) as i32; // lower lat
    let mut lat_u = f64::ceil(lat) as i32; // upper lat
    
    let mut r_adj = r;
    if lat_l < 0 {
        lat_l = 18; // rollunder
        // The sense of the fractional row is reversed when negative so fix it
        r_adj = 1.0 - r;
    }
    
    if lat_u > 18 {
        lat_u = 0; // rollover
    }
    
    let mut hour_l = f64::floor(hour) as i32; // lower hour
    let mut hour_u = f64::ceil(hour) as i32; // upper hour
    
    let mut c_adj = c;
    if hour_l < 0 {
        hour_l = 23; // rollunder
        // The sense of the fractional column is reversed when negative so fix it
        c_adj = 1.0 - c;
    }
    
    if hour_u > 23 {
        hour_u = 0; // rollover hour
    }
    
    // Determine the sunspot number index
    let ssn = if path.ssn < 50 {
        0
    } else if path.ssn <= 100 {
        1
    } else {
        2
    };
    
    // Find the neighbors
    let ll = path.fof2var[path.season as usize][hour_l as usize][lat_l as usize][ssn][decile];
    let lr = path.fof2var[path.season as usize][hour_u as usize][lat_l as usize][ssn][decile];
    let ul = path.fof2var[path.season as usize][hour_l as usize][lat_u as usize][ssn][decile];
    let ur = path.fof2var[path.season as usize][hour_u as usize][lat_u as usize][ssn][decile];
    
    let irc = bilinear_interpolation(ll, lr, ul, ur, r_adj, c_adj); // Interpolated value
    
    irc
}
