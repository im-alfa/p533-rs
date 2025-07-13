use crate::constants::*;
use crate::path_data::PathData;

/// Calculates the operational MUF from P.533-12 Section 3.7 "The path operational MUF"
/// This function computes operational MUF values for F2 and E layer modes and determines
/// the overall path operational MUF for different percentiles (50%, 10%, 90%)
pub fn calculate_muf_operational(path: &mut PathData) {
    // Only do this subroutine if the path is less than or equal to 9000 km if not exit
    if path.distance > 9000.0 {
        return;
    }
    
    // Initialize Table 1 ITU-R P.1240-1 "Ratio of the median operational MUF to the median
    // basic MUF for an F2-mode, Rop".
    // Dimensions: [power_index][season][day_night]
    let rop = [
        [
            [1.20, 1.30], // Winter: [Day, Night]
            [1.15, 1.25], // Equinox: [Day, Night] 
            [1.10, 1.20], // Summer: [Day, Night]
        ],
        [
            [1.15, 1.25], // Winter: [Day, Night]
            [1.20, 1.30], // Equinox: [Day, Night]
            [1.25, 1.35], // Summer: [Day, Night]
        ]
    ];
    
    // Determine if the time of interest is night or day. It is either Day, when Sunset is greater than Sunrise with local time in between
    // and it is not Night, if not then it is Night.
    let time = if path.cp[MP].ltime < path.cp[MP].sun.lss 
                && path.cp[MP].ltime > path.cp[MP].sun.lsr
                && !(path.cp[MP].ltime > path.cp[MP].sun.lss && path.cp[MP].ltime < path.cp[MP].sun.lsr) {
        DAY
    } else {
        NIGHT // Night
    };
    
    // Determine the EIRP index power
    let power = if path.eirp <= 30.0 {
        0
    } else {
        0
    };
    
    // Initialize the OPMUF extrema for the F2 Layer
    let mut opf2muf = 0.0;
    let mut opf2muf10 = 0.0;
    let mut opf2muf90 = 0.0; // F2 layer operational MUFs
    
    // 6 F2 Modes OPMUF.
    // Set the OPMUF for the ith F2 mode then determine if it are the largest for all existant F2 modes
    for i in 0..MAX_F2_MDS {
        if path.md_f2[i].bmuf != 0.0 {
            // The mode exists
            path.md_f2[i].opmuf = path.md_f2[i].muf50 * rop[power][path.season as usize][time];
            path.md_f2[i].opmuf10 = path.md_f2[i].opmuf * path.md_f2[i].deltau;
            path.md_f2[i].opmuf90 = path.md_f2[i].opmuf * path.md_f2[i].deltal;
            
            // Now assume that the same procedure that applies to the basic MUF and variability MUF applies to the Operation MUF.
            // As was done in the MUFBasic(), determine the largest of the OPMUF, OPMUF10 and OPMUF90.
            // Pick the MUF extrema for the path for the F2 layer.
            if path.md_f2[i].opmuf > opf2muf {
                opf2muf = path.md_f2[i].opmuf;
            }
            if path.md_f2[i].opmuf90 > opf2muf90 {
                opf2muf90 = path.md_f2[i].opmuf90;
            }
            if path.md_f2[i].opmuf10 > opf2muf10 {
                opf2muf10 = path.md_f2[i].opmuf10;
            }
        }
    }
    
    // Initialize the OPMUF extrema for the E Layers
    let mut opemuf = 0.0;
    let mut opemuf10 = 0.0;
    let mut opemuf90 = 0.0; // E layer operational MUFs
    
    // 3 E Modes OPMUF.
    // Set the OPMUF for the ith E mode then determine if it is the largest amongst all existant E modes
    for i in 0..MAX_E_MDS {
        if path.md_e[i].bmuf != 0.0 {
            // The mode exists
            path.md_e[i].opmuf = path.md_e[i].bmuf;
            path.md_e[i].opmuf10 = path.md_e[i].opmuf * path.md_e[i].deltau;
            path.md_e[i].opmuf90 = path.md_e[i].opmuf * path.md_e[i].deltal;
            
            // Now assume that the same procedure that applies to the basic MUF and variability MUF applies to the Operation MUF.
            // As was done in the MUFBasic(), determine the largest of the OPMUF, OPMUF10 and OPMUF90.
            // Pick the MUF extrema for the path for the E layer.
            if path.md_e[i].opmuf > opemuf {
                opemuf = path.md_e[i].opmuf;
            }
            if path.md_e[i].opmuf90 > opemuf90 {
                opemuf90 = path.md_e[i].opmuf90;
            }
            if path.md_e[i].opmuf10 > opemuf10 {
                opemuf10 = path.md_e[i].opmuf10;
            }
        }
    }
    
    // Now find the largest OPMUF of all existant modes
    path.opmuf = f64::max(opemuf, opf2muf); // largest OPMUF
    path.opmuf90 = f64::max(opemuf90, opf2muf90); // largest 90% OPMUF
    path.opmuf10 = f64::max(opemuf10, opf2muf10); // largest 10% OPMUF
}
