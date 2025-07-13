use crate::path_data::PathData;

/// Read ionospheric parameters from binary files
/// This reads Peter Suessman's ionospheric parameters from the program iongrid
/// The file format is binary with FORTRAN record structure
pub fn read_ion_parameters_bin(month: i32, fof2: &mut Vec<Vec<Vec<Vec<f32>>>>, m3kf2: &mut Vec<Vec<Vec<Vec<f32>>>>) -> bool {
    /*
     * ReadIonParametersBin() is a routine to read ionospheric parameters from a file into arrays necessary for the ITU-R P.533 
     *		calculation engine. All of the input data here that is "hard coded" will be passed presumably to the final version
     *		of the P.533 engine.
     *
     *	This routine reads Peter Suessman's ionospheric parameters from the program iongrid.
     *	The file that is read is a text file. Once the file is read, the result is stored in two arrays, foF2 and M3kF2, 
     *	that will be passed to the ITU HFProp engine to the propagation prediction
     *	The arrays are of the format
     *
     *			foF2[hour][longitude][latitude][SSN] & M3kF2[hour][longitude][latitude][SSN]
     *
     *				where
     *
     *					hour =		0 to 23
     *					longitude = 0 to 240 in 1.5-degree increments from 180 degrees West
     *					latitude =	0 to 120 in 1.5-degree increments from 90 degrees South
     *					SSN =		0 to 1 - 0 and 100 12-month smoothed sun spot number
     *
     *	The eccentricities of how the input file was created are based on P.1239 and Suessman's code, which is based
     *	on work of several administrations. These ionospheric parameter files are based on CCIR spherical harmonic coefficients from 
     *	the 1958 Geophysical year. Please refer to P.1239 for details on how to convert between the coefficients and foF2 and M(3000)F2
     */

    // At present the path structure is not used but is passed in so that it can select the correct map file based on month.
    // The dimensions of the array are fixed by Suessman's file generating program "iongrid".
    // Eventually it would be nice if these were not fixed values so that other resolutions could be used. 
    const HRS: usize = 24; // 24 hours
    const LNG: usize = 241; // 241 longitudes at 1.5-degree increments
    const LAT: usize = 121; // 121 latitudes at 1.5-degree increments
    const SSN: usize = 2; // 2 SSN (12-month smoothed sun spot numbers) high and low
    
    let num_fof2 = HRS * LNG * LAT * SSN;
    
    // Use include_bytes! to embed the file at compile time
    let data = match month {
        0 => include_bytes!("../static_data/ionos01.bin"),
        1 => include_bytes!("../static_data/ionos02.bin"),
        2 => include_bytes!("../static_data/ionos03.bin"),
        3 => include_bytes!("../static_data/ionos04.bin"),
        4 => include_bytes!("../static_data/ionos05.bin"),
        5 => include_bytes!("../static_data/ionos06.bin"),
        6 => include_bytes!("../static_data/ionos07.bin"),
        7 => include_bytes!("../static_data/ionos08.bin"),
        8 => include_bytes!("../static_data/ionos09.bin"),
        9 => include_bytes!("../static_data/ionos10.bin"),
        10 => include_bytes!("../static_data/ionos11.bin"),
        11 => include_bytes!("../static_data/ionos12.bin"),
        _ => return false,
    };
    
    // Skip the first 5 bytes of FORTRAN overhead
    let mut offset = 5;
    
    // Initialize arrays if they're empty
    if fof2.is_empty() {
        fof2.resize(HRS, Vec::new());
        for hour in fof2.iter_mut() {
            hour.resize(LNG, Vec::new());
            for lng in hour.iter_mut() {
                lng.resize(LAT, Vec::new());
                for lat in lng.iter_mut() {
                    lat.resize(SSN, 0.0);
                }
            }
        }
    }
    
    if m3kf2.is_empty() {
        m3kf2.resize(HRS, Vec::new());
        for hour in m3kf2.iter_mut() {
            hour.resize(LNG, Vec::new());
            for lng in hour.iter_mut() {
                lng.resize(LAT, Vec::new());
                for lat in lng.iter_mut() {
                    lat.resize(SSN, 0.0);
                }
            }
        }
    }
    
    // Read foF2 data
    for m in 0..SSN {
        for j in 0..LNG {
            for k in 0..LAT {
                for i in 0..HRS {
                    let buffer_offset = m * LNG * LAT * HRS + j * LAT * HRS + k * HRS + i;
                    let byte_offset = offset + buffer_offset * 4;
                    
                    if byte_offset + 4 <= data.len() {
                        let bytes = &data[byte_offset..byte_offset + 4];
                        let value = f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
                        fof2[i][j][k][m] = value;
                    }
                }
            }
        }
    }
    
    // Skip 10 bytes (5 bytes tail of foF2 record + 5 bytes header for M3kF2 record)
    offset += num_fof2 * 4 + 10;
    
    // Read M3kF2 data
    for m in 0..SSN {
        for j in 0..LNG {
            for k in 0..LAT {
                for i in 0..HRS {
                    let buffer_offset = m * LNG * LAT * HRS + j * LAT * HRS + k * HRS + i;
                    let byte_offset = offset + buffer_offset * 4;
                    
                    if byte_offset + 4 <= data.len() {
                        let bytes = &data[byte_offset..byte_offset + 4];
                        let value = f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
                        m3kf2[i][j][k][m] = value;
                    }
                }
            }
        }
    }
    
    true
}

/// Initialize the ionospheric data arrays for the current month only
pub fn initialize_ionospheric_data(path: &mut PathData) -> bool {
    read_ion_parameters_bin(path.month, &mut path.fof2, &mut path.m3kf2)
}

/// Read ionospheric parameters (legacy interface)
pub fn read_ion_parameters(path: &mut PathData) {
    let _ = initialize_ionospheric_data(path);
}
