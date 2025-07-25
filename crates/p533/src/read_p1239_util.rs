use crate::path_data::PathData;

/// Read the file "P1239-3 Decile Factors.txt", which is Table 2 and 3 in ITU-R P1239-2 (10/09)
/// The data in this file are the decile factors for within-the-month variations of foF2
pub fn read_p1239(path: &mut PathData) -> bool {
    /*
     * ReadP1239() - Read the file "P1239-2 Decile Factors.txt", which is Table 2 and 3 in ITU-R P1239-2 (10/09).
     *		The data in this file are the decile factors for within-the-month variations of foF2.
     *
     *			INPUT
     *				struct PathData *path
     *
     *			OUTPUT
     *				data is written into the array path.foF2var
     *
     */

    // Use include_bytes! to embed the file at compile time and handle encoding issues
    let data_bytes = include_bytes!("../static_data/P1239-3 Decile Factors.txt");
    let data = String::from_utf8_lossy(data_bytes);
    let lines: Vec<&str> = data.lines().collect();

    const SEASON: usize = 3; // 3 seasons
                             //		1) WINTER 2) EQUINOX 3) SUMMER
    const LAT: usize = 19; // 19 latitude by 5
                           //      0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90
    const SSN: usize = 3; // 3 SSN ranges
                          //		1) R12 < 50 2) 50 <= R12 <= 100 3) R12 > 100
    const DECILE: usize = 2; // 2 deciles
                             //	1) lower 2) upper

    // Initialize the foF2var array if it's empty
    if path.fof2var.is_empty() {
        path.fof2var.resize(SEASON, Vec::new());
        for season in path.fof2var.iter_mut() {
            season.resize(24, Vec::new()); // 24 hours
            for hour in season.iter_mut() {
                hour.resize(LAT, Vec::new());
                for lat in hour.iter_mut() {
                    lat.resize(SSN, Vec::new());
                    for ssn in lat.iter_mut() {
                        ssn.resize(DECILE, 0.0);
                    }
                }
            }
        }
    }

    let mut current_line = 2; // Skip first 2 lines

    // Now read numbers for the lower decile.
    for n in 0..DECILE {
        // 2 deciles lower and upper
        for i in 0..SEASON {
            // Three seasons
            for m in 0..SSN {
                // Three sunspot ranges
                // Read the next four lines of text.
                current_line += 4;

                for k in (0..LAT).rev() {
                    // 19 latitudes counting backward to make the indices correspond to increasing latitude
                    if current_line < lines.len() {
                        let line = lines[current_line];
                        let parts: Vec<&str> = line.split("    ").collect();

                        if parts.len() >= 25 {
                            // Expecting at least 25 parts (1 label + 24 values)
                            for l in 0..24 {
                                if let Ok(value) = parts[l + 1].parse::<f64>() {
                                    path.fof2var[i][l][k][m][n] = value;
                                }
                            }
                        }

                        current_line += 1;
                    }
                }
            }
        }
    }

    true
}
