use lazy_static::lazy_static;
use std::f64::consts::PI;

#[derive(Debug, Default, Clone, Copy)]
struct ManMadeNoiseParams {
    fa_m: f64, // Man-made noise
    du_m: f64, // Man-made noise upper decile
    dl_m: f64, // Man-made noise lower decile
}

#[derive(Debug, Default, Clone, Copy)]
struct GalacticNoiseParams {
    fa_g: f64, // Galactic noise
    du_g: f64, // Galactic noise upper decile
    dl_g: f64, // Galactic noise lower decile
}

#[derive(Debug, Default)]
pub struct AtmosphericStats {
    timeblock: u32, // Timeblock
    fa: f64,        // Atmospheric noise in dB above kT0b at 1 MHz
    sigma_fam: f64, // Standard deviation of values, Fam
    du: f64,        // Ratio of upper decile to median value, Fam
    sigma_du: f64,  // Standard deviations of values of Du
    dl: f64,        // Ratio of median value, Fam, to lower decile
    sigma_dl: f64,  // Standard deviation of values of Dl
}

#[derive(Debug, Default, Clone)]
pub struct AtmosphericCoefficients {
    fakp: Vec<Vec<Vec<f64>>>,
    fakabp: Vec<Vec<f64>>,
    fam: Vec<Vec<f64>>,
    dud: Vec<Vec<Vec<f64>>>,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct AtmosphericNoise {
    fa_a: f64, // Atmospheric noise
    du_a: f64, // Atmospheric noise upper decile
    dl_a: f64, // Atmospheric noise lower decile
}

#[derive(Debug, Default, Clone, Copy)]
pub struct NoiseParams {
    man_made_noise: ManMadeNoiseParams,
    galactic_noise: GalacticNoiseParams,
    atmospheric_noise: AtmosphericNoise,
    du_t: f64,  // Total noise upper decile
    dl_t: f64,  // Total noise lower decile
    fam_t: f64, // Total noise
}

#[derive(Debug, Default)]
pub enum TerrainCategory {
    #[default]
    City,
    Residential,
    Rural,
    Quietrural,
    Quiet,
    Noisy,
}

#[derive(Debug)]
pub struct FamStats {
    fa: f64,        // Atmospheric noise in dB above kT0b at 1 MHz
    sigma_fam: f64, // Standard deviation of values, Fam
    du: f64,        // Ratio of upper decile to median value, Fam
    sigma_du: f64,  // Standard deviations of values of Du
    dl: f64,        // Ratio of median value, Fam, to lower decile
    sigma_dl: f64,  // Standard deviation of values of Dl
}

lazy_static! {
    static ref ATMOSPHERIC_COEFFS: Vec<AtmosphericCoefficients> = (1..13)
        .map(|month| read_fam_dud(month))
        .collect::<Vec<AtmosphericCoefficients>>();
}

pub fn noise(
    hour: u32,
    rx_long: f64,
    rx_lat: f64,
    terrain_category: TerrainCategory,
    frequency: f64,
    month: u32,
) -> NoiseParams {
    let atmospheric_noise = calculate_atmospheric_noise(hour, rx_long, rx_lat, frequency, month);
    let man_made_noise = calculate_man_made_noise(terrain_category, frequency);
    let galactic_noise = calculate_galactic_noise(frequency);

    // Determine the combined noise according to
    // ITU-R P.372-10 Section 8 "The Combination of Noises from Several Sources"
    let find_decile = |sigma_a: f64, sigma_g: f64, sigma_m: f64| -> (f64, f64) {
        let c = 10.0 / 10.0_f64.log(std::f64::consts::E);

        let alpha_t = ((atmospheric_noise.fa_a / c) + (sigma_a.powi(2) / (c.powi(2) * 2.0))).exp()
            + ((galactic_noise.fa_g / c) + (sigma_g.powi(2) / (c.powi(2) * 2.0))).exp()
            + ((man_made_noise.fa_m / c) + (sigma_m.powi(2) / (c.powi(2) * 2.0))).exp();

        let beta_t = ((atmospheric_noise.fa_a / c) + sigma_a.powi(2) / (c.powi(2) * 2.0))
            .exp()
            .powi(2)
            * ((sigma_a / c).powi(2).exp() - 1.0)
            + ((galactic_noise.fa_g / c) + sigma_g.powi(2) / (c.powi(2) * 2.0))
                .exp()
                .powi(2)
                * (((sigma_g / c).powi(2)).exp() - 1.0)
            + ((man_made_noise.fa_m / c) + sigma_m.powi(2) / (c.powi(2) * 2.0))
                .exp()
                .powi(2)
                * (((sigma_m / c).powi(2)).exp() - 1.0);

        let gamma_t = (atmospheric_noise.fa_a / c).exp()
            + (galactic_noise.fa_g / c).exp()
            + (man_made_noise.fa_m / c).exp();

        let sigma_t = if (atmospheric_noise.du_a > 12.0)
            || (galactic_noise.du_g > 12.0)
            || (man_made_noise.du_m > 12.0)
        {
            c * (2.0 * (alpha_t / gamma_t).log(std::f64::consts::E)).sqrt()
        } else {
            c * ((1.0 + beta_t / alpha_t.powi(2)).log(std::f64::consts::E)).sqrt()
        };

        let fam_t = c * (alpha_t.log(std::f64::consts::E) - (sigma_t.powi(2) / (c.powi(2) * 2.0)));

        let d_t = 1.282 * sigma_t;

        (fam_t, d_t)
    };

    // Find the upper decile sigma_t
    let sigma_a = atmospheric_noise.du_a / 1.282;
    let sigma_g = 1.56_f64;
    let sigma_m = man_made_noise.du_m / 1.282;

    let (fam_tu, du_t) = find_decile(sigma_a, sigma_g, sigma_m);

    // Find the lower decile sigma_t
    let sigma_a = atmospheric_noise.dl_a / 1.282;
    let sigma_g = 1.56;
    let sigma_m = man_made_noise.dl_m / 1.282;

    let (fam_tl, dl_t) = find_decile(sigma_a, sigma_g, sigma_m);

    let fam_t = f64::min(fam_tl, fam_tu);

    NoiseParams {
        atmospheric_noise,
        man_made_noise,
        galactic_noise,
        du_t,
        dl_t,
        fam_t,
    }
}

pub fn calculate_atmospheric_noise(
    hour: u32,
    rx_lng: f64,
    rx_lat: f64,
    frequency: f64,
    month: u32,
) -> AtmosphericNoise {
    // calculate local receiver mean time
    // every 15 degrees of longitude, the local mean time is 1 hour different
    // https://www.cs4fn.org/mobile/owntimezone.php
    let mut local_receiver_mean_time = hour as i32 + (rx_lng / (15.0 * PI / 180.0)) as i32;

    if local_receiver_mean_time < 0 {
        local_receiver_mean_time += 24;
    } else if local_receiver_mean_time > 23 {
        local_receiver_mean_time -= 24;
    }

    // The atmospheric noise is determined by
    // i) finding the atmospheric noise at the current time
    // block at the receiver local mean time,
    //
    // ii) finding the noise for the "adjacent" time block and
    // then
    //
    // iii) iterating between the two noise values. There are 6 time blocks so the modulo 6
    // keeps the indexes inbounds
    let fs_now_timeblock = (local_receiver_mean_time / 4) % 6;
    let fs_now = get_fam_parameters(rx_lng, rx_lat, frequency, fs_now_timeblock as usize, month);

    let fs_adj_timeblock = (fs_now_timeblock + 1) % 6;
    let fs_adj = get_fam_parameters(rx_lng, rx_lat, frequency, fs_adj_timeblock as usize, month);

    // Interpolate is based on the local receiver mean time, lrxmt,
    // and the 4 hour timeblock.
    let slp = (local_receiver_mean_time % 4) as f64 / 4.0;

    AtmosphericNoise {
        fa_a: 10.0
            * (10.0_f64.powf(fs_now.fa / 10.0)
                + (10.0_f64.powf(fs_adj.fa / 10.0) - 10.0_f64.powf(fs_now.fa / 10.0)) * slp)
                .log10(),
        du_a: 10.0
            * (10.0_f64.powf(fs_now.du / 10.0)
                + (10.0_f64.powf(fs_adj.du / 10.0) - 10.0_f64.powf(fs_now.du / 10.0)) * slp)
                .log10(),
        dl_a: 10.0
            * (10.0_f64.powf(fs_now.dl / 10.0)
                + (10.0_f64.powf(fs_adj.dl / 10.0) - 10.0_f64.powf(fs_now.dl / 10.0)) * slp)
                .log10(),
    }
}

pub fn get_fam_parameters(lng: f64, lat: f64, frequency: f64, timeblock: usize, month: u32) -> FamStats {
    let month = (month - 1) as usize;

    // Find the atmospheric noise Fam (db above kT0b at 1 MHz)

    // limits for Fourier series
    let lm = 29;
    let ln = 15;

    // the longitude used here is the geographic east long (0 to 2*PI rads)
    // initialize the tepm, q, as half the graphic east long

    let q = if lng < 0.0 {
        (lng + 2.0 * PI) / 2.0
    } else {
        lng / 2.0
    };

    // calculate longitude series
    let mut zz = [0.0; 30];
    for j in 0..lm {
        let mut r = 0.0;
        for k in 0..ln {
            r += ((k as f64 + 1.0) * q).sin() * ATMOSPHERIC_COEFFS[month].fakp[j][k][timeblock];
        }

        zz[j] = r + ATMOSPHERIC_COEFFS[month].fakp[j][15][timeblock];
    }

    // Calculate the latitude series
    // Reuse the temp, q, as the latitude plus 90 degrees
    let q = lat + PI / 2.0;
    let mut r = 0.0;

    for j in 0..lm {
        r += ((j as f64 + 1.0) * q).sin() * zz[j];
    }

    // Final Fourier series calculation (Note the linear nomalization using fakabp values)
    let fam_1mhz = r
        + ATMOSPHERIC_COEFFS[month].fakabp[0][timeblock]
        + ATMOSPHERIC_COEFFS[month].fakabp[1][timeblock] * q;

    // Determine if the receiver latitude is in the northern or southern hemisphere
    // If lat < 0, then timeblock += 6
    let timeblock = if lat < 0.0 { timeblock + 6 } else { timeblock };

    // for K = 0 then U1 = -0.75
    // for K = 1 then U1 = U
    let mut u = vec![0.0; 2];
    u[0] = -0.75;
    u[1] = (8.0 * 2.0_f64.powf(frequency.log10()) - 11.0) / 4.0; // U = (8. * 2.**X - 11.)/4. where X = ALOG10(FREQ)

    // Please See Page 5
    // NBS Tech Note 318 Lucas and Harper
    // "A Numerical Representation of CCIR Report 322 High Frequeny (3-30 Mc/s) Atmospheric Radio Noise Data"
    let mut cz = 0.0;
    let mut pz = 0.0;
    let mut px = 0.0;

    for k in 0..2 {
        pz = u[k] * ATMOSPHERIC_COEFFS[month].fam[0][timeblock]
            + ATMOSPHERIC_COEFFS[month].fam[1][timeblock];
        px = u[k] * ATMOSPHERIC_COEFFS[month].fam[7][timeblock]
            + ATMOSPHERIC_COEFFS[month].fam[8][timeblock];

        for j in 2..7 {
            pz = u[k] * pz + ATMOSPHERIC_COEFFS[month].fam[j][timeblock];
            px = u[k] * px + ATMOSPHERIC_COEFFS[month].fam[j + 7][timeblock];
        }

        if k == 0 {
            cz = fam_1mhz * (2.0 - pz) - px;
        }
    }

    // Frequency variation of atmospheric noise
    let fa = cz * pz + px;

    //	Limit frequency to 20 MHz for Du, Dl, SigmaDu, SigmaDl
    //	because curves in ITU-R P.372 only go to 20 MHz
    let mut x = if frequency > 20.0 {
        20.0_f64.log10()
    } else {
        frequency.log10()
    };
    let mut v = vec![0.0; 5];

    for j in 0..5 {
        // Limit frequency to 10 MHz for SigmaFam
        // because curves in ITU-R P.372 only go to 10 MHz
        if j == 4 && frequency > 10.0 {
            x = 1.0;
        }

        let mut y = ATMOSPHERIC_COEFFS[month].dud[0][timeblock][j]; // Y = DUD(1,TIMEBLOCKINDX,I)

        for k in 1..5 {
            y = x * y + ATMOSPHERIC_COEFFS[month].dud[k][timeblock][j]; // Y = Y*X + DUD(J,TIMEBLOCKINDX,I)
        }

        v[j] = y;
    }

    FamStats {
        fa,
        du: v[0],
        dl: v[1],
        sigma_du: v[2],
        sigma_dl: v[3],
        sigma_fam: v[4],
    }
}

const CITY: (f64, f64, f64, f64) = (76.8, 27.7, 11.0, 6.7);
const RESIDENTIAL: (f64, f64, f64, f64) = (72.5, 27.7, 10.6, 5.3);
const RURAL: (f64, f64, f64, f64) = (67.2, 27.7, 9.2, 4.6);
const QUIET_RURAL: (f64, f64, f64, f64) = (53.6, 28.6, 9.2, 4.6);
const QUIET: (f64, f64, f64, f64) = (65.2, 29.1, 9.2, 4.6);
const NOISY: (f64, f64, f64, f64) = (83.2, 37.5, 11.0, 6.7);

fn calculate_man_made_noise(
    terrain_category: TerrainCategory,
    frequency: f64,
) -> ManMadeNoiseParams {
    let (c, d, du_m, dl_m) = match terrain_category {
        TerrainCategory::City => CITY,
        TerrainCategory::Residential => RESIDENTIAL,
        TerrainCategory::Rural => RURAL,
        TerrainCategory::Quietrural => QUIET_RURAL,
        TerrainCategory::Quiet => QUIET,
        TerrainCategory::Noisy => NOISY,
    };

    let fa_m = c - d * frequency.log10();

    ManMadeNoiseParams { fa_m, du_m, dl_m }
}

fn calculate_galactic_noise(frequency: f64) -> GalacticNoiseParams {
    let c = 52.0;
    let d = 23.0;

    return GalacticNoiseParams {
        du_g: 2.0,
        dl_g: 2.0,
        fa_g: c - d * frequency.log10(),
    };
}

pub fn read_fam_dud(month: u32) -> AtmosphericCoefficients {
    let lines = match month {
        1 => include_str!("../static_data/COEFF01W.txt"),
        2 => include_str!("../static_data/COEFF02W.txt"),
        3 => include_str!("../static_data/COEFF03W.txt"),
        4 => include_str!("../static_data/COEFF04W.txt"),
        5 => include_str!("../static_data/COEFF05W.txt"),
        6 => include_str!("../static_data/COEFF06W.txt"),
        7 => include_str!("../static_data/COEFF07W.txt"),
        8 => include_str!("../static_data/COEFF08W.txt"),
        9 => include_str!("../static_data/COEFF09W.txt"),
        10 => include_str!("../static_data/COEFF10W.txt"),
        11 => include_str!("../static_data/COEFF11W.txt"),
        12 => include_str!("../static_data/COEFF12W.txt"),
        _ => "",
    };

    // Skip lines
    // * Skip first header line. -> 1 lines
    // * Skip if2(10) & xf2(13,76,2) -> 400 lines
    // * Skip ifm3(10) & xfm3(9,49,2) -> 181 lines
    // * Skip ie(10) & xe(9,22,2) -> 84 lines
    // * Skip iesu(10) & xesu(5,55,2) -> 114 lines
    // * Skip ies(10) & xes(7,61,2) -> 175 lines
    // * Skip iels(10) & xels(5,55,2 -> 114 lines
    // * Skip ihpo1(10) & xhpo1(13,29,2) -> 155 lines
    // * Skip ihpo2(10) & xhpo2(9,55,2) -> 202 lines
    // * Skip ihp(10) & xhp(9,37,2) -> 138 lines
    let skip_lines = 1 + 400 + 181 + 84 + 114 + 175 + 114 + 155 + 202 + 138;

    let mut lines = lines.lines().skip(skip_lines);

    // ***
    // READ FAKP BLOCK
    // ***

    // skip header line
    lines.next();

    // Read fakp(29,16,6) from the table
    // 29 * 16 * 6 = 2784
    let mut a = vec![0.0; 2784];

    // Read the next 556 lines (the length of the fakp block)
    for n in 0..556 {
        let line = lines.next().unwrap();
        for x in 0..5 {
            let s = line.get(16 * x..16 * x + 16).unwrap();
            a[5 * n + x] = s.trim().parse::<f64>().unwrap();
        }
    }

    // Read the last partial line
    // It's the same as the previous loop but with 4 blocks instead of 5
    let line = lines.next().unwrap();
    for x in 0..4 {
        let s = line.get(16 * x..16 * x + 16).unwrap();
        a[5 * 556 + x] = s.trim().parse::<f64>().unwrap();
    }

    let mut fakp = vec![vec![vec![0.0; 6]; 16]; 29];

    // Reshape A into the Coeff structure
    for i in 0..6 {
        for j in 0..16 {
            for k in 0..29 {
                fakp[k][j][i] = a[16 * 29 * i + 29 * j + k];
            }
        }
    }

    // ***
    // READ FAKABP BLOCK
    // ***

    // Read fakabp(2,6) from the table
    // 2 * 6 = 12
    let mut a = vec![0.0; 12];

    // Skip header line
    lines.next();

    for n in 0..2 {
        let line = lines.next().unwrap();
        for x in 0..5 {
            let s = line.get(16 * x..16 * x + 16).unwrap();
            a[5 * n + x] = s.trim().parse::<f64>().unwrap();
        }
    }

    // read the last partial line
    // It's the same as the previous loop but with 2 blocks instead of 5
    let line = lines.next().unwrap();
    for x in 0..2 {
        let s = line.get(16 * x..16 * x + 16).unwrap();
        a[5 * 2 + x] = s.trim().parse::<f64>().unwrap();
    }

    // Reshape A into the Coeff structure
    let mut fakabp = vec![vec![0.0; 6]; 2];
    for j in 0..6 {
        for k in 0..2 {
            fakabp[k][j] = a[2 * j + k];
        }
    }

    // ***
    // READ DUD BLOCK
    // ***

    // Read dud(5,12,5) from the table
    // 5 * 12 * 5 = 300
    let mut a = vec![0.0; 300];

    // Skip header line
    lines.next();

    // Read 60 lines
    for n in 0..60 {
        let line = lines.next().unwrap();
        for x in 0..5 {
            let s = line.get(16 * x..16 * x + 16).unwrap();
            a[5 * n + x] = s.trim().parse::<f64>().unwrap();
        }
    }

    // Reshape A into the Coeff structure
    let mut dud = vec![vec![vec![0.0; 5]; 12]; 5];
    for i in 0..5 {
        for j in 0..12 {
            for k in 0..5 {
                dud[k][j][i] = a[12 * 5 * i + 5 * j + k];
            }
        }
    }

    // ***
    // READ FAM BLOCK
    // ***

    // Read fam(14,12) from the table
    // 14 * 12 = 168
    let mut a = vec![0.0; 168];

    // Skip header line
    lines.next();

    // Read 33 lines
    for n in 0..33 {
        let line = lines.next().unwrap();
        for x in 0..5 {
            let s = line.get(16 * x..16 * x + 16).unwrap();
            a[5 * n + x] = s.trim().parse::<f64>().unwrap();
        }
    }

    // Read the last partial line
    // It's the same as the previous loop but with 3 blocks instead of 5
    let line = lines.next().unwrap();
    for x in 0..3 {
        let s = line.get(16 * x..16 * x + 16).unwrap();
        a[5 * 33 + x] = s.trim().parse::<f64>().unwrap();
    }

    // Reshape A into the Coeff structure
    let mut fam = vec![vec![0.0; 12]; 14];
    for j in 0..12 {
        for k in 0..14 {
            fam[k][j] = a[14 * j + k];
        }
    }

    AtmosphericCoefficients {
        fakp,
        fakabp,
        fam,
        dud,
    }
}

impl NoiseParams {
    /// Get the total noise value
    pub fn get_total_noise(&self) -> f64 {
        self.fam_t
    }

    /// Get the total noise upper decile
    pub fn get_upper_decile(&self) -> f64 {
        self.du_t
    }

    /// Get the total noise lower decile
    pub fn get_lower_decile(&self) -> f64 {
        self.dl_t
    }

    /// Get the atmospheric noise component (FaA)
    pub fn get_atmospheric_noise(&self) -> f64 {
        self.atmospheric_noise.fa_a
    }

    /// Get the man-made noise component (FaM)
    pub fn get_man_made_noise(&self) -> f64 {
        self.man_made_noise.fa_m
    }

    /// Get the galactic noise component (FaG)
    pub fn get_galactic_noise(&self) -> f64 {
        self.galactic_noise.fa_g
    }

    /// Get the atmospheric noise deciles
    pub fn get_atmospheric_deciles(&self) -> (f64, f64) {
        (self.atmospheric_noise.du_a, self.atmospheric_noise.dl_a)
    }

    /// Get the man-made noise deciles  
    pub fn get_man_made_deciles(&self) -> (f64, f64) {
        (self.man_made_noise.du_m, self.man_made_noise.dl_m)
    }

    /// Get the galactic noise deciles
    pub fn get_galactic_deciles(&self) -> (f64, f64) {
        (self.galactic_noise.du_g, self.galactic_noise.dl_g)
    }
}
