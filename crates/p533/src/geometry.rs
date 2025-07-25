use crate::constants::*;
use crate::path_data::*;

/// Determines the distance in km between here and there on the great circle.
pub fn great_circle_distance(here: Location, there: Location) -> f64 {
    2.0 * R0
        * f64::asin(f64::sqrt(
            f64::powi(f64::sin((here.lat - there.lat) / 2.0), 2)
                + f64::cos(here.lat)
                    * f64::cos(there.lat)
                    * f64::powi(f64::sin((here.lng - there.lng) / 2.0), 2),
        ))
}

/// Determines the lat and long of a point on a great circle path of length distance from here to there.
/// The point is determined at fraction*distance from here.
pub fn great_circle_point(
    here: Location,
    there: Location,
    midpnt: &mut ControlPt,
    distance: f64,
    fraction: f64,
) {
    if distance != 0.0 {
        midpnt.distance = distance * fraction;
        let d = distance / R0;
        let a = f64::sin((1.0 - fraction) * d) / f64::sin(d);
        let b = f64::sin(fraction * d) / f64::sin(d);
        let x = a * f64::cos(here.lat) * f64::cos(here.lng)
            + b * f64::cos(there.lat) * f64::cos(there.lng);
        let y = a * f64::cos(here.lat) * f64::sin(here.lng)
            + b * f64::cos(there.lat) * f64::sin(there.lng);
        let z = a * f64::sin(here.lat) + b * f64::sin(there.lat);
        midpnt.l.lat = f64::atan2(z, f64::sqrt(x.powi(2) + y.powi(2)));
        midpnt.l.lng = f64::atan2(y, x);
    } else {
        midpnt.distance = 0.0;
        midpnt.l.lat = here.lat;
        midpnt.l.lng = here.lng;
    }
}

/// Conversion from geographic coordinates to geomagnetic coordinates (lat, lng)
pub fn geomagnetic_coords(here: Location, there: &mut Location) {
    // This was the pole location when the coefficients Lh were determined.
    // In 1955 the Geomagnetic North Pole was at 78.5N   69.2W,
    //              while the South Pole was at  78.5S  110.8E
    let geo_mag_n_pole = Location {
        lat: 78.5 * D2R,
        lng: -68.2 * D2R,
    };

    there.lat = f64::asin(
        f64::sin(here.lat) * f64::sin(geo_mag_n_pole.lat)
            + f64::cos(here.lat)
                * f64::cos(geo_mag_n_pole.lat)
                * f64::cos(here.lng - geo_mag_n_pole.lng),
    );
    there.lng = f64::asin(
        f64::cos(here.lat) * f64::sin(here.lng - geo_mag_n_pole.lng) / f64::cos(there.lat),
    );
}

/// Determines the bearing from here to there
pub fn bearing(here: Location, there: Location) -> f64 {
    // N and E are positive
    // S and W are negative

    let numerator = f64::sin(there.lng - here.lng) * f64::cos(there.lat);
    let denominator = f64::cos(here.lat) * f64::sin(there.lat)
        - f64::sin(here.lat) * f64::cos(there.lat) * f64::cos(there.lng - here.lng);

    let mut bearing = f64::atan2(numerator, denominator);

    bearing = (2.0 * PI + bearing) % (2.0 * PI);

    bearing
}
