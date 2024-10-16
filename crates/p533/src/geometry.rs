use crate::constants;
use crate::Location;

/// Determines the distance in km between here and there on the great circle
pub(crate) fn great_circle_distance(here: &Location, there: &Location) -> f64 {
    2.0 * 6371.0
        * (((here.latitude - there.latitude) / 2.0).sin().powi(2)
            + here.latitude.cos()
                * there.latitude.cos()
                * ((here.longitude - there.longitude) / 2.0).sin().powi(2))
        .sqrt()
        .asin()
}

pub(crate) struct CirclePoint {
    pub(crate) distance: f64,
    pub(crate) location: Location,
}

pub(crate) fn great_circle_point(
    here: &Location,
    there: &Location,
    distance: f64,
    fraction: f64,
) -> CirclePoint {
    match distance {
        dst if dst == 0.0 => CirclePoint {
            distance: 0.0,
            location: Location {
                latitude: here.latitude,
                longitude: here.longitude,
            },
        },
        _ => {
            let d = distance / constants::EARTH_RADIUS_0;
            let a = ((1. - fraction) * d).sin() / d.sin();
            let b = (fraction * d).sin() / d.sin();

            let x = a * (here.latitude).cos() * (here.longitude).cos()
                + b * (there.latitude).cos() * (there.longitude).cos();
            let y = a * (here.latitude).cos() * (here.longitude).sin()
                + b * (there.latitude).cos() * (there.longitude).sin();
            let z = a * (here.latitude).sin() + b * (there.latitude).sin();

            CirclePoint {
                distance: distance * fraction,
                location: Location {
                    latitude: z.atan2(x.powi(2) + y.powi(2)).sqrt(),
                    longitude: y.atan2(x),
                },
            }
        }
    }
}
