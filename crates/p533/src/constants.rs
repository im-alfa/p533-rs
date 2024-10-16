pub(crate) const EARTH_RADIUS_0: f64 = 6371.009;
pub(crate) const DBL_EPSILON: f64 = 2.2204460492503131E-016;

#[repr(usize)]
pub(crate) enum ControlPointIndexes {
    T1k = 0,
    Td02,
    MP,
    Rd02,
    R1k,
}
