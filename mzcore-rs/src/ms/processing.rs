///
/// Some parts of this file originates from [Sage](https://github.com/lazear/sage/blob/master/crates/sage/src/spectrum.rs)
/// Copyright (c) 2022 Michael Lazear
/// SPDX-License-Identifier: MIT
///

use crate::ms::spectrum::Peak;
use crate::ms::utils::{binary_search_slice, MassTolWindow};

/// Binary search followed by linear search to select the most intense peak within `tolerance` window
/// * `offset` - this parameter allows for a static adjustment to the lower and upper bounds of the search window.
///     Sage subtracts a proton (and assumes z=1) for all experimental peaks, and stores all fragments as monoisotopic
///     masses. This simplifies downstream calculations at multiple charge states, but it also subtly changes tolerance
///     bounds. For most applications this is completely OK to ignore - however, for exact similarity of TMT reporter ion
///     measurements with ProteomeDiscoverer, FragPipe, etc, we need to account for this minor difference (which has an impact
///     perhaps 0.01% of the time)
pub fn select_most_intense_peak(
    peaks: &[Peak],
    center: f64,
    tolerance: MassTolWindow,
    offset: Option<f64>,
) -> Option<&Peak> {
    let (lo, hi) = tolerance.bounds(center);
    let (lo, hi) = (
        lo + offset.unwrap_or_default(),
        hi + offset.unwrap_or_default(),
    );

    let (i, j) = binary_search_slice(peaks, |peak, query| peak.mz.total_cmp(query), lo, hi);

    let mut best_peak = None;
    let mut max_int = 0.0;
    for peak in peaks[i..j]
        .iter()
        .filter(|peak| peak.mz >= lo && peak.mz <= hi)
    {
        if peak.intensity >= max_int {
            max_int = peak.intensity;
            best_peak = Some(peak);
        }
    }
    best_peak
}