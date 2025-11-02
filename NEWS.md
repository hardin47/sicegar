# sicegar 0.3.0

* Added functionality so that the lower asymptote (h0) can be estimated freely, instead of forced to be zero.
To implement, toggle `use_h0 = TRUE` in `fitAndCategorize()` and `figureModelCurves()`.
* Set default value of `threshold_t0_max_int` to 1E10.
* Fixed typo in object assignment of `decisionList$test.sm_tmax_IntensityRatio`.
* Fixed typo so that now the asymptote ratio (`parameterDF$finalAsymptoteIntensityRatio_Estimate`) is also unnormalized, along with all other parameters.
* Updated **ggplot** functions to have current syntax.
* Edited grammar and spelling throughout.
