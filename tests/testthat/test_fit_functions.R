test_that('sigmoidal fit function', {
  time <- seq(1, 50, 1)
  intensity <- sigmoidalFitFormula(time, maximum = 4, slopeParam = 1, midPoint = 8)
  expect_equal(max(intensity), 4) # maximum is 4
  expect_equal(intensity[8], 2) # midpoint is at time = 8
})

test_that('double-sigmoidal fit function', {
  time <- seq(1, 50, 1)
  intensity <- doublesigmoidalFitFormula(time,
                                         finalAsymptoteIntensityRatio = .3,
                                         maximum = 4,
                                         slope1Param = 1,
                                         midPoint1Param = 7,
                                         slope2Param = 1,
                                         midPointDistanceParam = 8)
  expect_equal(max(intensity), 4) # maximum is 4
  expect_equal(intensity[length(intensity)]/max(intensity), .3) # ratio from max to final is 0.3
})

test_that('sigmoidal fit function with free h0', {
  time <- seq(1, 50, 1)
  intensity <- sigmoidalFitFormula_h0(time, maximum = 4, slopeParam = 1, midPoint = 8, h0 = 1)
  expect_equal(max(intensity), 4) # maximum is 4
  expect_equal(intensity[8], 2) # midpoint is at time = 8
})

test_that('double-sigmoidal fit function with free h0', {
  time <- seq(1, 50, 1)
  intensity <- doublesigmoidalFitFormula_h0(time,
                                         finalAsymptoteIntensityRatio = .3,
                                         maximum = 4,
                                         slope1Param = 1,
                                         midPoint1Param = 7,
                                         slope2Param = 1,
                                         midPointDistanceParam = 8,
                                         h0 = 1)
  expect_equal(max(intensity), 4) # maximum is 4
  expect_equal(intensity[length(intensity)]/max(intensity), .3) # ratio from max to final is 0.3
})
