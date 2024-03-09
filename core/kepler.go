package core

import "math"

const _DLA_DELTA = 1e-7 // precision for Kepler equation

// Solve Kepler equation to calculate ea, the eccentric anomaly,
// in elliptical motion given s (< 1), the eccentricity, and m, mean anomaly.
// All agular values are in radians.
func EccentricAnomaly(s, m, ea float64) float64 {
	dla := ea - (s * math.Sin(ea)) - m
	if math.Abs(dla) < _DLA_DELTA {
		return ea
	}

	dla = dla / (1 - (s * math.Cos(ea)))
	return EccentricAnomaly(s, m, ea-dla)
}

// Given s, eccentricity, and ea, eccentric anomaly, find true anomaly.
// All angular values are in radians.
func TrueAnomaly(s, ea float64) float64 {
	return 2 * math.Atan(math.Sqrt((1+s)/(1-s))*math.Tan(ea/2))
}
