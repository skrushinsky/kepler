// The main purpose is to convert between civil dates and Julian dates.
// Julian date (JD) is the number of days elapsed since mean UT noon of
// January 1st 4713 BC. This system of time measurement is widely adopted by
// the astronomers.
package sun

import (
	"math"

	"github.com/skrushinsky/kepler/core"
	"github.com/skrushinsky/scaliger/julian"
	"github.com/skrushinsky/scaliger/mathutils"
)

const ABERRATION = 5.69e-3 // aberration in degrees
const _PI2 = math.Pi * 2

// Controls  type of the result.
type SunOptions struct {
	// nutation in longitude, degrees
	dpsi float64
	// ignore light-time travel correction?
	ignoreLightTravel bool
	// Mean Longitude of the Sun, degrees
	meanLongitude float64
	// Mean Anomaly of the Sun, degrees
	meanAnomaly float64
}

// Mean longitude of the Sun, arc-degrees
func MeanLongitude(t float64) float64 {
	return mathutils.ReduceDeg(2.7969668e2 + 3.025e-4*t*t + mathutils.Frac360(1.000021359e2*t))
}

// Mean anomaly of the Sun, arc-degrees
func MeanAnomaly(t float64) float64 {
	return mathutils.ReduceDeg(3.5847583e2 - (1.5e-4+3.3e-6*t)*t*t + mathutils.Frac360(9.999736042e1*t))
}

func TrueGeocentric(t, ms, ls float64) (lsn float64, rsn float64) {
	ma := mathutils.Radians(ms)
	s := mathutils.Polynome(t, 1.675104e-2, -4.18e-5, -1.26e-7) // eccentricity
	ea := core.EccentricAnomaly(s, ma, ma)                      // eccentric anomaly
	nu := core.TrueAnomaly(s, ea)                               // true anomaly
	t2 := t * t

	calcPert := func(a, b float64) float64 {
		return mathutils.Radians(a + mathutils.Frac360(b*t))
	}
	a := calcPert(153.23, 6.255209472e1)            // Venus
	b := calcPert(216.57, 1.251041894e2)            // ?
	c := calcPert(312.69, 9.156766028e1)            // ?
	d := calcPert(350.74-1.44e-3*t2, 1.236853095e3) // Moon
	h := calcPert(353.4, 1.831353208e2)             // ?
	e := mathutils.Radians(231.19 + 20.2*t)         // inequality of long period

	// correction in orbital longitude
	dl := 1.34e-3*math.Cos(a) +
		1.54e-3*math.Cos(b) +
		2e-3*math.Cos(c) +
		1.79e-3*math.Sin(d) +
		1.78e-3*math.Sin(e)
	// correction in radius-vector
	dr := 5.43e-6*math.Sin(a) +
		1.575e-5*math.Sin(b) +
		1.627e-5*math.Sin(c) +
		3.076e-5*math.Cos(d) +
		9.27e-6*math.Sin(h)
	lsn = mathutils.ReduceDeg(mathutils.Degrees(nu) + ls - ms + dl)
	rsn = 1.0000002*(1-s*math.Cos(ea)) + dr
	return
}

func Apparent(jd float64, options SunOptions) core.EclipticPosition {
	t := (jd - julian.J1900) / julian.DAYS_PER_CENT
	// ms := MeanAnomaly(t)
	// ls := MeanLongitude(t)
	lsn, rsn := TrueGeocentric(t, options.meanAnomaly, options.meanLongitude)
	lsn += options.dpsi // correct for nutation
	lsn -= ABERRATION   // correct for aberration
	if !options.ignoreLightTravel {
		dt := 1.365 * rsn     // seconds
		lsn -= dt * 15 / 3600 // convert to degrees and substract
	}
	return core.EclipticPosition{Lambda: lsn, Delta: rsn}
}
