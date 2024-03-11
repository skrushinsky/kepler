package moon

import (
	"math"

	"github.com/skrushinsky/kepler/core"
	"github.com/skrushinsky/scaliger/julian"
	"github.com/skrushinsky/scaliger/mathutils"
)

var MoonOrbit = map[string][]float64{
	// Mean longitude
	"L": {
		218.3164477,
		481267.88123421,
		-0.0015786,
		1.0 / 538841,
		-(1.0 / 65194000),
	},
	// Mean elongation
	"D": {
		297.8501921,
		445267.1114034,
		-0.0018819,
		1.0 / 545868,
		-(1.0 / 113065000),
	},
	// Mean anomaly
	"M": {134.9633964, 477198.8675055, 0.0087414, 1.0 / 69699, -(1.0 / 14712000)},
	// Argument of latitude (mean distance of the Moon from its ascending node)
	"F": {
		93.272095,
		483202.0175233,
		-0.0036539,
		-(1.0 / 3526000),
		1.0 / 863310000,
	},
}

var _M = [...]float64{27.32158213, 365.2596407, 27.55455094, 29.53058868, 27.21222039, 6798.363307}

// Mean anomaly of the Sun
var SunOrbit = map[string][]float64{
	"M": {357.5291092, 35999.0502909, -0.0001536, 1.0 / 24490000},
}

var sin = math.Sin
var cos = math.Cos
var radians = mathutils.Radians
var polynome = mathutils.Polynome
var reduceDeg = mathutils.ReduceDeg

// Mean Lunar Node.
// t is a number of Julian centuries elapsed since 1900, Jan 0.5.
// Returns degrees.
func MeanLunarNode(t float64) float64 {
	return reduceDeg(polynome(t, 125.0445479, -1934.1362891, 0.0020754, 1.0/467441, 1.0/60616000))
}

// Longitude of Lunar Node, arc-degrees.
// jd is a Standard Julian Date.
// If mean is false (the default), the result refers to the true equinox of the date, otherwise
// to the mean equinox of date, so the result is equal to [MeanLunarNode].
func LunarNode(jd float64, mean bool) float64 {
	t := (jd - julian.J2000) / julian.DAYS_PER_CENT
	mn := polynome(t, 125.0445479, -1934.1362891, 0.0020754, 1.0/467441, 1.0/60616000)
	var nd float64
	if mean {
		nd = mn
	} else {
		assemble := func(terms ...float64) float64 { return radians(reduceDeg(polynome(t, terms...))) }
		assembleMoon := func(k string) float64 { return assemble(MoonOrbit[k]...) }
		d := assembleMoon("D")
		m := assembleMoon("M")
		f := assembleMoon("F")
		ms := assemble(SunOrbit["M"]...)
		nd = mn - 1.4979*sin(2*(d-f)) -
			0.1500*sin(ms) -
			0.1226*sin(2*d) +
			0.1176*sin(2*f) -
			0.0801*sin(2*(m-f))
	}
	return reduceDeg(nd)
}

// True position of the Moon.
// Given Julian Day, calculates Moon position, horizontal parallax (A.U.) and angular speed, degrees / 24h.
func TruePosition(jd float64) (pos core.EclipticPosition, parallax, motion float64) {
	djd := jd - julian.J1900
	t := djd / julian.DAYS_PER_CENT
	t2 := t * t
	m := core.Map(_M[:], func(x float64) float64 { return mathutils.Frac360(djd / x) })
	ld := 270.434164 + m[0] - (1.133e-3-1.9e-6*t)*t2  // Moon's mean longitude
	ms := 358.475833 + m[1] - (1.5e-4+3.3e-6*t)*t2    // mean anomaly of the Sun
	md := 296.104608 + m[2] + (9.192e-3+1.44e-5*t)*t2 // mean anomaly
	de := 350.737486 + m[3] - (1.436e-3-1.9e-6*t)*t2  // mean elongation
	f := 11.250889 + m[4] - (3.211e-3+3e-7*t)*t2      // mean distance of Moon from its ascending node
	n := 259.183275 - m[5] + (2.078e-3+2.2e-5*t)*t2   // longitude of Moon's asc. node
	a := radians(51.2 + 20.2*t)
	sa := sin(a)
	sn := sin(radians(n))
	b := 346.56 + (132.87-9.1731e-3*t)*t
	sb := 3.964e-3 * sin(radians(b))
	c := radians(n + 275.05 - 2.3*t)
	sc := sin(c)
	ld += 2.33e-4*sa + sb + 1.964e-3*sn
	ms -= 1.778e-3 * sa
	md += 8.17e-4*sa + sb + 2.541e-3*sn
	f += sb - 2.4691e-2*sn - 4.328e-3*sc
	de += 2.011e-3*sa + sb + 1.964e-3*sn
	e := 1 - (2.495e-3+7.52e-6*t)*t
	e2 := e * e
	ms = radians(ms)
	n = radians(n)
	de = radians(de)
	f = radians(f)
	md = radians(md)

	de2 := de + de
	de3 := de2 + de
	de4 := de2 + de2
	md2 := md + md
	md3 := md2 + md
	ms2 := ms + ms
	f2 := f + f
	f3 := f2 + f
	// ecliptic longitude
	l := 6.28875*sin(md) +
		1.274018*sin(de2-md) +
		6.58309e-1*sin(de2) +
		2.13616e-1*sin(md2) -
		e*1.85596e-1*sin(ms) -
		1.14336e-1*sin(f2) +
		5.8793e-2*sin(2*(de-md)) +
		5.7212e-2*e*sin(de2-ms-md) +
		5.332e-2*sin(de2+md) +
		4.5874e-2*e*sin(de2-ms) +
		4.1024e-2*e*sin(md-ms) -
		3.4718e-2*sin(de) -
		e*3.0465e-2*sin(ms+md) +
		1.5326e-2*sin(2*(de-f)) -
		1.2528e-2*sin(f2+md) -
		1.098e-2*sin(f2-md) +
		1.0674e-2*sin(de4-md) +
		1.0034e-2*sin(md3) +
		8.548e-3*sin(de4-md2) -
		e*7.91e-3*sin(ms-md+de2) -
		e*6.783e-3*sin(de2+ms) +
		5.162e-3*sin(md-de) +
		e*5e-3*sin(ms+de) +
		3.862e-3*sin(de4) +
		e*4.049e-3*sin(md-ms+de2) +
		3.996e-3*sin(2*(md+de)) +
		3.665e-3*sin(de2-md3) +
		e*2.695e-3*sin(md2-ms) +
		2.602e-3*sin(md-2*(f+de)) +
		e*2.396e-3*sin(2*(de-md)-ms) -
		2.349e-3*sin(md+de) +
		e2*2.249e-3*sin(2*(de-ms)) -
		e*2.125e-3*sin(md2+ms) -
		e2*2.079e-3*sin(ms2) +
		e2*2.059e-3*sin(2*(de-ms)-md) -
		1.773e-3*sin(md+2*(de-f)) -
		1.595e-3*sin(2*(f+de)) +
		e*1.22e-3*sin(de4-ms-md) -
		1.11e-3*sin(2*(md+f)) +
		8.92e-4*sin(md-de3) -
		e*8.11e-4*sin(ms+md+de2) +
		e*7.61e-4*sin(de4-ms-md2) +
		e2*7.04e-4*sin(md-2*(ms+de)) +
		e*6.93e-4*sin(ms-2*(md-de)) +
		e*5.98e-4*sin(2*(de-f)-ms) +
		5.5e-4*sin(md+de4) +
		5.38e-4*sin(4*md) +
		e*5.21e-4*sin(de4-ms) +
		4.86e-4*sin(md2-de) +
		e2*7.17e-4*sin(md-ms2)

	pos.Lambda = mathutils.ReduceDeg(ld + l)

	// ecliptic latitude
	g := 5.128189*sin(f) +
		.280606*sin(md+f) +
		.277693*sin(md-f) +
		.173238*sin(de2-f) +
		.055413*sin(de2+f-md) +
		.046272*sin(de2-f-md) +
		.032573*sin(de2+f) +
		.017198*sin(md2+f) +
		.009267*sin(de2+md-f) +
		.008823*sin(md2-f) +
		e*.008247*sin(de2-ms-f) +
		.004323*sin(2*(de-md)-f) +
		.0042*sin(de2+f+md) +
		e*.003372*sin(f-ms-de2) +
		e*.002472*sin(de2+f-ms-md) +
		e*.002222*sin(de2+f-ms) +
		e*.002072*sin(de2-f-ms-md) +
		e*.001877*sin(f-ms+md) +
		.001828*sin(de4-f-md) -
		e*.001803*sin(f+ms) -
		.00175*sin(f3) +
		e*.00157*sin(md-ms-f) -
		.001487*sin(f+de) -
		e*.001481*sin(f+ms+md) +
		e*.001417*sin(f-ms-md) +
		e*.00135*sin(f-ms) +
		.00133*sin(f-de) +
		.001106*sin(f+md3) +
		.00102*sin(de4-f) +
		.000833*sin(f+de4-md) +
		.000781*sin(md-f3) +
		.00067*sin(f+de4-md2) +
		.000606*sin(de2-f3) +
		.000597*sin(2*(de+md)-f) +
		e*.000492*sin(de2+md-ms-f) +
		.00045*sin(2*(md-de)-f) +
		.000439*sin(md3-f) +
		.000423*sin(f+2*(de+md)) +
		.000422*sin(de2-f-md3) -
		e*.000367*sin(ms+f+de2-md) -
		e*.000353*sin(ms+f+de2) +
		.000331*sin(f+de4) +
		e*.000317*sin(de2+f-ms+md) +
		e2*.000306*sin(2*(de-ms)-f) -
		.000283*sin(md+f3)
	w1 := .0004664 * cos(n)
	w2 := .0000754 * cos(c)
	pos.Beta = g * (1 - w1 - w2)

	// horizontal parallax
	parallax = .950724 +
		.051818*cos(md) +
		.009531*cos(de2-md) +
		.007843*cos(de2) +
		.002824*cos(md2) +
		.000857*cos(de2+md) +
		e*.000533*cos(de2-ms) +
		e*.000401*cos(de2-md-ms) +
		e*.00032*cos(md-ms) -
		.000271*cos(de) -
		e*.000264*cos(ms+md) -
		.000198*cos(f2-md) +
		.000173*cos(md3) +
		.000167*cos(de4-md) -
		e*.000111*cos(ms) +
		.000103*cos(de4-md2) -
		.000084*cos(md2-de2) -
		e*.000083*cos(de2+ms) +
		.000079*cos(de2+md2) +
		.000072*cos(de4) +
		e*.000064*cos(de2-ms+md) -
		e*.000063*cos(de2+ms-md) +
		e*.000041*cos(ms+de) +
		e*.000035*cos(md2-ms) -
		.000033*cos(md3-de2) -
		.00003*cos(md+de) -
		.000029*cos(2*(f-de)) -
		e*.000029*cos(md2+ms) +
		e2*.000026*cos(2*(de-ms)) -
		.000023*cos(2*(f-de)+md) +
		e*.000019*cos(de4-ms-md)

	// distance from Earth in A.U.
	pos.Delta = 8.794 / (parallax * 3600)

	// angular speed
	motion = 13.176397 +
		1.434006*cos(md) +
		0.280135*cos(de2) +
		0.251632*cos(de2-md) +
		0.097420*cos(md2) -
		0.052799*cos(f2) +
		0.034848*cos(de2+md) +
		0.018732*cos(de2-ms) +
		0.010316*cos(de2-ms-md) +
		0.008649*cos(ms-md) -
		0.008642*cos(f2+md) -
		0.007471*cos(ms+md) -
		0.007387*cos(de) +
		0.006864*cos(md2+md) +
		0.006650*cos(de4-md) +
		0.003523*cos(de2+md2) +
		0.003377*cos(de4-md2) +
		0.003287*cos(de4) -
		0.003193*cos(ms) -
		0.003003*cos(de2+ms) +
		0.002577*cos(md-ms+de2) -
		0.002567*cos(f2-md) -
		0.001794*cos(de2-md2) -
		0.001716*cos(md-f2-de2) -
		0.001698*cos(de2+ms-md) -
		0.001415*cos(de2+f2) +
		0.001183*cos(md2-ms) +
		0.001150*cos(de+ms) -
		0.001035*cos(de+md) -
		0.001019*cos(f2+md2) -
		0.001006*cos(ms+md2)

	return
}
