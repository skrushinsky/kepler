package sun

import (
	"testing"

	"github.com/skrushinsky/scaliger/julian"
	"github.com/skrushinsky/scaliger/mathutils"
	"github.com/skrushinsky/scaliger/nutequ"
)

const _DELTA = 1e-4

type _TestCase struct {
	djd float64
	l   float64
	r   float64
	ap  float64
}

var cases = [...]_TestCase{
	{
		djd: 30916.5, // 24 Aug 1984 00:00
		l:   151.01309547440778,
		r:   1.010993800005251,
		ap:  151.0035132296576,
	},
	{
		djd: 30819.10833333333, // 18 May 1984 14:36
		l:   57.83143688493146,
		r:   1.011718488789592,
		ap:  57.82109236581925,
	},
	{
		djd: 28804.5, // 12 Nov 1978 00:00
		l:   229.2517039627867,
		r:   0.9898375,
		ap:  229.2450957063683,
	},
	{
		djd: 33888.5, // 1992, Oct. 13 0h
		l:   199.90600618015975,
		r:   .9975999344847888,
		ap:  199.9047664927989, // Meeus: 199.90734722222223
	},
}

func TestTrueGeocentric(t *testing.T) {
	for _, test := range cases {
		tperiod := test.djd / julian.DAYS_PER_CENT
		lsn, rsn := TrueGeocentric(tperiod, MeanAnomaly(tperiod), MeanLongitude(tperiod))
		if !mathutils.AlmostEqual(lsn, test.l, _DELTA) {
			t.Errorf("Expected: %f, got: %f", test.l, lsn)
		}
		if !mathutils.AlmostEqual(rsn, test.r, _DELTA) {
			t.Errorf("Expected: %f, got: %f", test.r, rsn)
		}
	}
}

func TestApparent(t *testing.T) {
	for _, test := range cases {
		tperiod := test.djd / julian.DAYS_PER_CENT
		jd := test.djd + julian.J1900
		dpsi, _ := nutequ.Nutation(jd)
		opts := ApparentSunOptions{
			meanAnomaly:       MeanAnomaly(tperiod),
			meanLongitude:     MeanLongitude(tperiod),
			ignoreLightTravel: true,
			dpsi:              dpsi,
		}
		pos := Apparent(jd, opts)
		if !mathutils.AlmostEqual(pos.Lambda, test.ap, _DELTA) {
			t.Errorf("Expected: %f, got: %f", test.ap, pos.Lambda)
		}
		if !mathutils.AlmostEqual(pos.Delta, test.r, _DELTA) {
			t.Errorf("Expected: %f, got: %f", test.r, pos.Delta)
		}
		if !mathutils.AlmostEqual(pos.Beta, 0, _DELTA) {
			t.Errorf("Expected: 0, got: %f", pos.Beta)
		}
	}
}
