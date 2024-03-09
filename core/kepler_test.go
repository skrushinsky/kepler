package core

import (
	"testing"

	"github.com/skrushinsky/scaliger/mathutils"
)

const _DELTA = 1e-4

type _TestCase struct {
	m  float64
	s  float64
	ea float64
	ta float64
}

var cases = [...]_TestCase{
	{
		m:  3.5208387374141448,
		s:  0.016718,
		ea: 3.5147440476661806,
		ta: -2.774497552017826,
	},
	{
		m:  0.763009079752865,
		s:  0.965,
		ea: 1.7176273861066755,
		ta: 2.9122563898777387,
	},
}

func TestEccentricAnomaly(t *testing.T) {
	for _, test := range cases {
		ea := EccentricAnomaly(test.s, test.m, test.m)
		if !mathutils.AlmostEqual(ea, test.ea, _DELTA) {
			t.Errorf("Expected: %f, got: %f", test.ea, ea)
		}
	}
}

func TestTrueAnomaly(t *testing.T) {
	for _, test := range cases {
		ta := TrueAnomaly(test.s, test.ea)
		if !mathutils.AlmostEqual(ta, test.ta, _DELTA) {
			t.Errorf("Expected: %f, got: %f", test.ta, ta)
		}
	}
}
