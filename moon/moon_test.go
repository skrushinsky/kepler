package moon

import (
	"testing"

	"github.com/skrushinsky/kepler/core"
	"github.com/skrushinsky/scaliger/julian"
	"github.com/skrushinsky/scaliger/mathutils"
)

const _DELTA = 1e-4

func TestMoonPosition(t *testing.T) {
	type _TestCase struct {
		djd      float64
		pos      core.EclipticPosition
		parallax float64
		motion   float64
	}
	var cases = [...]_TestCase{
		{
			djd:      -1.000050e+04,
			pos:      core.EclipticPosition{Lambda: 253.85478, Beta: -0.35884, Delta: 0.002475},
			parallax: 0.98681,
			motion:   14.073505,
		},
		{
			djd:      -7.000500e+03,
			pos:      core.EclipticPosition{Lambda: 183.03298, Beta: -5.10613, Delta: 0.0025318451878263725},
			parallax: 0.96482,
			motion:   13.614904285991807,
		},
		{
			djd:      -4.000500e+03,
			pos:      core.EclipticPosition{Lambda: 114.49714, Beta: 0.29899, Delta: 0.002661387458557927},
			parallax: 0.91786,
			motion:   12.284203442108854,
		},
		{
			djd:      -1.000500e+03,
			pos:      core.EclipticPosition{Lambda: 46.33258, Beta: 5.03904, Delta: 0.0027150753763781643},
			parallax: 0.89971,
			motion:   11.86016463804351,
		},
		{
			djd:      1.999500e+03,
			pos:      core.EclipticPosition{Lambda: 340.74811, Beta: -0.76686, Delta: 0.002665042330735118},
			parallax: 0.91660,
			motion:   12.137096046101872,
		},
		{
			djd:      4.999500e+03,
			pos:      core.EclipticPosition{Lambda: 273.11888, Beta: -5.22297, Delta: 0.0026145243283283597},
			parallax: 0.93431,
			motion:   12.706509343283184,
		},
		{
			djd:      7.999500e+03,
			pos:      core.EclipticPosition{Lambda: 198.76809, Beta: 0.13467, Delta: 0.0025060386483159646},
			parallax: 0.97476,
			motion:   13.79049510733593,
		},
		{
			djd:      1.099950e+04,
			pos:      core.EclipticPosition{Lambda: 123.17331, Beta: 5.01217, Delta: 0.002393311553917354},
			parallax: 1.02067,
			motion:   15.25014893599962,
		},
		{
			djd:      1.399950e+04,
			pos:      core.EclipticPosition{Lambda: 50.40519, Beta: 0.59539, Delta: 0.002440903741394359},
			parallax: 1.00077,
			motion:   14.567332957243783,
		},
		{
			djd:      1.699950e+04,
			pos:      core.EclipticPosition{Lambda: 336.88148, Beta: -5.04905, Delta: 0.0025896311772431097},
			parallax: 0.94329,
			motion:   13.015006558327384,
		},
		{
			djd:      1.999950e+04,
			pos:      core.EclipticPosition{Lambda: 266.43192, Beta: -1.18331, Delta: 0.0026726946153555506},
			parallax: 0.91398,
			motion:   12.05705112860313,
		},
		{
			djd:      2.299950e+04,
			pos:      core.EclipticPosition{Lambda: 200.91657, Beta: 5.13843, Delta: 0.00270357511434672},
			parallax: 0.90354,
			motion:   11.883519914105939,
		},
		{
			djd:      2.599950e+04,
			pos:      core.EclipticPosition{Lambda: 134.05765, Beta: 0.87204, Delta: 0.0026941433274419316},
			parallax: 0.90670,
			motion:   11.945823078908266,
		},
		{
			djd:      2.899950e+04,
			pos:      core.EclipticPosition{Lambda: 64.16216, Beta: -4.94147, Delta: 0.0025731373392409293},
			parallax: 0.94934,
			motion:   13.2314091077357,
		},
		{
			djd:      3.199950e+04,
			pos:      core.EclipticPosition{Lambda: 354.53313, Beta: -0.77311, Delta: 0.0024513561792419898},
			parallax: 0.99650,
			motion:   14.398538661582212,
		},
		{
			djd:      3.499950e+04,
			pos:      core.EclipticPosition{Lambda: 280.10165, Beta: 5.06817, Delta: 0.002455022531559789},
			parallax: 0.99501,
			motion:   14.431229034360273,
		},
		{
			djd:      3.799950e+04,
			pos:      core.EclipticPosition{Lambda: 201.62149, Beta: 2.25573, Delta: 0.0025070947174279036},
			parallax: 0.97435,
			motion:   13.731560363493482,
		},
		{
			djd:      4.099950e+04,
			pos:      core.EclipticPosition{Lambda: 128.41649, Beta: -4.51661, Delta: 0.0025554365866768364},
			parallax: 0.95591,
			motion:   13.279315343224834,
		},
		{
			djd:      4.399950e+04,
			pos:      core.EclipticPosition{Lambda: 61.54198, Beta: -2.45092, Delta: 0.0026505164857345337},
			parallax: 0.92162,
			motion:   12.374443595332336,
		},
		{
			djd:      4.699950e+04,
			pos:      core.EclipticPosition{Lambda: 353.93133, Beta: 4.49791, Delta: 0.00271630014162966},
			parallax: 0.89930,
			motion:   11.857387239871063,
		},
	}

	for _, test := range cases {
		pos, _, motion := TruePosition(test.djd + julian.J1900)
		if !mathutils.AlmostEqual(pos.Lambda, test.pos.Lambda, _DELTA) {
			t.Errorf("Expected Lambda: %f, got: %f", test.pos.Lambda, pos.Lambda)
		}
		if !mathutils.AlmostEqual(pos.Beta, test.pos.Beta, _DELTA) {
			t.Errorf("Expected Beta: %f, got: %f", test.pos.Beta, pos.Beta)
		}
		if !mathutils.AlmostEqual(pos.Delta, test.pos.Delta, _DELTA) {
			t.Errorf("Expected Delta: %f, got: %f", test.pos.Delta, pos.Delta)
		}
		if !mathutils.AlmostEqual(motion, test.motion, _DELTA) {
			t.Errorf("Expected Motion: %f, got: %f", test.motion, motion)
		}
	}
}

func TestMeanLunarNode(t *testing.T) {
	got := LunarNode(2438792.99027777778, true)
	exp := 80.31173473979322
	if !mathutils.AlmostEqual(got, exp, _DELTA) {
		t.Errorf("Expected Mean Lunar Node: %f, got: %f", exp, got)
	}
}

func TestTrueLunarNode(t *testing.T) {
	got := LunarNode(2438792.99027777778, false)
	exp := 81.86652882901491
	if !mathutils.AlmostEqual(got, exp, _DELTA) {
		t.Errorf("Expected True Lunar Node: %f, got: %f", exp, got)
	}
}
