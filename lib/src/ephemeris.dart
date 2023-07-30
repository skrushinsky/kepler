import 'dart:math';
import 'package:vector_math/vector_math_64.dart';

import 'package:scaliger/mathutils.dart';
import 'sun.dart' as sun;
import 'moon.dart' as moon;
import 'ephemeris/planets.dart';
import 'package:scaliger/nutequ.dart';
import 'ephemeris/pert.dart';

/// Ecliptic posiion of a celestial body
///
/// * [lambda] : geocentric ecliptic longitude, arc-degrees
/// * [beta]: geocentric ecliptic latitude, arc-degrees
/// * [delta] : distance from the Earth, A
class EclipticPosition {
  final double lambda;
  final double beta;
  final double delta;

  /// Constructor
  const EclipticPosition(this.lambda, this.beta, this.delta);
}

class Ephemeris {
  final double _djd; // days since Jan 1900, 0
  final double _t; // time in centuries since Jan 1900, 0
  final double _ms; // Sun mean anomaly in arc-degrees

  final ({
    double dpsi,
    double deps
  }) _nut; // nutation in longitude and obliquity
  // double? _dpsi; // nutation in longitude, arc-degrees
  // double? _deps; // nutation in eclipic obliquity, arc-degrees
  final double _eps; // obliquity of the elipic, arc-degrees

  // lsn: longitude of Sun, arc-degrees
  // rsn: Sun-Earth distance, AU
  final ({double lsn, double rsn}) _sunGeo;

  final bool _apparent;
  final bool _trueNode;

  final Map<String, EclipticPosition> _positions = {};
  final Map<String, double> _dailyMotions = {};
  final Map<String, double> _meanAnomalies = {};
  final Map<String, double> _eccentricities = {};

  Ephemeris? _prev;
  Ephemeris? _next;

  // _positions = Map();

  Ephemeris(this._djd, this._apparent, this._trueNode, this._nut, this._eps,
      this._sunGeo, this._t, this._ms);

  factory Ephemeris.forDJD(
    double djd, {
    apparent = false,
    trueNode = true,
  }) {
    final t = djd / 36525;
    final nut = nutation(t);
    final eps = obliquity(djd, deps: nut.deps);
    final sunGeo = sun.trueGeocentric(t, ms: sun.meanAnomaly(t));
    final ms = sun.meanAnomaly(t);
    return Ephemeris(djd, apparent, trueNode, nut, eps, sunGeo, t, ms);
  }

  /// Days from Jan 0.5 1900
  double get djd => _djd;

  /// Centuries from Jan 0.5 1900
  double get t => _t;

  ///  Apparent flag
  bool get apparent => _apparent;

  /// True Node flag
  bool get trueNode => _trueNode;

  /// Sun mean anomaly
  double get ms => _ms;

  ///  /// Nutation in longitude and obliquity, arc-degrees
  ({double dpsi, double deps}) get nut => _nut;

  /// Obliquity of the ecliptic, arc-degrees
  double? get eps => _eps;

  /// True geocentric position of the Sun:
  /// longitude, arc-degrees and  Sun-Earth distance, AU
  ({double lsn, double rsn}) get sunGeo => _sunGeo;

  /// Ephemeris instance 12h before
  Ephemeris get prev {
    _prev ??=
        Ephemeris.forDJD(_djd - 0.5, apparent: _apparent, trueNode: _trueNode);
    return _prev!;
  }

  /// Ephemeris instance 12h ahead
  Ephemeris get next {
    _next ??=
        Ephemeris.forDJD(_djd + 0.5, apparent: _apparent, trueNode: _trueNode);
    return _next!;
  }

  EclipticPosition _calculateSun() {
    double x = _sunGeo.lsn;
    if (apparent) {
      // nutation and aberration
      x = reduceDeg(x + _nut.dpsi - 5.69e-3);
      // XXX: Peter Duffett-Smith does not mention light-time travel correction
      // in case of the Sun.
      // // light travel
      // final lt = 1.365 * rsn; // seconds
      // x -= lt * 15 / 3600;
    }
    return EclipticPosition(x, 0.0, sunGeo.rsn);
  }

  EclipticPosition _calculateMoon() {
    EclipticPosition pos;
    final m = moon.truePosition(djd);
    var lambda = m.lambda;
    if (apparent) {
      lambda = reduceDeg(lambda + _nut.dpsi);
    }
    pos = EclipticPosition(lambda, m.beta, m.delta);
    _dailyMotions["Moon"] = m.motion;

    return pos;
  }

  /// Mean Anomaly of a planet [pla] in radians.
  /// Once calculated, it is saved in a cache.
  /// [dt] named parameter is a time corection necessary when calculating *true*
  /// (light-time corrected) planetary positions.
  double meanAnomaly(Planet pla, {dt = 0}) {
    double ma;

    if (_meanAnomalies.containsKey(pla.name)) {
      ma = _meanAnomalies[pla.name]!;
    } else {
      ma = radians(pla.orbit.meanAnomaly(t));
      _meanAnomalies[pla.name] = ma;
    }
    return ma - radians(dt * pla.orbit.DM);
  }

  /// Eccentricity of a planet [pla].
  /// Once calculated, it is saved in a cache.
  double eccentricity(Planet pla) {
    double ec;
    if (_eccentricities.containsKey(pla.name)) {
      ec = _meanAnomalies[pla.name]!;
    } else {
      ec = pla.orbit.EC.assemble(t);
      _eccentricities[pla.name] = ec;
    }
    return ec;
  }

  Map<String, double> _planetHelio(pla,
      {double? lg,
      double? s,
      double? sa,
      double? ph,
      double? nd,
      double? ic,
      double dt = 0.0}) {
    final m = meanAnomaly(pla, dt: dt);
    Map<PertType, double> Function() calcPert;
    switch (pla.id) {
      case PlanetId.Mercury:
        calcPert = () {
          final ve = meanAnomaly(Planet.forId(PlanetId.Venus), dt: dt);
          final ju = meanAnomaly(Planet.forId(PlanetId.Jupiter), dt: dt);
          return pla.calculatePerturbations([m, ve, ju]);
        };
        break;
      case PlanetId.Venus:
        calcPert = () {
          final ju = meanAnomaly(Planet.forId(PlanetId.Jupiter), dt: dt);
          final sm = radians(ms); // Sun, radians
          return pla.calculatePerturbations([t, sm, m, ju]);
        };
        break;
      case PlanetId.Mars:
        calcPert = () {
          final ve = meanAnomaly(Planet.forId(PlanetId.Venus), dt: dt);
          final ju = meanAnomaly(Planet.forId(PlanetId.Jupiter), dt: dt);
          final sm = radians(ms); // Sun, radians
          return pla.calculatePerturbations([sm, ve, m, ju]);
        };
        break;
      case PlanetId.Pluto:
        calcPert = () {
          return pla.calculatePerturbations(<double>[]);
        };
        break;
      default:
        calcPert = () {
          var args = List<double>.from([t, s]);
          return pla.calculatePerturbations(args);
        };
        break;
    }

    return pla.heliocentric(t,
        s: s,
        sa: sa,
        ph: ph,
        nd: nd,
        ic: ic,
        ma: m,
        re: sunGeo.rsn,
        lg: lg,
        pertFunc: calcPert);
  }

  EclipticPosition _calculatePlanet(Planet pla) {
    final lg = radians(sunGeo.lsn) + pi;
    final s = eccentricity(pla);

    final sa = pla.orbit.SA;
    final ph = radians(pla.orbit.PH.assemble(t));
    final nd = radians(pla.orbit.ND.assemble(t));
    final ic = radians(pla.orbit.IN.assemble(t));

    final h1 = _planetHelio(pla, lg: lg, s: s, sa: sa, ph: ph, nd: nd, ic: ic);
    final h2 = _planetHelio(pla,
        lg: lg,
        s: s,
        sa: sa,
        ph: ph,
        nd: nd,
        ic: ic,
        dt: h1['rho']! * 5.775518e-3);

// ll, rpd, lpd, sin(psi), cpsi, rho

    final geo = pla.geocentric(
        lg: lg,
        rsn: _sunGeo.rsn,
        lpd: h2['lpd']!,
        rpd: h2['rpd']!,
        cpsi: h2['cpsi']!,
        spsi: h2['spsi']!,
        ll: h2['ll']!,
        apparent: apparent,
        dpsi: _nut.dpsi);

    return EclipticPosition(geo[0], geo[1], h1['rho']!);
  }

  /// Given an object name, return its geocentric ecliptic position.
  ///
  /// Latitide and longitude are corrected for light time.
  /// These are the apparent values as seen from the center of the Earth
  /// at the given instant. If Ephemeris instance 'apparent' flag is set
  /// to 'true' (default), then corrections for nutation and aberration
  /// are also applied.
  EclipticPosition geocentricPosition(String name) {
    if (_positions.containsKey(name)) {
      return _positions[name]!;
    }

    EclipticPosition pos;
    switch (name) {
      case "Sun":
        pos = _calculateSun();
        break;
      case "Moon":
        pos = _calculateMoon();
        break;
      case "LunarNode":
        pos =
            EclipticPosition(moon.lunarNode(djd, trueNode: trueNode), 0.0, 0.0);
        break;
      default:
        Planet pla = Planet.forName(name);
        pos = _calculatePlanet(pla);
    }

    _positions[name] = pos;
    return pos;
  }

  /// Daily motion of a celestial body, arc-deg
  double dailyMotion(String name) {
    if (_dailyMotions.containsKey(name)) {
      return _dailyMotions[name]!;
    }

    if (name == "Moon") {
      // as a side effect _calculateMoon() will initialize _dailyMotions cache
      _calculateMoon();
      return dailyMotion(name);
    }

    final x0 = prev.geocentricPosition(name).lambda;
    final x1 = next.geocentricPosition(name).lambda;
    final delta = x1 - x0;
    _dailyMotions[name] = delta;
    return delta;
  }
}
