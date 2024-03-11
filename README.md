# Kepler

Library of core routines for practical astronomy. It is named in honor of 
*Johannes Kepler (1571-1630)*, German astronomer, mathematician, astrologer and natural philosopher.

- [Kepler](#kepler)
  - [Features](#features)
  - [Quick Start](#quick-start)
  - [Usage](#usage)
    - [Sun and Moon](#sun-and-moon)
    - [Planets](#planets)
    - [Utilities](#utilities)
  - [See also](#see-also)
  - [How to contribute](#how-to-contribute)
  - [Sources](#sources)


## Features

* Accurate positions of Sun, Moon and the planets, including _Pluto_.
* Time of solstices, equinoxes and lunations.

## Quick Start

```console

$ go get github.com/skrushinsky/kepler

```

## Usage

### Sun and Moon

* `sun.TrueGeocentric(t, ms, ls float64) (lsn float64, rsn float64)` calculates true geocentric longitude of the Sun for the mean equinox of date and the Sun-Earth distance.
* `sun.Apparent(jd float64, options ApparentSunOptions) core.EclipticPosition` apparent geocentric ecliptical longitude of the Sun.
* `sun.MeanLongitude(t float64) float64` Mean longitude of the Sun.
* `sun.MeanAnomaly(t float64) float64` Mean anomaly of the Sun. 
* `moon.MeanLunarNode(t float64) float64` Mean Lunar Node. 
* `moon.LunarNode(jd float64, mean bool) float64` Mean or True Lunar Node.
* `moon.TruePosition(jd float64) (pos core.EclipticPosition, parallax, motion float64)` Moon position, horizontal parallax and daily motion for mean equinox of date.

### Planets

TODO

### Utilities

* `core.EccentricAnomaly(s, m, ea float64) float64` solves Kepler equation.
* `core.TrueAnomaly(s, ea float64) float64` Given **s**, eccentricity, and **ea**, eccentric anomaly, finds true anomaly.
* `core.Map(data []float64, f func(float64) float64) []float64` applies **f** function to each element of **data** slice.

## See also

[Library of date/time manipulation routines for practical astronomy](https://github.com/skrushinsky/scaliger)


## How to contribute

You may contribute to the project by many different ways, starting from refining and correcting its documentation,
especially if you are a native English speaker, and ending with improving the code base. Any kind of testing and suggestions is welcome.

You may follow the standard Github procedures or, in case you are not comfortable with them, just send your suggestions
to the author by other means.

## Sources

The formulae were adopted from the following sources:

* _Peter Duffett-Smith, "Astronomy With Your Personal Computer", Cambridge University Press, 1997_
* _Jean Meeus, "Astronomical Algorithms", 2d edition, Willmann-Bell, 1998_
* _J.L.Lawrence, "Celestial Calculations", The MIT Press, 2018_


