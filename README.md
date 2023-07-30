

# Scaliger

Library of core routines for practical astronomy. It is named in honor of 
*Johannes Kepler (1571-1630)*, German astronomer, mathematician, astrologer, natural philosopher and writer on music.

## Features

* Accurate positions of Sun, Moon and the planets, including _Pluto_.
* Time of solstices, equinoxes and lunations.

## Getting started

Add to pubspec.yaml of your project

```yaml
dependencies:
  ...
  scaliger:
    git:
      url: https://github.com/skrushinsky/kepler.git
      ref: master
  ...
```

Then run:

```console
$ dart pub update
```


## Usage

```dart
// Apparent position of the Sun
import 'package:kepler/sun.dart';

final djd = 33888.5; // 1992, Oct. 13 0h in days since epoch 1900,0
final (lambda, delta) = apparent(djd); // 199.9, .9975999344847888
```

```dart
// Apparent position of the Mars
import 'package:kepler/ephemeris.dart';

final eph = Ephemeris.forDJD(30830.5); // 1984 May 30
final pos = eph.geocentricPosition("Mercury");

print('lambda: ${pos.lambda}, beta: ${pos.beta}, delta: ${pos.delta}')
```


## Unit tests

```console
$ dart test ./test
```

## Additional information

TODO

## Sources

The formulae were adopted from the following sources:

* _Peter Duffett-Smith, "Astronomy With Your Personal Computer", Cambridge University Press, 1997_
* _Jean Meeus, "Astronomical Algorithms", 2d edition, Willmann-Bell, 1998_
* _J.L.Lawrence, "Celestial Calculations", The MIT Press, 2018_


## How to contribute

You may contribute to the project by many different ways, starting from refining and correcting its documentation,
especially if you are a native English speaker, and ending with improving the code base. Any kind of testing and
suggestions are welcome.

You may follow the standard Github procedures or, in case you are not comfortable with them, just send your suggestions
to the author by other means.