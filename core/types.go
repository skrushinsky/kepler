package core

// Position of a celestial body on the Ecliptic plane
type EclipticPosition struct {
	// celestial longitude, degrees
	Lambda float64
	// celestial latitude, degrees
	Beta float64
	// distance from Earth, A.U.
	Delta float64
}
