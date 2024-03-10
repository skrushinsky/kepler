package core

func Map(data []float64, f func(float64) float64) []float64 {

	res := make([]float64, 0, len(data))

	for _, e := range data {
		res = append(res, f(e))
	}

	return res
}
