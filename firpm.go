// Copyright (C) 2014 David Turnbull
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

package main

// #cgo CFLAGS: -O3
// #include "remez.c"
import "C"

import (
	"bytes"
	"errors"
	"fmt"
	"math"
	"os/exec"
	"runtime"
	"time"
)

func main() {
	runCTest()
	RunGoTest()
	runRubyTest()
}

func runCTest() {
	var h [104]float64
	weights := []float64{10, 1, 3, 1, 20}
	desired := []float64{0, 1, 0, 1, 0}
	bands := []float64{0, 0.05, 0.1, 0.15, 0.18, 0.25, 0.3, 0.36, 0.41, 0.5}
	t0 := time.Now()
	C.remez((*C.double)(&h[0]), (C.int)(len(h)), (C.int)(len(weights)), (*C.double)(&bands[0]), (*C.double)(&desired[0]), (*C.double)(&weights[0]), C.BANDPASS)
	t1 := time.Now()
	fmt.Printf("C: %v\n", t1.Sub(t0))
	// fmt.Println(h)
}

func RunGoTest() {
	h := make([]float64, 104)
	edge := []float64{0, 0.05, 0.1, 0.15, 0.18, 0.25, 0.3, 0.36, 0.41, 0.5}
	fx := []float64{0, 0, 1, 1, 0, 0, 1, 1, 0, 0}
	wtx := []float64{10, 1, 3, 1, 20}
	t0 := time.Now()
	if err := Firpm(h, edge, fx, wtx, BANDPASS, 0, 0); err != nil {
		fmt.Println(err)
	}
	t1 := time.Now()
	fmt.Printf("Go: %v\n", t1.Sub(t0))
	// fmt.Println(h)
}

func runRubyTest() {
	_, file, _, _ := runtime.Caller(0)
	file = file[:len(file)-2] + "rb"
	cmd := exec.Command("ruby", file)
	var stdout bytes.Buffer
	cmd.Stdout = &stdout
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	err := cmd.Run()
	if err != nil {
		fmt.Printf("Error running Ruby test:\n")
		fmt.Printf(stderr.String())
	} else {
		fmt.Printf(stdout.String())
	}
}

func Firpm(h, edge, fx, wtx []float64, jtype, lgrid, itrmax int) error {
	f := new(firpm)
	f.h, f.edge, f.fx, f.wtx = h, edge, fx, wtx
	f.jtype, f.lgrid, f.itrmax = jtype, lgrid, itrmax
	return f.run()
}

const (
	// jtypes
	BANDPASS       = 1
	DIFFERENTIATOR = 2
	HILBERT        = 3
	// math shortcuts
	Pi  = math.Pi
	Pi2 = math.Pi * 2
)

type firpm struct {
	// input parameters
	h, edge, fx, wtx     []float64
	jtype, lgrid, itrmax int
	// runtime variables
	neg, nodd     bool
	nbands, nfcns int
	iext, fext    []int
	grid, wt, des []float64
	e, x, y, ad   []float64
}

func (f *firpm) run() error {

	f.nbands = len(f.wtx)
	if f.nbands < 1 || len(f.fx) != f.nbands*2 || len(f.edge) != f.nbands*2 {
		return errors.New("input slice lengths do not make sense")
	}
	if f.jtype < 1 || f.jtype > 3 {
		return errors.New("unknown jtype")
	}
	if f.lgrid <= 0 {
		f.lgrid = 16
	}
	if f.itrmax <= 0 {
		f.itrmax = 40
	}
	f.neg = true
	if f.jtype == BANDPASS {
		f.neg = false
	}
	nfilt := len(f.h)
	f.nodd = nfilt%2 == 1
	f.nfcns = nfilt / 2
	if f.nodd && !f.neg {
		f.nfcns++
	}
	gridlen := int(0)
	for i := 0; i < f.nbands; i++ {
		gridlen += int(2.0*float64(f.nfcns*f.lgrid)*(f.edge[2*i+1]-f.edge[2*i]) + 0.5)
	}
	f.grid = make([]float64, gridlen)
	f.wt = make([]float64, gridlen)
	f.des = make([]float64, gridlen)
	f.e = make([]float64, gridlen)
	f.x = make([]float64, f.nfcns+1)
	f.y = make([]float64, f.nfcns+1)
	f.ad = make([]float64, f.nfcns+1)
	f.iext = make([]int, f.nfcns+1)
	f.fext = make([]int, f.nfcns*2)
	taps := make([]float64, f.nfcns+1)

	delf := 0.5 / float64(f.lgrid*f.nfcns)
	lowf := delf
	if !f.neg || f.edge[0] > delf {
		lowf = f.edge[0]
	}

	for j, band := 0, 0; band < f.nbands; band++ {
		f.grid[j] = f.edge[2*band]
		if band > 0 {
			lowf = f.edge[2*band]
		}
		highf := f.edge[2*band+1]
		k := int((highf-lowf)/delf + 0.5)
		for i := 0; i < k; i++ {
			f.des[j] = f.fx[2*band] + float64(i)*(f.fx[2*band+1]-f.fx[2*band])/float64(k-1)
			f.wt[j] = f.wtx[band]
			f.grid[j] = lowf
			lowf += delf
			j++
		}
		f.grid[j-1] = highf
	}
	if f.neg && f.grid[gridlen-1] > (0.5-delf) && f.nodd {
		f.grid[gridlen-1] = 0.5 - delf
	}

	if f.jtype == DIFFERENTIATOR {
		for i := 0; i < gridlen; i++ {
			if f.des[i] > 0.0001 {
				f.wt[i] = f.wt[i] / f.grid[i]
			}
		}
	}

	if f.neg {
		if f.nodd {
			for i := 0; i < gridlen; i++ {
				c := math.Sin(Pi2 * f.grid[i])
				f.des[i] /= c
				f.wt[i] *= c
			}
		} else {
			for i := 0; i < gridlen; i++ {
				c := math.Sin(Pi * f.grid[i])
				f.des[i] /= c
				f.wt[i] *= c
			}
		}
	} else {
		if !f.nodd {
			for i := 0; i < gridlen; i++ {
				c := math.Cos(Pi * f.grid[i])
				f.des[i] /= c
				f.wt[i] *= c
			}
		}
	}

	for i := 0; i <= f.nfcns; i++ {
		f.iext[i] = i * (gridlen - 1) / f.nfcns
	}

	itr := 0
	for ; itr < f.itrmax; itr++ {
		f.calcParams()
		f.calcError()
		f.search()
		if f.done() {
			break
		}
	}
	f.calcParams()

	for i, c := 0, 0.0; i <= nfilt/2; i++ {
		if f.neg {
			if f.nodd {
				c = math.Sin(Pi2 * float64(i) / float64(nfilt))
			} else {
				c = math.Sin(Pi * float64(i) / float64(nfilt))
			}
		} else {
			if f.nodd {
				c = 1
			} else {
				c = math.Cos(Pi * float64(i) / float64(nfilt))
			}
		}
		taps[i] = f.computeA(float64(i)/float64(nfilt)) * c
	}

	m := (float64(nfilt) - 1) / 2
	if f.neg {
		if f.nodd {
			for n := 0; n < nfilt; n++ {
				val := taps[0]
				x := Pi2 * (float64(n) - m) / float64(nfilt)
				for k := 1; k <= int(m); k++ {
					val += 2 * taps[k] * math.Sin(x*float64(k))
				}
				f.h[n] = val / float64(nfilt)
			}
		} else {
			for n := 0; n < nfilt; n++ {
				val := taps[0]
				x := Pi2 * (float64(n) - m) / float64(nfilt)
				for k := 1; k <= (nfilt/2 - 1); k++ {
					val += 2 * taps[k] * math.Sin(x*float64(k))
				}
				f.h[n] = val / float64(nfilt)
			}
		}
	} else {
		if f.nodd {
			for n := 0; n < nfilt; n++ {
				val := taps[0]
				x := Pi2 * (float64(n) - m) / float64(nfilt)
				for k := 1; k <= int(m); k++ {
					val += 2 * taps[k] * math.Cos(x*float64(k))
				}
				f.h[n] = val / float64(nfilt)
			}
		} else {
			for n := 0; n < nfilt; n++ {
				val := taps[0]
				x := Pi2 * (float64(n) - m) / float64(nfilt)
				for k := 1; k <= (nfilt/2 - 1); k++ {
					val += 2 * taps[k] * math.Cos(x*float64(k))
				}
				f.h[n] = val / float64(nfilt)
			}
		}
	}

	if itr == f.itrmax {
		return errors.New("Maximum iterations reached. Results may not be valid.")
	}
	return nil
}

func (f *firpm) calcParams() {
	for i := 0; i <= f.nfcns; i++ {
		f.x[i] = math.Cos(Pi2 * f.grid[f.iext[i]])
	}
	ld := (f.nfcns-1)/15 + 1
	for i := 0; i <= f.nfcns; i++ {
		denom := 1.0
		xi := f.x[i]
		for j := 0; j < ld; j++ {
			for k := j; k <= f.nfcns; k += ld {
				if k != i {
					denom *= 2 * (xi - f.x[k])
				}
			}
		}
		if math.Abs(denom) < 0.00001 {
			denom = 0.00001
		}
		f.ad[i] = 1 / denom
	}
	numer, denom := 0.0, 0.0
	sign := 1.0
	for i := 0; i <= f.nfcns; i++ {
		numer += f.ad[i] * f.des[f.iext[i]]
		denom += sign * f.ad[i] / f.wt[f.iext[i]]
		sign = -sign
	}
	delta := numer / denom
	sign = 1
	for i := 0; i <= f.nfcns; i++ {
		f.y[i] = f.des[f.iext[i]] - sign*delta/f.wt[f.iext[i]]
		sign = -sign
	}
}

func (f *firpm) calcError() {
	for i := range f.e {
		f.e[i] = f.wt[i] * (f.des[i] - f.computeA(f.grid[i]))
	}
}

func (f *firpm) computeA(freq float64) float64 {
	numer, denom := 0.0, 0.0
	xc := math.Cos(Pi2 * freq)
	for i := 0; i <= f.nfcns; i++ {
		c := xc - f.x[i]
		if math.Abs(c) < 1.0e-7 {
			numer = f.y[i]
			denom = 1
			break
		}
		c = f.ad[i] / c
		denom += c
		numer += c * f.y[i]
	}
	return numer / denom
}

func (f *firpm) search() {
	k := 0
	if ((f.e[0] > 0) && (f.e[0] > f.e[1])) ||
		((f.e[0] < 0) && (f.e[0] < f.e[1])) {
		f.fext[k] = 0
		k++
	}
	j := len(f.e) - 1
	for i := 1; i < j; i++ {
		if ((f.e[i] >= f.e[i-1]) && (f.e[i] > f.e[i+1]) && (f.e[i] > 0)) ||
			((f.e[i] <= f.e[i-1]) && (f.e[i] < f.e[i+1]) && (f.e[i] < 0)) {
			f.fext[k] = i
			k++
		}
	}
	if ((f.e[j] > 0.0) && (f.e[j] > f.e[j-1])) ||
		((f.e[j] < 0.0) && (f.e[j] < f.e[j-1])) {
		f.fext[k] = j
		k++
	}
	for extra := k - (f.nfcns + 1); extra > 0; extra-- {
		up := f.e[f.fext[0]] > 0
		l := 0
		alt := true
		for j := 1; j < k; j++ {
			if math.Abs(f.e[f.fext[j]]) < math.Abs(f.e[f.fext[l]]) {
				l = j
			}
			if up && f.e[f.fext[j]] < 0 {
				up = false
			} else if !up && f.e[f.fext[j]] > 0 {
				up = true
			} else {
				alt = false
				break
			}
		}
		if alt && extra == 1 {
			if math.Abs(f.e[f.fext[k-1]]) < math.Abs(f.e[f.fext[0]]) {
				l = f.fext[k-1]
			} else {
				l = f.fext[0]
			}
		}
		for j := l; j < k; j++ {
			f.fext[j] = f.fext[j+1]
		}
		k--
	}
	for i := 0; i <= f.nfcns; i++ {
		f.iext[i] = f.fext[i]
	}
}

func (f *firpm) done() bool {
	min := math.Abs(f.e[f.iext[0]])
	max := min
	for i := 1; i <= f.nfcns; i++ {
		current := math.Abs(f.e[f.iext[i]])
		if current < min {
			min = current
		}
		if current > max {
			max = current
		}
	}
	return ((max - min) / max) < 0.0001
}
