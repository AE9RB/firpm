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
  "fmt"
  "time"
  "os/exec"
  "bytes"
  "runtime"
)

func run_c_test() {
  var h [104]float64  
  weights := []float64{10,1,3,1,20}
  desired := []float64{0,1,0,1,0}
  bands := []float64{0,0.05,0.1,0.15,0.18,0.25,0.3,0.36,0.41,0.5}
  
  // Ignored test run
  C.remez((*C.double)(&h[0]), (C.int)(len(h)), (C.int)(len(weights)), (*C.double)(&bands[0]), (*C.double)(&desired[0]), (*C.double)(&weights[0]), C.BANDPASS)
  // Actual run
  t0 := time.Now()
  C.remez((*C.double)(&h[0]), (C.int)(len(h)), (C.int)(len(weights)), (*C.double)(&bands[0]), (*C.double)(&desired[0]), (*C.double)(&weights[0]), C.BANDPASS)
  t1 := time.Now()
  fmt.Printf("C: %v\n", t1.Sub(t0))
  
  // fmt.Println(h)
  
}

func run_ruby_test() {
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

func main() {
  
  run_c_test()
  run_ruby_test()  
  
}
