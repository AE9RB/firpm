#!/usr/bin/env swift -Ounchecked

//  Parks-McClellan algorithm for FIR filter design
//
//  Copyright (c) 2016 David Turnbull
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Library General Public
//  License as published by the Free Software Foundation; either
//  version 2 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Library General Public License for more details.
//
//  You should have received a copy of the GNU Library General Public
//  License along with this library; if not, write to the Free
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// Swift version of "remez.c" by Jake Janovetz (janovetz@uiuc.edu).
// Not an exact port; contains various fixes and enhancements.


import Foundation
#if os(Linux)
    import Glibc
#else
    import Darwin.C
#endif


public enum Filter {
    case BANDPASS
    case DIFFERENTIATOR
    case HILBERT
}


public func firpm(numTaps:Int,
    bands:[Double],
    desired:[Double],
    weights:[Double],
    filter:Filter
    ) -> [Double]
{

    let MAXITERATIONS = 40
    let GRIDDENSITY = 16

    let taps = computeTaps(numTaps: numTaps, filter: filter)

    let grid = computeGrid(
        taps: taps,
        density: GRIDDENSITY,
        bands: bands,
        desired: desired,
        weights: weights
    )

    var E = [Double](count: grid.count, repeatedValue: 0)
    var Ext = [Int](count: taps.extrema+1, repeatedValue: 0)
    for i in 0..<Ext.count {
        Ext[i] = i * (E.count - 1) / (Ext.count - 1)
    }

    var A = AmplitudeResponse(taps: taps, grid: grid)
    A.calcParms(Ext)

    var iter = 0
    while iter < MAXITERATIONS {
        iter += 1
        A.calcParms(Ext)

        /// Calculate Error function
        for i in 0..<grid.count
        {
            let a = A[grid[i].G]
            E[i] = grid[i].W * (grid[i].D - a)
        }

        if search(E, Ext: &Ext) {
            // print("Iterations: \(iter)")
            A.calcParms(Ext)
            return A.freqSample(taps)
        }
    }
    fatalError("Exceeds maximum iteration count.")
}


private enum Symmetry {
    case NEGATIVE
    case POSITIVE
}


private typealias TapsType = (
    count: Int,
    extrema: Int,
    isOdd: Bool,
    symmetry: Symmetry,
    filter: Filter
)


private func computeTaps(
    numTaps numTaps: Int,
    filter: Filter) -> TapsType
{
    let isOdd = numTaps % 2 == 1
    var symmetry:Symmetry

    var extrema = numTaps/2
    if (filter == .BANDPASS) {
        if isOdd {
            extrema += 1
        }
        symmetry = .POSITIVE
    } else {
        symmetry = .NEGATIVE
    }

    return (
        count:numTaps,
        isOdd:isOdd,
        extrema:extrema,
        symmetry:symmetry,
        filter:filter
    )
}


private typealias GridType = (G:Double, D:Double, W:Double)


private func computeGrid(
    taps taps: TapsType,
    density: Int,
    bands: [Double],
    desired: [Double],
    weights: [Double]) -> [GridType]
{
    // Predict dense grid size in advance for memory allocation
    //  .5 is so we round up, not truncate
    var gridsize = 0
    for i in 0..<bands.count/2
    {
        var gridinc = Double(2 * (taps.extrema) * density)
        gridinc *= bands[2*i+1] - bands[2*i]
        gridinc += 0.5
        gridsize += Int(gridinc)
    }
    if (taps.symmetry == .NEGATIVE)
    {
        gridsize -= 1
    }

    var grid = [GridType](count: gridsize, repeatedValue: (0,0,0))

    // Create the dense grid of frequencies from the specified bands.
    // Also creates the Desired Frequency Response function (D[]) and
    // the Weight function (W[]) on that dense grid
    var lowf: Double
    let delf = 0.5 / Double(density*(taps.extrema))

    if (taps.symmetry == .NEGATIVE) && (delf > bands[0]) {
        lowf = delf
    } else {
        lowf = bands[0]
    }

    var j = 0
    for band in 0..<bands.count/2 {
        grid[j].G = bands[2*band]
        if band != 0 {
            lowf = bands[2*band]
        }
        let highf = bands[2*band + 1]
        let k = Int((highf - lowf)/delf + 0.5)
        for i in 0..<k {
            var dj = Double(i)
            dj *= desired[2*band+1] - desired[2*band]
            dj /= Double(k-1)
            dj += desired[2*band]
            grid[j].D = dj
            grid[j].W = weights[band]
            grid[j].G = lowf
            lowf += delf
            j += 1
        }
        grid[j-1].G = highf
    }

    // Similar to above, if odd symmetry, last grid point can't be .5
    //  - but, if there are even taps, leave the last grid point at .5
    if taps.isOdd && taps.symmetry == .NEGATIVE &&
        grid[gridsize-1].G > (0.5 - delf)
    {
        grid[gridsize-1].G = 0.5-delf
    }

    // For Differentiator: (fix grid)
    if taps.filter == .DIFFERENTIATOR {
        for i in 0..<gridsize {
            if (grid[i].D > 0.0001) {
                grid[i].W = grid[i].W / grid[i].G
            }
        }
    }

    // For odd or Negative symmetry filters, alter the
    // D[] and W[] according to Parks McClellan
    if taps.symmetry == .POSITIVE {
        if !taps.isOdd {
            for i in 0..<gridsize {
                let c = cos(M_PI * grid[i].G)
                grid[i].D /= c
                grid[i].W *= c
            }
        }
    } else {
        if taps.isOdd {
            for i in 0..<gridsize {
                let c = sin(2 * M_PI * grid[i].G)
                grid[i].D /= c
                grid[i].W *= c
            }
        } else {
            for i in 0..<gridsize {
                let c = sin(M_PI * grid[i].G)
                grid[i].D /= c
                grid[i].W *= c
            }
        }
    }

    return grid
}


// AmplitudeResponse is A[]
private struct AmplitudeResponse {

    let grid:[GridType]
    var ad:[Double] // 'b' in Oppenheim & Schafer
    var x:[Double]
    var y:[Double]  // 'C' in Oppenheim & Schafer


    init(taps:TapsType, grid:[GridType]) {
        self.grid = grid
        ad = [Double](count: taps.extrema+1, repeatedValue: 0)
        x = [Double](count: taps.extrema+1, repeatedValue: 0)
        y = [Double](count: taps.extrema+1, repeatedValue: 0)
    }


    mutating func calcParms(Ext:[Int]) {
        // Find x[]
        for i in 0..<x.count {
            x[i] = cos(2 * M_PI * grid[Ext[i]].G )
        }

        // Calculate ad[] - Oppenheim & Schafer eq 7.132
        // Skips around to avoid round errors
        let ld = (Ext.count-2) / 15 + 1
        for i in 0..<ad.count {
            var denom = 1.0
            let xi = x[i]
            for j in 0..<ld {
                for k in j.stride(to: Ext.count, by: ld) {
                    if (k != i) {
                        denom *= 2 * (xi - x[k])
                    }
                }
            }
            if abs(denom) < 0.00001 {
                denom = 0.00001
            }
            ad[i] = 1 / denom
        }

        // Calculate delta  - Oppenheim & Schafer eq 7.131
        var numer = 0.0
        var denom = 0.0
        var sign = 1.0
        for i in 0..<ad.count {
            numer += ad[i] * grid[Ext[i]].D
            denom += sign * ad[i]/grid[Ext[i]].W
            sign = -sign
        }
        let delta = numer/denom
        sign = 1

        // Calculate y[]  - Oppenheim & Schafer eq 7.133b
        for i in 0..<y.count {
            y[i] = grid[Ext[i]].D - sign * delta/grid[Ext[i]].W
            sign = -sign
        }
    }

    /// Using values calculated in CalcParms, ComputeA calculates the
    /// actual filter response at a given frequency.  Uses
    /// eq 7.133a from Oppenheim & Schafer.
    subscript(freq:Double) -> Double
    {
        var denom = 0.0
        var numer = 0.0
        let xc = cos(M_PI*2 * freq)
        for i in 0..<ad.count
        {
            var c = xc - x[i]
            if (abs(c) < 1.0e-7)
            {
                numer = y[i]
                denom = 1
                break
            }
            c = ad[i]/c
            denom += c
            numer += c*y[i]
        }
        return numer/denom
    }


    /// Construct the filter from 'A[]'
    func freqSample(taps:TapsType) -> [Double]
    {
        if taps.symmetry == .POSITIVE {
            return freqSamplePositive(taps)
        }
        return freqSampleNegative(taps)
    }


    func freqSamplePositive(taps:TapsType) -> [Double]
    {
        let N = Double(taps.count)
        var t = [Double](count: x.count, repeatedValue: 0)
        for i in 0...taps.count/2 {
            let ii = Double(i)
            var c = self[ii / N]
            if !taps.isOdd {
                c *= cos(M_PI * ii / N)
            }
            t[i] = c
        }

        let M = Double(taps.count-1) / 2.0
        var h = [Double](count:taps.count, repeatedValue: 0)
        if taps.isOdd {
            for n in 0..<taps.count {
                var val = t[0]
                let x = 2 * M_PI * (Double(n) - M) / N
                for k in 1...Int(M) {
                    val += 2.0 * t[k] * cos(x*Double(k))
                }
                h[n] = val / N
            }
        } else {
            for n in 0..<taps.count {
                var val = t[0]
                let x = M_PI*2 * (Double(n) - M) / N
                for k in 1..<taps.count/2 {
                    val += 2.0 * t[k] * cos(x*Double(k))
                }
                h[n] = val / N
            }
        }
        return h
    }


    func freqSampleNegative(taps:TapsType) -> [Double]
    {
        let N = Double(taps.count)
        var t = [Double](count: x.count, repeatedValue: 0)
        for i in 0...taps.count/2 {
            let ii = Double(i)
            var c = self[ii / N]
            if taps.isOdd {
                c *= sin(2 * M_PI * ii / N)
            } else {
                c *= sin(M_PI * ii / N)
            }
            t[i] = c
        }

        let M = Double(taps.count-1) / 2.0
        var h = [Double](count:taps.count, repeatedValue: 0)
        if taps.isOdd {
            for n in 0..<taps.count {
                var val = 0.0
                let x = 2 * M_PI * (Double(n) - M) / N
                for k in 1...Int(M) {
                    val += 2.0 * t[k] * sin(x*Double(k))
                }
                h[n] = val / N
            }
        } else {
            for n in 0..<taps.count {
                var val = t[taps.count/2] * sin(M_PI * (Double(n) - M))
                let x = 2 * M_PI * (Double(n) - M) / N
                for k in 1...taps.count/2-1 {
                    val += 2.0 * t[k] * sin(x*Double(k))
                }
                h[n] = val / N
            }
        }
        return h
    }

}


/// Searches for the maxima/minima of the error curve.
/// Return true when converged on a solution.
/// If more than r+1 extrema are found, it uses the following heuristic
/// (thanks Chris Hanson):
/// 1) Adjacent non-alternating extrema deleted first.
/// 2) If there are more than one excess extrema, delete the
///    one with the smallest error.  This will create a non-alternation
///    condition that is fixed by 1).
/// 3) If there is exactly one excess extremum, delete the smaller
///    of the first/last extremum
private func search(E:[Double], inout Ext:[Int]) -> Bool
{
    var found = 0
    var foundExt = [Int](count: Ext.count*2, repeatedValue: 0)

    // Check for extremum at 0.
    if E[0] > 0.0 && E[0] > E[1] || E[0] < 0.0 && E[0] < E[1] {
        foundExt[found] = 0
        found += 1
    }

    // Check for extrema inside dense grid
    for i in 1..<(E.count-1) {
        if E[i] >= E[i-1] && E[i] > E[i+1] && E[i] > 0.0 ||
            E[i] <= E[i-1] && E[i] < E[i+1] && E[i] < 0.0 {
                foundExt[found] = i
                found += 1
        }
    }

    // Check for extremum at 0.5
    let i = E.count-1
    let j = i - 1
    if E[i] > 0.0 && E[i] > E[j] || E[i] < 0.0 && E[i] < E[j] {
        foundExt[found] = j
        found += 1
    }

    // Remove extra extremals
    for extra in (found - Ext.count).stride(to: 0, by: -1) {
        var up:Bool
        if (E[foundExt[0]] > 0.0) {
            up = true   // first one is a maxima
        } else {
            up = false  // first one is a minima
        }

        var l = 0
        var alt = true
        for j in 1..<found
        {
            if (abs(E[foundExt[j]]) < abs(E[foundExt[l]])) {
                l = j         // new smallest error
            }
            if ((up) && (E[foundExt[j]] < 0.0)) {
                up = false    // switch to a minima
            } else if ((!up) && (E[foundExt[j]] > 0.0)) {
                up = true     // switch to a maxima
            } else {
                alt = false
                break         // Ooops, found two non-alternating
            }                 // extrema.  Delete smallest of them
        }  // if the loop finishes, all extrema are alternating

        // If there's only one extremal and all are alternating,
        if alt && extra == 1 {
            if abs(E[foundExt[found-1]]) < abs(E[foundExt[0]]) {
                l = foundExt[found-1]  // Delete last extremal
            } else {
                l = foundExt[0]        // Delete first extremal
            }
        }

        if (l < found) {
            for j in l..<found {
                foundExt[j] = foundExt[j+1]
            }
        }
        found -= 1
    }

    // Test for convergence and update Ext
    Ext[0] = foundExt[0]
    var minimum = abs(E[Ext[0]])
    var maximum = minimum
    for i in 1..<Ext.count {
        Ext[i] = foundExt[i]
        let current = abs(E[Ext[i]])
        minimum = min(minimum, current)
        maximum = max(maximum, current)
    }
    return (maximum-minimum)/maximum < 0.0001
}



let bands:[Double] = [0.0,0.05, 0.1,0.15, 0.18,0.25, 0.3,0.36, 0.41,0.5]
let desired:[Double] = [0,0,1,1,0,0,1,1,0,0]
let weights:[Double] = [10,1,3,1,20]
let filter:Filter = .BANDPASS

// Do a little work to wake up CPU and prime cache
for _ in 0...25 {
    firpm(104, bands: bands, desired: desired, weights: weights, filter: filter)
}

let startTime = CFAbsoluteTimeGetCurrent()
let h = firpm(104, bands: bands, desired: desired, weights: weights, filter: filter)
let timeElapsed = CFAbsoluteTimeGetCurrent() - startTime

print("Time elapsed: \(timeElapsed) s")
