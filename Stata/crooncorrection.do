set more off
clear all
set matsize 11000

/*

Project    : Resources, Croon's Measurement Error Correction

Description: this .do file corrects for measurement error the OLS estimates of a linear
			 regression after considering as regressors factor scores which are 
			 contaminated with measurement error

Basics: declare the factor scores, o
                the set of other controls, x
                the outcome, y
                estimated covariance of the factor scores, cmfs
                estimated b, first on factors then on x, last constant (column)
                output is bcorr
		
This version: 10/15/2014

This .do file: Jorge Luis Garcia
This project : CEHD

*/

// set seed (even if you do not simulate data, help for the resampling)
set seed 2

// construct fake data
// comment in if you want to test

clear all
// generate data
// 100 observations 
set obs 150
gen id = _n  

// w, weights
gen w = rnormal(100,10)
replace w = w/100

// x1, x2: regressors
gen x1 = rnormal(0,1)
local a = 0
local b = 4
gen x2 = floor((`b'-`a'+1)*runiform() + `a') 

// z1, z2, z3: items to factor analyze
local a = 0
local b = 1
gen z1 = floor((`b'-`a'+1)*runiform() + `a') 
gen z2 = rnormal(0,1)
gen z3 = rnormal(0,1)
gen z4 = rnormal(0,1)
gen z5 = rnormal(0,1)

// y, outcome
gen y = rnormal(20,5)

// perform factor analysis
factor z1 z2 z3 z4 z5, pf fac(2)
// rotate with quartimin
rotate, quartimin
// estimated covariance matrix of factor scores
matrix cmfs = e(r_Phi)
predict f1 f2

// obtain beta
reg y f1 f2 x1 x2
matrix b = e(b)'

// declare data
// individual identifier
global id id  

// weights
global w w

// items to factor analyze
global z z1 z2 z3

// regressors
gen ones = 1
global x x1 x2

// factors
global f f1 f2

// estimated covariance matrix of factor scores
global cmfs cmfs

// betas
global b b

// dimension of factors
local nf = 2
// dimensions of factors + controls + constant
local nx   = 5
local nxm1 = `nx' - 1 

// Croon's Correction
correlate $f $x [aw=$w]

// data covariance matrix
matrix cfx = r(C)

// covariance matrix with estimated covariance matrix for factors
matrix cefc = cfx

foreach row of numlist 1(1)`nf' {
	foreach col of numlist 1(1)`nf' {
		matrix cefc[`row',`col']  = $cmfs[`row',`col'] 
	}
}

// beta with no intercept
matrix bnc = $b[1..`nxm1',1]

// define correction matrix
matrix A = inv(cefc)*cfx

// correct
matrix bcorr = inv(A)*bnc

// add constant
matrix bcorr = [bcorr \ b[`nx',1]]
