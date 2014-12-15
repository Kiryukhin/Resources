set more off
clear all
set matsize 11000

/*

Project    : Resources, Small Sample Inference (Permutation Based)

Description: this .do file tests mean differences across two groups
			 using small sample inference (permutation based)

Basics: declare the variable indicating the orbits on observed variables, o
				the weights, w
                the sets of regressors you want to clean for, x
                the treatment indicator, z
                the outcomes, y
		
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

// weights
gen w = rnormal(100,10)
replace w = w/100

// orbits, uniform [a,b]
local a = 0
local b = 4
generate o = floor((`b'-`a'+1)*runiform() + `a')

// treatment, uniform [0,1] 
local a = 0
local b = 1
generate z = floor((`b'-`a'+1)*runiform() + `a')

// x, normal standard
gen x=rnormal(0,1)

// y, outcome
foreach num of numlist 1(1)10 {
	gen y`num' = rnormal(20,5)
}

// declare data
// individual identifier
global id id  

// weights
global w w

// orbits
global o o 

// treatment indicator
global z z

// continuous x
global x x

// declare number of outcomes
local nO = 10
// declare outcomes
global y y1-y10

// permute the treatmente variable
// set number of resamples
local B = 50

// clean the outcome variable for x
foreach var of varlist $y {
	reg `var' $x [iw=$w]
	predict `var'_cres, resid
}

// save current data as a tempfile
rename z z_0
tempfile data 
save   "`data'", replace

// permute the treatment indicator B times and store it
foreach num of numlist 1(1)`B'{
	preserve
	keep z o
	sample 99.9999, by(o)
	gen id = _n
	drop o
	rename z z_`num' 
	
	tempfile data_`num'
	save   "`data_`num''", replace
	restore
}

// merge each permuted treatment indicator
foreach num of numlist 1(1)`B'{
	merge 1:1 id using "`data_`num''"
	keep if _merge == 3
	drop _merge
}

// capture mean difference and naive p-value of one and two tailes tests
matrix test = J(1,5,.)
matrix colnames test = meandiff naiveonetail naivetwotail ponetail ptwotail

foreach var of varlist $y { 
	reg `var'_cres z_0 [iw=$w]
	matrix b = e(b)
	local meandiff`var' = b[1,1]
	local p1`var' = 1 - normal(abs(`meandiff`var''))
	local p2`var' = 2*(1 - normal(abs(`meandiff`var'')))
	matrix naive`var' = [`meandiff`var'',`p1`var'',`p2`var'']
	matrix colnames naive`var' = meandiff onetailp twotailp
	matrix rownames naive`var' = naive
}

foreach var of varlist $y {
	// capture the mean difference for all permutations
	matrix pmeandiff`var' = [.]
	matrix colnames pmeandiff`var' = meandiff
	foreach num of numlist 1(1)`B'{
		reg `var'_cres z_`num' [iw=$w]
		matrix b = e(b)
		local  meandiff_`num'  = b[1,1]
		matrix p`num'`var' = [`meandiff_`num'']
		matrix colnames p`num'`var' = meandiff
		matrix rownames p`num'`var' = p`num'`var'
		mat_rapp pmeandiff`var' : pmeandiff`var' p`num'`var'
	}
	matrix pmeandiff`var' = pmeandiff`var'[2...,1]

	// test
	// export vector of means to data
	preserve
	clear
	svmat pmeandiff`var'
	rename pmeandiff`var'1 pmeandiff`var'  
	summ pmeandiff`var'

	// force to null distribution
	replace pmeandiff`var' = (pmeandiff`var' - r(mean))
	// absolute value for two tailed test
	gen apmeandiff`var' = abs(pmeandiff`var')

	// gen meandiff as variable for easy comparison
	gen meandiff`var' = `meandiff`var''

	// gen no reject indicators
	// for absolute value
	gen     rameandiff`var' = 1
	replace rameandiff`var' = 0 if  apmeandiff`var' < abs(meandiff`var')

	// for level value
	gen     rmeandiff`var' = 1
	replace rmeandiff`var' = 0 if   pmeandiff`var' < meandiff`var'

	// get p-values
	// two-tailed
	summ rameandiff`var' 
	local pap`var' = r(mean)

	// one-tailed
	summ rmeandiff`var'
	local plp`var' = r(mean)
	matrix test`var' = [naive`var',`plp`var'',`pap`var'']
	matrix colnames test`var' = meandiff naiveonetail naivetwotail ponetail ptwotail
	matrix rownames test`var' = test`var'
	
	// append to test matrix
	mat_rapp test : test test`var' 
	restore
}

// clean matrix of first starting row
matrix test = test[2...,1...]


// go back to initial data
drop *_cres
rename z_0 z
drop z_*
