set more off
clear all
set matsize 11000

/*

Project    : Resources, Small Sample Inference (Permutation Based)

Description: this .do file tests mean differences across two groups
			 using small sample inferences (permutation based)

Basics: declare the variable indicating the orbits on observed variables, O
                the sets of regressors you want to clean for, X
                the treatment indicator, Z
                the outcome, Y
		
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
gen y = rnormal(20,5)

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

// permute the treatmente variable
// set number of resamples
local B = 50

// clean the outcome variable for x
reg y x
predict yc, resid

// save current data as a tempfile
rename z z_0
tempfile data 
save   "`data'", replace

// permute the treatment indicator B times and store it
foreach num of numlist 1(1)`B'{
	preserve
	keep z o
	bsample, strata(o)
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
reg yc z_0 [iw=w]
matrix b = e(b)
local meandiff = b[1,1]
local p1 = 1 - normal(`meandiff')
local p2 = 2*(1 - normal(`meandiff'))
matrix naive = [`meandiff',`p1',`p2']
matrix colnames naive = meandiff onetailp twotailp
matrix rownames naive = naive

// capture the mean difference for all permutations
matrix pmeandiff = [.]
matrix colnames pmeandiff = meandiff
foreach num of numlist 1(1)`B'{
	reg yc z_`num' [iw=w]
	matrix b = e(b)
	local  meandiff_`num'  = b[1,1]
	matrix p`num' = [`meandiff_`num'']
	matrix colnames p`num' = meandiff
	matrix rownames p`num' = p`num'
	mat_rapp pmeandiff : pmeandiff p`num'
}
matrix pmeandiff = pmeandiff[2...,1]

// test
// export vector of means to data
preserve
clear
svmat pmeandiff
rename pmeandiff1 pmeandiff  
summ pmeandiff

// force to null distribution
replace pmeandiff = (pmeandiff - r(mean))
// absolute value for two tailed test
gen apmeandiff = abs(pmeandiff)

// gen meandiff as variable for easy comparison
gen meandiff = `meandiff'

// gen no reject indicators
// for absolute value
gen     rameandiff = 1
replace rameandiff = 0 if  apmeandiff < meandiff

// for level value
gen     rmeandiff = 1
replace rmeandiff = 0 if   pmeandiff < meandiff

// get p-values
// two-tailed
summ rameandiff 
local pap = r(mean)

// one-tailed
summ rmeandiff
local plp = r(mean)
matrix test = [naive,`plp',`pap']
matrix colnames test = meandiff naiveonetail naivetwotail ponetail ptwotail
matrix rownames test = test 
restore

// go back to initial data
drop yc
rename z_0 z
drop z_*

