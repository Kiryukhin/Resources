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
                the outcomes, y --careful! the code reverses this outcomes
                *I recommend use this through a temporary file, not the actual
                 data file
		
This version: 12/19/2014

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

// estimation starts here
// declare number of outcomes
local nO = 10
// declare outcomes
global y y1-y10
// declare labels for outcomes
# delimit
global ylabel y1 y2 y3 y4 y5
              y6 y7 y8 y9 y10; 
# delimit cr

// set number of resamples
local B = 1000

// reverse and clean the outcome variable for x
foreach var of varlist $y {
	reg `var' $z [iw=$w]
	matrix b`var' = e(b)
	local b`var' = b`var'[1,1]
	gen reverse`var' = 1
	replace reverse`var' = -1 if `b`var'' < 0
	replace `var' = `var'*reverse`var'
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
// create matrices sorted according to p-values
// one tail
preserve
clear
svmat test
gen on = _n
sort test2
mkmat *, mat(testonesidedp)
restore

// two tails
preserve
clear
svmat test
gen on = _n
sort test3
mkmat *, mat(testtwosidedp)
restore

// stepdown correction for permutation p-value
// all permutation estmiates in a matrix
matrix pmeandiff = J(`B',1,.)
foreach var of varlist $y {
	matrix pmeandiff = [pmeandiff,pmeandiff`var']
}

// sort outcomes according to p-value size
matrix pmeandiff = pmeandiff[1...,2...]'
matrix pmeandiff = [test,pmeandiff]
preserve
clear
svmat pmeandiff
sort  pmeandiff2
// dop all test info
drop pmeandiff1-pmeandiff5

// to data
mkmat *, matrix(pmeandiff)
clear
matrix pmeandiff = pmeandiff'
svmat pmeandiff

// impose the null
foreach num of numlist 1(1)`nO'{
	summ pmeandiff`num'
	replace pmeandiff`num' = (pmeandiff`num' - r(mean))
}


matrix sdtt = [.]
matrix colnames sdtt = sdtt
matrix sdot = [.]
matrix colnames sdot = sdot

foreach num of numlist 1(1)`nO' {
	capture drop pmeandiff ipmeandiff_two ipmeandiff_one
	egen pmeandiff = rowmax(pmeandiff`num'-pmeandiff`nO')
	
	// two tails
	gen     ipmeandiff_two = 1
	replace ipmeandiff_two = 0 if  abs(pmeandiff) < abs(testtwosidedp[`num',1]) 
	summ ipmeandiff_two 
	matrix sd`num'tt = r(mean)
	matrix colnames sd`num'tt = sdtt
	mat_rapp sdtt : sdtt sd`num'tt
	
	
	// one tail
	gen		ipmeandiff_one = 1
	replace ipmeandiff_one = 0 if pmeandiff < testonesidedp[`num',1] 
	summ ipmeandiff_one
	matrix sd`num'ot = r(mean)
	matrix colnames sd`num'ot = sdot
	mat_rapp sdot : sdot sd`num'ot
	
}

matrix sdtt = sdtt[2...,1...]
matrix sdot = sdot[2...,1...]
restore

// from test merge meandiff and naive pvals
preserve
clear
svmat test
keep test1-test3
gen on = _n
rename test1 meandiff
rename test2 naiveonetail
rename test3 naivetwotail
tempfile test
save "`test'", replace
restore


// from testonesidedp merge one tail permutation and sd
preserve 
clear
matrix testonesidedp = [testonesidedp,sdot]
svmat testonesidedp 
keep testonesidedp4 testonesidedp6 testonesidedp7
rename testonesidedp4 permonetail
rename testonesidedp6 on
rename testonesidedp7 peronetailsd
tempfile testoneper
save "`testoneper'", replace
restore

// from testtwosidedp merge one tail permutation and sd
preserve 
clear
matrix testtwosidedp = [testtwosidedp,sdtt]
svmat testtwosidedp 
keep testtwosidedp5 testtwosidedp6 testtwosidedp7
rename testtwosidedp5 permtwotail
rename testtwosidedp6 on
rename testtwosidedp7 pertwotailsd
tempfile testtwoper
save "`testtwoper'", replace
restore

// generate table with all the p-values
preserve
use "`test'", clear
merge 1:1 on using "`testoneper'"
tab _merge
drop if _merge != 3
drop _merge
merge 1:1 on using "`testtwoper'"
tab _merge
drop if _merge != 3
drop _merge on
aorder
mkmat *, mat(testfmatrix)
matrix rownames testfmatrix = $ylabel  
restore

// output matrix
/*
#delimit
outtable using yourfile, 
mat(testfmatrix) replace nobox center f(%9.3f);
#delimit cr
*/


// go back to initial data
drop *_cres
drop reverse*
rename z_0 z
drop z_*
