

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
		
This version: 2/9/2015

This .do file: Jorge Luis Garcia, Andrés Hojman
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
label var y`num' "Label of Outcome `num'"
}
global y y1 y2 y3 y4 y5 y6 y7 y8 y9 y10
// declare data: input your dataset here

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
// declare outcomes
global y $y

*This section creates automatic labels for the rows (variables) based on original labels
*Careful: STATA matrices can only have 26 characters as rownames (including spaces)

local labels `" a "'
foreach var in $y{
local l`var' : variable label `var'
local l`var' `"`"`"`l`var''"' "'"'    // careful with modifying this! A space can make a difference
local labels `labels' `l`var''
di `" `labels' "' //just to check everything is going well
}
local labels: subinstr local labels `"a"' `""'
di `" `labels' "' //just to check everything is going well
global ylabel  `labels' 

// declare labels for outcomes only if you are not using the automatic labels
*# delimit
*global ylabel  ;
*# delimit cr;

*Reverse: this will change the one-sided tests later. Put here a list of "socially bad" variables
foreach var in {
local `var'reverse 1
}

// set number of resamples
local B = 10

* - - - - - - - - - - - - 
local nO: word count $y
// reverse and clean the outcome variable for x
foreach var of varlist $y $z {
	// reg `var' $z [iw=$w]
	// matrix b`var' = e(b)
	// local b`var' = b`var'[1,1]
	// gen reverse`var' = 1
	// replace reverse`var' = -1 if `b`var'' < 0
	// replace `var' = `var'*reverse`var'
	reg `var' $x //[iw=$w]
	predict `var'_cres, resid
}

// save current data as a tempfile
rename ${z}_cres z_0

drop id
gen id = _n //this is necessary for the id's in the original data to match with the new fake id's

tempfile data 
save   "`data'", replace

// permute the treatment indicator B times and store it
foreach num of numlist 1(1)`B'{
	preserve
	keep z_0 $o
	sample 99.9999, by($o)
	gen id = _n
	if "$o"!="" drop $o
	rename z_0 z_`num' 
	
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
matrix test = J(1,10,.)
matrix colnames test = Cmean Csd Tmean Tsd t diff naive1 naive2 perm1 perm2

foreach var of varlist $y { 
	reg `var'_cres z_0 //[iw=$w]
	matrix b = e(b)
	local meandiff`var' = b[1,1]
	di "Diff: `meandiff`var''"
	matrix V = e(V)
	local meandiff`var't = (`meandiff`var'')/sqrt(V[1,1])

	egen miss`var'=rowmiss(`var'_cres z_0) //ignores weighting=0
	
	sum `var' if miss`var'==0 & z==0
	local Csd`var'  	= r(sd)
	local Cmean`var'  	= r(mean)

	sum `var' if miss`var'==0 & z==1
	local Tsd`var'		= r(sd)
	local Tmean`var'	= r(mean)

	if "``var'reverse'"=="1"	local p1`var' = normal(`meandiff`var't')
	else						local p1`var' = 1 - normal(`meandiff`var't')
	local p2`var' = 2*(1 - normal(abs(`meandiff`var't')))
	matrix naive`var' = [`Cmean`var'',`Csd`var'',`Tmean`var'',`Tsd`var'',`meandiff`var'',`meandiff`var't',`p1`var'',`p2`var'']
	matrix colnames naive`var' = Cmean Csd Tmean Tsd diff t naive1 naive2
	matrix rownames naive`var' = naive
	}

foreach var of varlist $y {
	// capture the mean difference for all permutations
	matrix pmeandiff`var' = [.]
	matrix colnames pmeandiff`var' = meandiff
	foreach num of numlist 1(1)`B'{
		reg `var'_cres z_`num' //[iw=$w]
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

	// absolute value for two tailed test
	gen apmeandiff`var' = abs(pmeandiff`var')

	// gen meandiff as variable for easy comparison
	gen meandiff`var' = `meandiff`var''

	// gen no reject indicators
	// for absolute value (2-tails)
	gen     rameandiff`var' = .
	replace rameandiff`var' = 1 if  apmeandiff`var' > abs(meandiff`var')
	replace rameandiff`var' = 0 if  apmeandiff`var' < abs(meandiff`var')

	// for level value (1-tail)
	gen     rmeandiff`var' = .
	if "``var'reverse'"=="1"{
	replace rmeandiff`var' = 1 if   pmeandiff`var' < meandiff`var'
	replace rmeandiff`var' = 0 if   pmeandiff`var' > meandiff`var'
	}
	else{
	replace rmeandiff`var' = 1 if   pmeandiff`var' > meandiff`var'
	replace rmeandiff`var' = 0 if   pmeandiff`var' < meandiff`var'
	}


	// get p-values
	// two-tailed
	summ rameandiff`var' 
	local pap`var' = r(mean)

	// one-tailed
	summ rmeandiff`var'
	local plp`var' = r(mean)
	
	matrix test`var' = [naive`var',`plp`var'',`pap`var'']
	matrix colnames test`var' = Cmean Csd Tmean Tsd diff t naive1 naive2 perm1 perm2
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
svmat test, names(col)
gen on = _n
sort naive1
mkmat *, mat(testonesidedp)
restore

// two tails
preserve
clear
svmat test, names(col)
gen on = _n
sort naive2
mkmat *, mat(testtwosidedp)
restore


// stepdown correction for permutation p-value  //
// - - - - - - - - - - - - - - - - - - - - - - -//

// all permutation estimates in a matrix
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

matrix sdtt = [.]
matrix colnames sdtt = sdtt
matrix sdot = [.]
matrix colnames sdot = sdot

foreach num of numlist 1(1)`nO' {
	capture drop pmeandiff ipmeandiff_two ipmeandiff_one
	egen pmeandiff = rowmax(pmeandiff`num'-pmeandiff`nO')
	
	matrix list testtwosidedp
	asd	

	
	// two tails
	gen     ipmeandiff_two = 1
	replace ipmeandiff_two = 0 if  abs(pmeandiff) < abs(testtwosidedp[`num',6]) 
	summ ipmeandiff_two 
	matrix sd`num'tt = r(mean)
	matrix colnames sd`num'tt = sdtt
	mat_rapp sdtt : sdtt sd`num'tt
	
	
	// one tail
	gen		ipmeandiff_one = 1
	replace ipmeandiff_one = 0 if pmeandiff < testonesidedp[`num',6] 
	summ ipmeandiff_one
	matrix sd`num'ot = r(mean)
	matrix colnames sd`num'ot = sdot
	mat_rapp sdot : sdot sd`num'ot
	
}

matrix sdtt = sdtt[2...,1...]
matrix sdot = sdot[2...,1...]
restore

// from test merge meandiff and naive pvals
//[AH: I think this just puts test into data and could be done before, but gen on is necessary]
preserve
clear
svmat test, names(col)
gen on = _n
tempfile test
save "`test'", replace
restore
*use "`test'", clear
*use `test', clear
*asd

// from testonesidedp merge one tail permutation and sd
// columns to keep: on sd2
preserve 
clear
matrix testonesidedp = [testonesidedp,sdot]
svmat testonesidedp 
/*
keep testonesidedp4 testonesidedp6 testonesidedp7
rename testonesidedp4 permonetail
rename testonesidedp6 on
rename testonesidedp7 peronetailsd
*/
//AH: trying with this:
keep testonesidedp10 testonesidedp11 
rename testonesidedp10 on
rename testonesidedp11 peronetailsd

//AH: this was from before: 
tempfile testoneper
save "`testoneper'", replace
restore

// from testtwosidedp merge one tail permutation and sd
// columns to keep: perm2 on sd2
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
matrix list testfmatrix
asd
global edulabel   `"`"Years High School "' "' `"`"Years Community college "' "'  `"`"Years High School "' "' `"`"Years Community college "' "'  `"`"Years High School "' "' 
macro list
di `" $edulabel "'
di `" $ylabel "'
matrix rownames testfmatrix =  `" $edulabel "'
matrix list testfmatrix
matrix rownames testfmatrix =  `" $ylabel "'
matrix list testfmatrix
restore

// output matrix
matrix testfmatrix = testfmatrix[1...,1..5]

outtable using "yourtableintext", ///
 mat(testfmatrix) replace nobox center f(%9.3f)


// go back to initial data
drop *_cres
// drop reverse*
rename z_0 z_cres
drop z_*
