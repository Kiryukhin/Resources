

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
/*
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
*/
// declare data

// individual identifier
global id id  

// weights
*global w w

// orbits
global o 

// treatment indicator
global z treat

// continuous x
global x f_home3 m_work3 sb3 ses

// estimation starts here
// declare outcomes
global y 	pbi_acad6 pbi_acad7 pbi_acad8 pbi_acad9 pbi_acad69 ///
			pbi_cond6 pbi_cond7 pbi_cond8 pbi_cond9 pbi_cond69 /// 
			pbi_emo6 pbi_emo7 pbi_emo8 pbi_emo9 pbi_emo69 ///
			pbi_dep6 pbi_dep7 pbi_dep8 pbi_dep9 pbi_dep69 ///
			pbi_beh6 pbi_beh7 pbi_beh8 pbi_beh9  pbi_beh69 ///
			yrs_emo6 yrs_emo7 yrs_emo8 yrs_emo9 yrs_emo69 ///
			yrs_acad6 yrs_acad7 yrs_acad8 yrs_acad9 yrs_acad69 ///
			yrs_vrb6 yrs_vrb7 yrs_vrb8 yrs_vrb9 yrs_vrb69 ///
			yrs_mom6 yrs_mom7 yrs_mom8 yrs_mom9 yrs_mom69 ///
			yrs_soc6 yrs_soc7 yrs_soc8 yrs_soc9 yrs_soc69

local labels `" a "'
foreach var in $y{
local l`var' : variable label `var'
local l`var' `"`"`"`l`var''"' "'"'
local labels `labels' `l`var''
di `" `labels' "'
}
local labels: subinstr local labels `"a"' `""'
di `" `labels' "'

// declare labels for outcomes
# delimit
*global ylabel `" `labels' "' ;
global ylabel  `labels'  ;
# delimit cr;

*Reverse
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
gen id = _n


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
matrix test = J(1,5,.)
matrix colnames test = diff naive1 naive2 perm1 perm2

foreach var of varlist $y { 
	reg `var'_cres z_0 //[iw=$w]
	matrix b = e(b)
	local meandiff`var'  = b[1,1]
	di "`meandiff`var''"

	matrix V = e(V)
	local meandiff`var't = (`meandiff`var'')/sqrt(V[1,1])
	
	
	if "``var'reverse'"=="1"	local p1`var' = normal(`meandiff`var't')
	else						local p1`var' = 1 - normal(`meandiff`var't')
	local p2`var' = 2*(1 - normal(abs(`meandiff`var't')))
	matrix naive`var' = [`meandiff`var'',`p1`var'',`p2`var'']
	matrix colnames naive`var' = diff naive1 naive2
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
	matrix colnames test`var' = diff naive1 naive2 perm1 perm2
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

outtable using "/Users/andreshojman/Desktop/RodrigoTE/RodrigoTE1", ///
 mat(testfmatrix) replace nobox center f(%9.3f)


// go back to initial data
drop *_cres
// drop reverse*
rename z_0 z_cres
drop z_*
