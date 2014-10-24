set more off
clear all
set matsize 11000

/*

Project    : Resources, Diagonal Factor Analysis

Description: this .do file produces diagonal factor analysis for an arbitrary full rank
             measurement system (the code crashes if the measurement system is not full rank)

Basics: declare the measures Z, and the individual indetifier id
        obtain a temporal file called F to be merged with the original
		data set through the provided id

		
This version: 07/01/2014

This .do file: Jorge Luis Garcia
This project : CEHD

*/

// construct fake data
// comment in if you want to test
/*
clear all
// 100 observations 
set obs 100
gen id = _n  

// weights
*gen w = floor((100-1+1)*runiform() + 1)
gen w = rnormal(100,10)

// x1, uniform [a,b]
local a = 0
local b = 1
gen x1 = `a'+(`b'-`a')*runiform()

// x2, uniform [c,d] 
local c = 2
local d = 3
gen x2 = `c'+(`d'-`c')*runiform()

// epsilon, normal standard
gen x3=rnormal(0,1)

// declare data
global Z x1 x2 x3  

// declare individual identifier
global id id  

// declare weights
global w w
*/

// diagonal factor analysis
// declare measurement system as a global variable list Z
// declare individual identifier as a global variable list id
// declare weights as a global variable w (ones if no weights)

// observations
qui des $Z
local N = r(N) 

// get rid of missings to avoid matrix inversion issues
foreach var of varlist $Z {
	drop if `var' == .
}

// declare variables and standardize

// standardize
foreach var of varlist $Z {
	summ `var' [iw = $w] 
	scalar mean`var' = r(mean)
	scalar sd`var' = r(sd)
	
	replace `var' = (`var' - mean`var')/sd`var'
}

// standardized data to matrix	
mkmat $Z, matrix(Z)
// identifier to matrix
mkmat $id, matrix(id)

// destandardize to leave original data untouched	
foreach var of varlist $Z {
	summ `var' [iw = $w] 
	replace `var' = `var'*sd`var'+mean`var'
}

// clear and start
// correlation matrix
correlate $Z [aw = $w]
matrix C = r(C)

// drop weights
drop w

// stop if measures are not linearly indepdent
matrix list C
scalar cC = colsof(C) 
matrank C, r(rC)

di cC " number of measures"
di rC " are linearly independent"
di cC " should be identical to " rC " for the code to finish with success"

// calculate eigenvalues and define the matrix of
// absolute values of eigenvalues (all R given pos def)
matrix symeigen vector eig = C
matrix eig = eig'
drop *

svmat eig
rename eig1 eig
replace eig = -eig if eig < 0 
gsort -eig
gen neig = _n 

summ neig 
local maxe = r(max)
// scree plot
#delimit
twoway (scatter eig neig, msymbol(circle) mfcolor(white) mlcolor(gs0) mlpattern(solid) msize(medium) connect(l) lcolor(gs0) lwidth(medium) 
                          yline(1, lcolor(gs8) lpattern(dash)))
	, ytitle(Eigenvalue (Absolute Value), size(small)) xtitle(Factor, size(small))
	  xlabel(1[1]`maxe', grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14)) 
      graphregion(color(white)) 
	  plotregion(fcolor(white));
#delimit cr 

// define number of factors
drop if eig < 1
des

// account for cases with a singe measure
local NE = max(r(N),1)
di "number of factors is`NE'"

// square of correlation matrix
clear
svmat C
des
local K = r(k)

matrix P = J(`K',1,0)

forvalues k = 1(1)`NE'{
	matrix C_2 = J(`K',`K',0)

	forvalues l = 1(1)`K' {
		forvalues j = 1(1)`K' {
			matrix C_2[`l',`j']= C[`l',`j']*C[`l',`j']
		}
	}

	// add square sum of columns (by row) to correlation matrix and sort
	matrix col2sum = [.] 
	forvalues i = 1(1)`K' {
		local colsum = 0
		forvalues j = 1(1)`K' {
			local colsum = `colsum' + C_2[`j',`i']
		}
		
		matrix col2sum_`i' = `colsum'
		matrix colnames col2sum_`i' = c`i' 
		mat_capp col2sum : col2sum col2sum_`i'
	
	}

	matrix col2sum = col2sum[1...,2...]
	matrix C = [C\col2sum]
	matrix CT = C'
	
	clear
	svmat CT
	
	// keep the row (which is the column since I transpose)
	// with max c2
	local l = `K' + 1
	summ CT`l'
	local maxc = r(max) - .00001                           
	drop if CT`l' < `maxc'
	drop CT`l'
	mkmat *, matrix(p)
	matrix p = p'
	
	matrix C = C[1..`K',1...]
	matrix C = C - p*p'
	
	matrix P = [P,p]
}	

// cut the matrix of loadings
matrix P = P[1...,2...]

// construct factors and save in a temporal file
// together with individual identifier
matrix F = Z*P*(inv(P'*P))
svmat F 
keep F*
gen n = _n

tempfile F 
save "`F'", replace

clear
svmat id
rename id1 $id
gen n = _n

merge 1:1 n using "`F'"
tab _merge
drop n _merge

save "`F'", replace
