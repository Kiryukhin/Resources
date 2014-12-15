set more off
clear all
set matsize 11000

/*

Project    : Resources, Quantiles Comparison Plots

Description: this .do file generates plots establishing how one distribution fits other
			 distribution

Basics: declare the variable to be analyzes, f, the category, c, and 
		allow for three categories: 1,2,3. 1 is the main category. 

		
This version: 12/01/2014

This .do file: Jorge Luis Garcia
This project : CEHD

*/


// construct fake data
// comment in if you want to test
clear all
set obs 1000
gen id = _n  

// generate relevant variables
gen f = rnormal(100,15)

// generate three categories
generate cat = floor((3-1+1)*runiform() + 1)


// starts here when not test
// declare variables
// variable to be analyzed
global f f
// category analysis
global cat cat

// label main category
global label1 W 

// label reminder categories
global label2 B
global label3 H

// declare quantiles
local   q = 10
local qm1 = `q' - 1
// obtain quantiles of the main distribution
_pctile $f if $cat == 1, nq(`q')

foreach num of numlist 1(1)`qm1' {
	gen fqtile`num' = r(r`num')
}

foreach numy of numlist 2 3 {
	foreach num of numlist `qm1'(-1)1 {
		gen qindicator`numy'`num' = 1 if $f != . & $cat == `numy'
		replace qindicator`numy'`num' = `num' + 1 if $f >= fqtile`num' & $cat == `numy' & $f != .
	}
	egen qindicator`numy' = rowmax(qindicator`numy'*)
	drop qindicator`numy'`qm1'-qindicator`numy'1
}


# delimit
twoway (histogram qindicator2, start(1) discrete fraction color(gs10)  barwidth(.75) yline(.1, lcolor(gs5) lpattern(dash)))
	   (histogram qindicator3, start(1) discrete fraction fcolor(none) barwidth(.75) lcolor(black)),
	   legend(label(1 ${label2}) label(2 ${label3}))
	   xtitle(Quantiles in the ${label1} Distribution) ytitle(Fractiion)
	   xlabel(1[1]`q', grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14))
	   graphregion(color(white)) plotregion(fcolor(white));
delimit cr
