set more off
clear all
set matsize 11000

/*

Project    : Resources, Quantiles Comparison Plots

Description: this .do file generates plots establishing how one distribution fits other
			 distribution

Basics: declare the measures Z, and the individual indetifier id
        obtain a temporal file called F to be merged with the original
		data set through the provided id

		
This version: 12/01/2014

This .do file: Jorge Luis Garcia
This project : CEHD

*/
/*
// construct fake data
// comment in if you want to test

clear all
// 100 observations 
set obs 100
gen id = _n  

// generate relevant variables
gen w = rnormal(100,15)
gen b = rnormal(100,15)
gen h = rnormal(100,15)
*/ 

// declare variables
// main distribution
global w w
// main distribution label
global wlabel W 
// input distributions to analyze
global bh b h
// label input distributions
global blabel B
global hlabel H

// create plot
// declare quantiles
local   q = 10
local qm1 = `q' - 1
// obtain quantiles of the main distribution
_pctile $w, nq(`q')

foreach num of numlist 1(1)`qm1' {
	gen wqtile`num' = r(r`num')
}

foreach var of varlist $bh {
	foreach num of numlist `qm1'(-1)1 {
		gen qindicator`var'`num' = 1 if `var' != .
		replace qindicator`var'`num' = `num' + 1 if `var' >= wqtile`num' 
	}
	egen qindicator`var' = rowmax(qindicator`var'*)
	drop qindicator`var'`qm1'-qindicator`var'1
}


# delimit
twoway (histogram qindicatorb, start(1) discrete fraction color(gs10)  barwidth(.75) yline(.1, lcolor(gs5) lpattern(dash)))
	   (histogram qindicatorh, start(1) discrete fraction fcolor(none) barwidth(.75) lcolor(black)),
	   legend(label(1 ${blabel}) label(2 ${hlabel}))
	   xtitle(Quantiles in the ${wlabel} Distribution) ytitle(Fractiion)
	   xlabel(1[1]`q', grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14))
	   graphregion(color(white)) plotregion(fcolor(white));
delimit cr
