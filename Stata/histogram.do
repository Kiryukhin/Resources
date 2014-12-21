set more off
clear all
set matsize 11000

/*
Project    : Resources, Histogram

Description: 
Basics:      
		
This version: 12/19/2014

This .do file: Sneha Elango
This project : CEHD

*/

/*
// construct fake data
// comment in if you want to test
clear all
// generate data
set obs 10
gen id = _n  

// weights
gen w = rnormal(100,10)
replace w = w/100

// y, outcome
foreach num of numlist 1(1)3 {
	gen y`num' = floor(rnormal(20,5))
}

// x,y labels
global y1label y1
global y2label y2
global y3label y3*/

# delimit
twoway (histogram y1, start(1) discrete fraction color(gs10)  barwidth(.75) yline(.1, lcolor(gs5) lpattern(dash)))
	   (histogram y2, start(1) discrete fraction fcolor(none) barwidth(.75) lcolor(black)),
	   legend(label(1 ${label2}) label(2 ${label3}))
	   xtitle(Y) ytitle(Fractiion)
	   xlabel(, grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14))
	   graphregion(color(white)) plotregion(fcolor(white));
#delimit cr
	graph export histogram'.eps, replace