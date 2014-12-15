set more off
clear all
set matsize 11000

/*
Project    : Resources, Single Axis Scatter Plot

Description: scatter plot various dependent variables over 
             x with a single vertical axis
             the format could be easily adapted for any twoway type plot
Basics:      
		
This version: 12/15/2014

This .do file: Jorge Luis Garcia
This project : CEHD

*/

// construct fake data
// comment in if you want to test
clear all
// generate data
set obs 10
gen id = _n  

// weights
gen w = rnormal(100,10)
replace w = w/100

// x simple ordered index
local b = 10
generate x = _n + 10

// y, outcome
foreach num of numlist 1(1)3 {
	gen y`num' = rnormal(20,5)
}

// x,y labels
global  xlabel x

// one axis scatter plot
#delimit
twoway (kdensity y1 [aw = w], lwidth(medthick) lpattern(solid) lcolor(gs0))
       (kdensity y2 [aw = w], lwidth(medthick) lpattern(solid) lcolor(gs8))
        , 
		  legend(label(1 $y1label) label(2 $y2label) label(3 $y3label) size(small))
		  xlabel(, grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14))
		  xtitle($xlabel) ytitle(Density, size(small))
		  graphregion(color(white)) plotregion(fcolor(white));
#delimit cr 
