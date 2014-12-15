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
global y1label y1
global y2label y2
global y3label y3
global  xlabel x

// one axis scatter plot
#delimit
twoway (scatter y1 x [aw = w], msymbol(circle)   mfcolor (gs0) mlcolor(gs0) connect(l) lwidth(medthick) lpattern(solid) lcolor(gs0))
       (scatter y2 x [aw = w], msymbol(triangle) mfcolor (gs4) mlcolor(gs4) connect(l) lwidth(medthick) lpattern(solid) lcolor(gs4))
       (scatter y3 x [aw = w], msymbol(square)   mfcolor (gs8) mlcolor(gs8) connect(l) lwidth(medthick) lpattern(dash)  lcolor(gs8))
        , 
		  legend(label(1 $y1label) label(2 $y2label) label(3 $y3label) size(small))
		  xlabel(, grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14))
		  xtitle($xlabel) ytitle(, size(small))
		  graphregion(color(white)) plotregion(fcolor(white));
#delimit cr 
