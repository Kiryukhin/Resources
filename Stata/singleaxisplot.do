set more off
clear all
set matsize 11000

/*
Project    : Resources, Single Axis Scatter Plot

Description: scatter plot various dependent variables over 
             x with a single vertical axis
             the format could be easily adapted for any twoway type plot
Basics:      
		
This version: 10/24/2014

This .do file: Jorge Luis Garcia
This project : CEHD

*/


// plot
#delimit
twoway (scatter y1 x, msymbol(circle) mfcolor (gs0) mlcolor(gs0) connect(l) lwidth(medthick) lpattern(solid) lcolor(gs0))
       (scatter y2 x, msymbol(triangle) mfcolor (gs4) mlcolor(gs4) connect(l) lwidth(medthick) lpattern(solid) lcolor(gs4))
       (scatter y3 x, msymbol(square) mfcolor (gs8) mlcolor(gs8) connect(l) lwidth(medthick) lpattern(dash)  lcolor(gs8))
        , 
		  legend(label(1 y1) label(2 y2) label(3 y3) size(small))
		  xlabel(, grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14))
		  xtitle(y) ytitle(y, size(small))
		  graphregion(color(white)) plotregion(fcolor(white));
#delimit cr 
graph export singleaxisscatter.eps, replace
