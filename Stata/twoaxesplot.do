set more off
clear all
set matsize 11000

/*
Project    : Resources, Two Axes Scatter Plot

Description: scatter plot two dependent variables over x with 
             two different vertical axes
             the format could be easily adapted for any twoway type plot

Basics:      
		
This version: 10/24/2014

This .do file: Jorge Luis Garcia
This project : CEHD
*/

// two axes plot
#delimit
twoway (scatter y1 x, msymbol(triangle) mfcolor (gs4) mlcolor(gs4) connect(l) lwidth(medthick) lpattern(solid) lcolor(gs4) yaxis(1))
       (scatter y2 x, msymbol(square)   mfcolor (gs8) mlcolor(gs8) connect(l) lwidth(medthick) lpattern(dash)  lcolor(gs8) yaxis(2))
        , 
		  legend(label(1 y1) label(2 y2) size(small))
		  xlabel(, grid glcolor(gs14)) ylabel(, angle(h) glcolor(gs14))
		  xtitle(x) ytitle(y1, axis(1)) ytitle(y2, axis(2))
		  graphregion(color(white)) plotregion(fcolor(white));
#delimit cr 
graph export twoaxesscatter.eps, replace
