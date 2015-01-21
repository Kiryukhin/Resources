set more off
clear all
set matsize 11000

/*
Project    : Resources, Segmented bar graph

Description: 
Basics:      
		
This version: 12/19/2014

This .do file: Sneha Elango
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

// category
gen cat = floor(rnormal()*3)
// y, outcome
gen total = rnormal(100,10)
gen y1 = total-rnormal(70,20)
gen y2 = total-y1-rnormal(30,5)
gen y3 = total-y1-y2

gen py1=y1/total
gen py2=y2/total
gen py3=y3/total

// x,y labels
global y1label y1
global y2label y2
global y3label y3

//stacked bar graph

# delimit 
 graph bar py1 py2 py3, over(cat) stack blabel(cat1 cat2) bar(1, color(gs3)) bar(2, color(gs6)) bar(3, color(gs9)) bargap(0)
		
	   legend()
	   ytitle(Percentage)
	   ylabel(, angle(h) glcolor(gs14))
	   graphregion(color(white)) plotregion(fcolor(white));
#delimit cr
	// graph export stackedbar.eps, replace
	

	
// horizontal bargraph

# delimit 
 graph hbar py1 py2 py3, over(cat)  blabel(cat1 cat2) bar(1, color(gs3)) bar(2, color(gs6)) bar(3, color(gs9)) bargap(0)
	   legend()
	   ytitle(Percentage)
	   ylabel(, angle(h) glcolor(gs14))
	   graphregion(color(white)) plotregion(fcolor(white));
#delimit cr
	// graph export hbar.eps, replace

	
# delimit 
 twoway( bar py2 id, vertical bargap(2) barwidth(.75) fcolor(gs4) lcolor(gs16)),

	   legend(off)
		  xlabel(, valuelabel labsize(vsmall) angle(45))
		  ylabel(, angle(h) glcolor(gs14)) 
		  xtitle("") ytitle(Percentage, size(small))
	   ytitle(Percentage)
	   ylabel(, angle(h) glcolor(gs14))
	   graphregion(color(white)) plotregion(fcolor(white));
#delimit cr
	// graph export twowaybar.eps, replace
	
	
