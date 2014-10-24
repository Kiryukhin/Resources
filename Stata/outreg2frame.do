/*
Project    : Resources, Outreg2 Frame

Description: outreg frame to account for usual options
Basics:      est1 est2... are saves through 
 			 estimates store after each estimation
	
This version: 10/24/2014

This .do file: Jorge Luis Garcia
This project : CEHD
*/

#delimit
outreg2 [est1 est2 ... estn] using yourfile, replace tex(frag) 
		alpha(.01, .05, .10) sym (***, **, *) dec(5) par(se) r2
		ti(Your Title);
#delimit cr
