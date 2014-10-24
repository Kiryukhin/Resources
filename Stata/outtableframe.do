/*
Project    : Resources, Outtable Frame

Description: outtable frame to output a table
			 from a matrix in Stata
Basics:      
		
This version: 10/24/2014

This .do file: Jorge Luis Garcia
This project : CEHD
*/

#delimit
outtable using yourfile, 
mat(yourmatrux) replace nobox center f(%9.3f);
#delimit cr
