version 12.0
set more off
clear all
set matsize 11000

/*

Project : Resources

Description: this .do file creates and plots Laspeyres Decompositions
This version: February 1 , 2015

This .do file: Jorge L. Garcia
This project : CEHD
*/


// set environment variables
global erc:         env erc
global projects:    env projects
global klmshare:    env klmshare
global klmmexico:   env klmMexico
global googledrive: env googledrive

// input example data
insheet using test.csv, clear
gen w = 1

// declare weights
global w w
// output variable
global output output
// outcome indicator 
global outcome outcome
// treatment indicator
global treatment treat
// category indicator
global category  male
// decomposition partitions
global partition1 sb3 sb4
global partition2 cons5 cons6
global partition3 wt0

// outcome labels
label define outcomeslabel 1 "Idle at 30" 2 "HS Grad at 30"
label values $outcome outcomeslabel

// run Lasypeyres Decomposition
levelsof $outcome
local outcomes   = r(levels)
levelsof $category
local categories = r(levels)

foreach out of numlist `outcomes' {
	foreach cat of numlist `categories'{
		# delimit
		oaxaca $output $partition1 $partition2 $partition3 [iw=$w] if $category == `cat' & $outcome == `out', 
			   by($treatment) swap 
			   detail(partition1: $partition1, partition2: $partition2, partition3: $partition3); 
		# delimit cr
		
			   matrix estimation`cat'`out' = r(table)
			   
			   local   treat_cat`cat'_o`out' = round(estimation`cat'`out'[1,3], .01)
			   local pvtreat_cat`cat'_o`out' = round(estimation`cat'`out'[4,3], .01)
			   
			   local   endp1_cat`cat'_o`out' = estimation`cat'`out'[1,4]
			   local pvendp1_cat`cat'_o`out' = estimation`cat'`out'[4,4]
			   
			   local   endp2_cat`cat'_o`out' = estimation`cat'`out'[1,5]
			   local pvendp2_cat`cat'_o`out' = estimation`cat'`out'[4,5]
			   
			   local   endp3_cat`cat'_o`out' = estimation`cat'`out'[1,6]
			   local pvendp3_cat`cat'_o`out' = estimation`cat'`out'[4,6]
			   
			   local   techp1_cat`cat'_o`out' = estimation`cat'`out'[1,8]
			   local pvtechp1_cat`cat'_o`out' = estimation`cat'`out'[4,8]
			   
			   local   techp2_cat`cat'_o`out' = estimation`cat'`out'[1,9]
			   local pvtechp2_cat`cat'_o`out' = estimation`cat'`out'[4,9]
			   
			   local   techp3_cat`cat'_o`out' = estimation`cat'`out'[1,10]
			   local pvtechp3_cat`cat'_o`out' = estimation`cat'`out'[4,10]
				
			   local   resid_cat`cat'_o`out' = estimation`cat'`out'[1,11]
			   local pvresid_cat`cat'_o`out' = estimation`cat'`out'[4,11]			   
	}
}  


// partitions as proportions of the treatment effect
// allocate all relevant estimates in a matrix
# delimit
global allestimates treat treatpvalue endowp1 endowp1pval endowp2 endowp2pval 
	   endowp3 endowp3pval techp1 techp1pval techp2 techp2pval
	   techp3 techp3pval resid residpval $outcome $category;
# delimit cr 
matrix restimates = J(1,18,.)
matrix colnames restimates = $allestimates

foreach out of numlist `outcomes' {
	foreach cat of numlist `categories'{
		# delimit
		foreach val in endp1_cat endp2_cat endp3_cat techp1_cat 
				       techp2_cat techp3_cat resid_cat {;
		# delimit cr
			local `val'`cat'_o`out' = ``val'`cat'_o`out''/`treat_cat`cat'_o`out''
		
		# delimit
		matrix restimates_cat`cat'_o`out' = [`treat_cat`cat'_o`out'',`pvtreat_cat`cat'_o`out'',`endp1_cat`cat'_o`out'',
		              `pvendp1_cat`cat'_o`out'',`endp2_cat`cat'_o`out'',`pvendp2_cat`cat'_o`out'',
		              `endp3_cat`cat'_o`out'',`pvendp3_cat`cat'_o`out'',`techp1_cat`cat'_o`out'',
		              `pvtechp1_cat`cat'_o`out'',`techp2_cat`cat'_o`out'',
		              `pvtechp2_cat`cat'_o`out'',`techp3_cat`cat'_o`out'',
		              `pvtechp3_cat`cat'_o`out'',
		              `resid_cat`cat'_o`out'',
		              `pvresid_cat`cat'_o`out'',`out',`cat'];
		matrix colnames restimates_cat`cat'_o`out' = $allestimates; 
		# delimit cr
		}
		mat_rapp restimates : restimates restimates_cat`cat'_o`out'
	}
}
matrix restimates = restimates[2...,1...]

// plot
clear 
svmat restimates, names(col)
// duplicate variable to create different color when p-value is less than .1
foreach var in endowp1 endowp2 endowp3 techp1 techp2 techp3 resid {
	gen `var'_2 = `var'
	replace `var'    = . if `var'pval >  .10 | abs(`var') < .1
	replace `var'_2  = . if `var'pval <= .10 | abs(`var'_2) < .1
}

global partition1label Partition 1
global partition2label Partition 2
global partition3label Partition 3

foreach cat of numlist `categories' {
	# delimit
	global t1 `treat_cat`cat'_o1'; global p1 `pvtreat_cat`cat'_o1'; 
	global t2 `treat_cat`cat'_o2'; global p2 `pvtreat_cat`cat'_o2';		
	# delimit cr

	# delimit
	graph hbar endowp1 endowp1_2 techp1 techp1_2 endowp2 endowp2_2 techp2 
			   techp2_2 endowp3 endowp3_2 techp3 techp3_2 resid resid_2
	           
	           if $category == `cat', over(outcome, lab(labsize(small))) stack
	ylabel(, angle(h) glcolor(gs14) labsize(small)) 
	bar(1,  color(maroon*.5))            				bar(2, color(maroon)) 
	bar(3,  color(navy*.5))              				bar(4, color(navy)) 
	bar(5,  color(yellow*.5))            				bar(6, color(yellow))
	bar(7,  color(purple*.5))            				bar(8, color(purple))
	bar(9,  color(black*.5))             				bar(10, color(black))
	bar(11, color(dkorange*.5))                         bar(12, color(dkorange))
	bar(13, fcolor(black*.3) lcolor(black) lpattern(dash)) bar(14, fcolor(black*.3) lcolor(black) lpattern(solid))
	
	text(0 59 "Treatment Effect: ${t2} (p-value: ${p2})", place(e) size(vsmall))
	text(0 4  "Treatment Effect: ${t1} (p-value: ${p1})", place(e) size(vsmall))	 
	
	blabel(bar, position(center) format(%3.2f) color(gs14) size(tiny))
	legend(span rows(7) cols(2) label(1  "{&Delta} ${partition1label} Levels") label(2  "if p-value < .10") 
								label(3  "{&Delta} ${partition1label} Tech")   label(4  "if p-value < .10") 
								label(5  "{&Delta} ${partition2label} Levels") label(6  "if p-value < .10")
	                    		label(7  "{&Delta} ${partition2label} Tech")   label(8  "if p-value < .10")
	                    		label(9  "{&Delta} ${partition3label} Levels") label(10 "if p-value < .10")
	                    		label(11 "{&Delta} ${partition3label} Tech")   label(12 "if p-value < .10") 
	                    		label(13 "{&Delta} Residual")                  label(14 "if p-value < .10") 
	                    size(small)) 
    title(, size(small) position(12))
    graphregion(color(white)) plotregion(fcolor(white));
    # delimit cr
}
