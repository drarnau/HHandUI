/*
Choi & Valladares-Esteban (2018)

This codes reads Monte-Carlo simulation data from the model and generates results
*/

cls
clear all
set more off, perm
set trace off

// Output dir
global dir_output = "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/"

// Working directory
global dir_work = "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/code/HHandUI/"

// Specify variable and experiment
global myvar = "valuevf"
global myexp = "benchmark"

// Create empty data set where to store results for all experiments
// Each observation is one experiment
cd $dir_work
foreach v in "b_0" "b_bar" "mu" {
	gen `v' = .
	}

foreach v in "mean" "gini" "p10" "p25" "p50" "p75" "p90" {
	forval t = 1/3 {
		gen `v'_HH`t' = .
		}
	}

save output, replace
clear

// Iterate over experiments
forval e = 1/8 {
	local dir_aux = "$dir_work" + "experiment" + "`e'" + "/"

	cd `dir_aux'
	// Load experiment info
	import delim using UI.txt

	global b_0 = v1[1]
	global b_bar = v1[2]
	global mu = v1[3]
	clear

	// Load aggregate variables
	import delim using aggvars.txt
	global r = v1[1]
	global tau = v1[2]
	global KLratio = v1[3]
	global average_z = v1[4]
	clear


	// Read csv files
	forval f = 1/3 {
		preserve
		import delimited simulation_`f'.csv, varnames(1) clear
		gen HHtype = `f'
		save temp, replace
		restore
		append using temp
		}
	erase temp.dta

	// Label HHtype variable
	label var HHtype "Household type"
	label define lab_hhtype  1 "Single Men" 2 "Single Women" 3 "Married Household"
	lab val HHtype lab_hhtype

	// Compute VF moments
	forval t = 1/3 {
		su $myvar if HHtype == `t', d
		foreach v in "mean" "p10" "p25" "p50" "p75" "p90" {
			gen `v'_HH`t' = `r(`v')'
		}
	}

	// Compute gini VF
	// TBC

	// Keep just one observation and relevant variables
	keep if _n == 1
	keep mean* p*
	drop period

	// Add variables read in txt files
	foreach v in "b_0" "b_bar" "mu" "r" "tau" "KLratio" "average_z" {
		gen `v' = ${`v'}
		}

	// Append to output.dat
	local f = "$dir_work" + "output.dta"
	ap using `f'
	cd $dir_work
	save output, replace
	clear
	}

// Open output dataset
cd $dir_work
use output.dta

// Normalise all variables with as percentage change with respect to benchmark
forval t = 1/3 {
	foreach v in "mean" "p10" "p25" "p50" "p75" "p90" "gini" {
		local bm = `v'_HH`t'[6]
		replace `v'_HH`t' = ((`v'_HH`t' - `bm')/`bm')*100
		}
	}

// Label all variables
foreach v in "mean" "p10" "p25" "p50" "p75" "p90" "gini" {
	label var `v'_HH1 "Single Men"
	label var `v'_HH2 "Single Women"
	label var `v'_HH3 "Married"
	}

label var b_0 "Replacement ratio"
label var b_bar "Benefits cap"
label var mu "Duration"
label var tau "Income tax"
label var KLratio "K to L ratio"
label var average_z "Average z employed"
label var r "Interest rate"

// Plot all variables
set scheme plotplainblind

foreach v in "mean" "p10" "p25" "p50" "p75" "p90" {
	local file_aux = "$dir_output" + "$myexp" + "_" + "$myvar" + "_" + "`v'" + ".png"
	forval t = 1/3{
		qui su `v'_HH`t'
		local max`t' = `r(max)'
		local min`t' = `r(min)'
		}
		local ub = ceil(max(`max1',`max2',`max3'))
		local lb = floor(min(`min1',`min2',`min3'))

	if abs(`ub') + abs(`lb') > 2  {
		line `v'_HH1 `v'_HH2 `v'_HH3 b_0, ymtick(`lb'(1)`ub')
		}
	else {
		line `v'_HH1 `v'_HH2 `v'_HH3 b_0
		}
	gr export `file_aux', replace
}


// Erase output file
erase output.dta
