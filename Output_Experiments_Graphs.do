/*
Choi & Valladares-Esteban (2018)

This codes reads Monte-Carlo simulation data from the model and generates results
*/

cls
clear all
set more off, perm
set trace off

// Output dir
global dir_output = "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/Graphs/"

// Working directory
global dir_work = "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/code/HHandUI/"

// Specify variables and experiment
global myvars = "valuevf consumption"
global qs = 4 // Number of quantiles for conditional wealth and income
global myexp = "SingleMen"
global nexp = 32 // Number of experiments
global nbm = 8 // Number of benchmark experiment

// Create empty data set where to store results for all experiments
// Each observation is one experiment
cd $dir_work
foreach v in "b_0" "b_bar" "mu" {
	gen `v' = .
	}

foreach v of global myvars {
	foreach d in "mean" "p50"  {
		gen `v'_`d'_HH = .
	}
	foreach cvar in "wealth" "income" {
		forval myq = 1/$qs {
			gen `v'_`cvar'q`myq'_HH = .
		}
	}
}

save output_$myexp, replace
clear

// Iterate over experiments
forval e = 1/$nexp {
	local dir_aux = "$dir_work" + "$myexp" + "`e'" + "/"

	cd `dir_aux'
	// Load experiment info
	import delim using UI.txt

	global b_0 = v1[1]
	global b_bar = v1[2]
	global mu = v1[3]
	clear

	// Load aggregate variables
	import delim using aggvars.txt
	global r = ((v1[1]+1)^12)-1
	global tau = v1[2]
	global KLratio = v1[3]
	global average_z = v1[4]
	clear

	// Load model otuput
	import delim using model.txt
	local linia = 1
	forval td = 1/3 {
	forval tw = 1/3 {
		global TR`td'`tw' = v1[`linia']
		local linia = `linia' + 1
	}
	}
	foreach r in "ER" "UR" {
		global `r' = v1[`linia']
		local linia = `linia' + 1
	}
	clear

	// Read csv file
	import delimited simulation_1.csv, varnames(1) clear
	gen income = labourincome  + benefitsreceived

	// Compute myvars moments
	foreach v of global myvars {
		su `v', d
		foreach d in "mean" "p50" {
			gen `v'_`d'_HH = `r(`d')'
		}
		foreach cvar in "wealth" "income" {
			xtile aux = `cvar', n(4)
			forval myq = 1/$qs {
				su `v' if aux == `myq'
				if (`r(N)') > 0 {
					gen `v'_`cvar'q`myq'_HH = `r(mean)'
				}
				else {
				  gen `v'_`cvar'q`myq'_HH = .
				}

			}
			drop aux
		}
	}

	// Compute gini
	fastgini wealth
	gen gwealth_HH = `r(gini)'

	// Keep just one observation and relevant variables
	keep if _n == 1
	keep *_HH

	// Add variables read in txt files
	foreach v in "b_0" "b_bar" "mu" "r" "tau" "KLratio" "average_z" {
		gen `v' = ${`v'}
		}

	forval td = 1/3 {
	forval tw = 1/3 {
		gen TR`td'`tw' = ${TR`td'`tw'}
	}
	}
	foreach r in "ER" "UR" {
		gen `r' = ${`r'}
	}

	// Append to output.dat
	local f = "$dir_work" + "output_" + "$myexp" + ".dta"
	ap using `f'
	cd $dir_work
	save output_$myexp, replace
	clear
	}


// Open output dataset
cd $dir_work
use output_$myexp.dta

// Normalise all variables as percentage change with respect to benchmark
sort b_0
foreach v of var *_HH* {
	local bm = `v'[$nbm]
	gen wrtb_`v' = ((`v' - `bm')/`bm')*100
	gen cev_`v' = exp((`v' - `bm')*(1-0.9950))-1
	}

foreach v of var tau* KLratio* TR* ER* UR* {
	local bm = `v'[$nbm]
	gen wrtb_`v' = ((`v' - `bm')/`bm')*100
	}

// Rates and probabilities as percentage
foreach v of var TR* ER* UR* {
	replace `v' = `v' * 100
}

label var b_0 "Replacement ratio"
label var b_bar "Benefits cap"
label var mu "Duration"
label var tau "Income tax"
label var wrtb_tau "Income tax"
label var KLratio "K to L ratio"
label var wrtb_KLratio "K to L ratio"
label var average_z "Average z employed"
label var r "Interest rate"


// Plot all variables
set scheme plotplainblind

// Variables graphed alone
foreach v of var *_HH* *tau* *KLratio* *TR* *ER* *UR* {
	local file_aux = "$dir_output" + "$myexp" + "_" + "`v'" + ".eps"
	twoway (scatter `v' b_0) ///
	(mspline `v' b_0, lcolor(black)),  legend(off)
	gr export `file_aux', replace
}


// Erase output file
// erase output_$myexp.dta
