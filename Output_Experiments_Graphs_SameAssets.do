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
global myvars = "valuevf"
global qs = 4 // Number of quantiles for conditional wealth and income
global myexp = "benchmark"
global nexp = 32 // Number of experiments
global nbm = 8 // Number of benchmark experiment

// Create empty data set where to store results for all experiments
// Each observation is one experiment
cd $dir_work
foreach v in "b_0" "b_bar" "mu" {
	gen `v' = .
	}

forval t = 1/4 {
foreach v of global myvars {
	foreach d in "mean" "p50"  {
		gen `v'_`d'_HH`t' = .
	}
	}
}

save output_SameAssets, replace
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

	// Read csv files
	forval f = 1/3 {
		preserve
		import delimited simulation_`f'.csv, varnames(1) clear
		gen HHtype = `f'
		gen income = labourincome  + benefitsreceived
		save temp, replace
		restore
		append using temp
		}

	// Married with same assets distribution as singles
	preserve
	count if HHtype == 3
	sample `r(N)' if HHtype <= 2, count // Make sure sample of singles not bigger than married
	xtile room = wealth if HHtype <= 2, n(10)
	levelsof room, local(rooms)
	foreach r of local rooms {
		su wealth if room == `r' & HHtype <= 2
		replace room = `r' if HHtype == 3 & wealth >= `r(min)' & wealth <= `r(max)'
	}
	drop if room == .

	// For rooms with less married than singles
	bysort room: egen N_singles = total(HHtype <= 2)
	bysort room: egen N_married = total(HHtype == 3)
	bysort room: gen double u = runiform()
	bysort room u: gen needed = N_singles - N_married if N_singles > N_married
	su needed
	while `r(N)' > 0 {
		bysort room: gen copiar = 1 if _n <= needed & needed != . & HHtype == 3
		expand 2 if copiar == 1, generate(nova)
		drop N_* u needed copiar nova
		bysort room: egen N_singles = total(HHtype <= 2)
		bysort room: egen N_married = total(HHtype == 3)
		bysort room: gen double u = runiform()
		sort room u
		bysort room: gen needed = N_singles - N_married if N_singles > N_married
		su needed
	}
	drop N_* u needed

	// For rooms with more married than singles
	bysort room: egen N_singles = total(HHtype <= 2)
	bysort room: egen N_married = total(HHtype == 3)
	bysort room: gen double u = runiform()
	bysort room u: gen needed = N_married - N_singles if N_married > N_singles
	su needed
	while `r(N)' > 0 {
		bysort room: gen borrar = 1 if _n <= needed & needed != . & HHtype == 3
		drop if borrar == 1
		drop N_* u needed borrar
		bysort room: egen N_singles = total(HHtype <= 2)
		bysort room: egen N_married = total(HHtype == 3)
		bysort room: gen double u = runiform()
		sort room u
		bysort room u: gen needed = N_married - N_singles if N_married > N_singles
		su needed
	}
	drop N_* u needed

	keep if HHtype == 3
	replace HHtype = 4
	save temp, replace
	restore
	append using temp

	erase temp.dta

	// Compute myvars moments
	forval t = 1/4 {
	foreach v of global myvars {
		su `v' if HHtype == `t', d
		foreach d in "mean" "p50" {
			gen `v'_`d'_HH`t' = `r(`d')'
		}
	}
	}

	// Keep just one observation and relevant variables
	keep if _n == 1
	keep *_HH*

	// Add variables read in txt files
	foreach v in "b_0" "b_bar" "mu" "r" "tau" "KLratio" "average_z" {
		gen `v' = ${`v'}
		}

	// Append to output.dat
	local f = "$dir_work" + "output_" + "SameAssets.dta"
	ap using `f'
	cd $dir_work
	save output_SameAssets, replace
	clear
	}


// Open output dataset
cd $dir_work
use output_SameAssets.dta

// Normalise all variables as percentage change with respect to benchmark
sort b_0
foreach v of var *_HH* tau* KLratio* {
	local bm = `v'[$nbm]
	gen wrtb_`v' = ((`v' - `bm')/`bm')*100
	gen cev_`v' = exp((`v' - `bm')*(1-0.9950))-1
	}

// Label all variables
foreach v of var *_HH1 {
	label var `v' "Single Men"
}
foreach v of var *_HH2 {
	label var `v' "Single Women"
}
foreach v of var *_HH3 {
	label var `v' "Married"
}
foreach v of var *_HH4 {
	label var `v' "Married, same assets"
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

foreach tp in "" "wrtb_" "cev_" {
foreach v of global myvars {
	foreach s in "mean" "p50" {
		local file_aux = "$dir_output" + "SameAssets" + "_" + "`tp'" + "`v'" + "_" + "`s'" + ".eps"
		forval t = 1/4 {
			qui su `tp'`v'_`s'_HH`t'
			local max`t' = `r(max)'
			local min`t' = `r(min)'
			}
			local ub = ceil(max(`max1',`max2',`max3',`max4'))
			local lb = floor(min(`min1',`min2',`min3',`max4'))

		if abs(`ub') + abs(`lb') > 2  {
			twoway (scatter `tp'`v'_`s'_HH1 `tp'`v'_`s'_HH2 `tp'`v'_`s'_HH3 `tp'`v'_`s'_HH4 b_0) ///
			(mspline `tp'`v'_`s'_HH1 b_0, lcolor(black)) ///
			(mspline `tp'`v'_`s'_HH2 b_0, lcolor(gray)) ///
			(mspline `tp'`v'_`s'_HH3 b_0, lcolor(ltblue)) ///
			(mspline `tp'`v'_`s'_HH4 b_0, lcolor(green)), ymtick(`lb'(1)`ub') legend(order(1 2 3 4))
			}
		else {
			twoway (scatter `tp'`v'_`s'_HH1 `tp'`v'_`s'_HH2 `tp'`v'_`s'_HH3 `tp'`v'_`s'_HH4 b_0) ///
			(mspline `tp'`v'_`s'_HH1 b_0, lcolor(black)) ///
			(mspline `tp'`v'_`s'_HH2 b_0, lcolor(gray)) ///
			(mspline `tp'`v'_`s'_HH3 b_0, lcolor(ltblue)) ///
			(mspline `tp'`v'_`s'_HH4 b_0, lcolor(green)),  legend(order(1 2 3 4))
			}
		gr export `file_aux', replace
	}
}
}




// Erase output file
// erase output_$myexp.dta
