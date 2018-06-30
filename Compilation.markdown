# Compilation instructions
1. Navigate to the folder where the `.f90` files are located. In my machine: `cd ~/Dropbox/Choi_Valladares_2015/QEresubmission/code/HHandUI`
2. Create an executable file with gfortran:
  - *Normal*: `gfortran -g -fcheck=all -Wall Globals.f90 Utils.f90 BellmanEq.f90 ExpectedValues.f90 Main.f90 -o main.out`
  - Optimised: `gfortran Globals.f90 Utils.f90 BellmanEq.f90 ExpectedValues.f90 Main.f90 -O3 -o main.out`
3. Execute: `./main`

## [gfortran flags](http://faculty.washington.edu/rjl/classes/am583s2013/notes/gfortran_flags.html)
- `-g`Generates extra debugging information usable by GDB.
- `-Wall`: Short for “warn about all,” this flag tells gfortran to generate warnings about many common sources of bugs, such as having a subroutine or function with the same name as a built-in one, or passing the same variable as an intent(in) and an intent(out) argument of the same subroutine. In spite of its name, this does not turn all possible -W options on.
