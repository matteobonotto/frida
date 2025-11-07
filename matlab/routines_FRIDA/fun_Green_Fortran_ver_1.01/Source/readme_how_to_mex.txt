HOW TO INSTALL INTEL FORTRAN COMPILER AND MEX A .F90 ROUTINE



0. run mex -v -setup FORTRAN on MATLAB command window and see what's the answer


1. Install the correct version of Intel Parallel Studio XE and Microsoft Visual Studio Community.
INPORTANT: - choose the correct combination of Parallel Studio and Vsual Studio
	   - I got troubles for distributions of VS other than Community


2. (maybe not necessary, to be checked)
- create the directory usr/local/bin
- copy/paste all the content of C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.4.245\windows\bin into usr/local/bin
- create the directory usr/local/compiler/lib
- copy/paste all the content of C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.4.245\windows\compiler\lib into usr/local/compiler/lib


3. (still maybe not necessary) from MATLAB command window run
- setenv('IFORT_COMPILER19','C:\usr\local')
- setenv('PATH', [getenv('PATH') 'C:\usr\local\bin']);


4. (to be done only in case of error "Error using mex error #5149" after step 6): go to C:\Program Files\MATLAB\R2019a\bin\win64\mexopts and edit intel_fortran_##_vs20##.xml. Change the line
COMPFLAGS="/nologo /fpp /Qprec /fixed /MD /fp:source /assume:bscc $INCLUDE  $COMPDEFINES"
to (remove the /fixed option)
COMPFLAGS="/nologo /fpp /Qprec /MD /fp:source /assume:bscc $INCLUDE  $COMPDEFINES"


5. tell MATLAB to choose the ifort compiler: from MATLAB command window run (arrange intel_fortran_##_vs20## depending your version of Intel Parallel Studio XE and Microsoft Visual Studio)
mex -setup:'C:\Program Files\MATLAB\R2019a\bin\win64\mexopts\intel_fortran_19_vs2017.xml' FORTRAN


6. mex your routine, assuming that such routine (e.g. fun_thicksolenoid.f90) is written in a suitable way. From MATLAB command window run (-v = verbose, it gives a lot of explanations)
mex -v -largeArrayDims -output fun_thicksolenoid_f subfields.f90 data_structure.f90 fun_thicksolenoid.f90

NOTE: if 6. gives the error "Error using mex error #5149 ..." go to .4 and follow all the steps 