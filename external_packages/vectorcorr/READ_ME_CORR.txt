                   Description of operation of script Vector_Corr

Operation of Vector_Corr is the same as for any other MATLAB script.  

Start MATLAB as usual, and set the Current Directory to that in which the 30
scripts and functions (copied from the IAMG website) of Vector_Corr are
located.  The MATLAB-supplied functions and the MATLAB Toolbox also must
be accessible for execution.  

Transformations or operations on the data may be required (e.g., convert from
time to angles; convert azimuths from the range [-180, 180] to [0, 360]).  If so,
read the data file and perform the required operations with standard MATLAB
commands. 

Start the script by typing Vector_Corr on the command line in the Command
Window; then press Enter.  The following sequence of windows allows setting
controls for operating the script.  The script checks for errors and omissions in
filling out the window entries.  If an error is detected, a message will appear in
the Command Window and the GUI window will remain on the screen.     

1.  A window appears that 
A) allows entering a brief description of the job or data being analyzed.  This is
printed with calculations and on generated figures.
B) allows the user to specify if the data to be analyzed are to be read from an
External File, or if they are already in MATLAB as a Data Array.

2.  If an External File was specified above (1.B), then a standard window opens
that allows the user to specify the name of the text file.  See below for
information on the data file.

3a.  If an External File was specified in 1.B, then a window opens that requests
information about the file.  This contains entries as follows:
A) Number of header records in the file to skip (0, 1, 2, ...).
B) Number of columns in the data set; two or more of these will contain variables
(azimuths and/or linear) to be analyzed.
C) Specify if azimuthal variables are in degrees or radians.
Following are 12 sets of three entries.  These sets are labeled X1 - X12 and
consist of controls that specify up to 12 variables to be analyzed.  Entries D), E),
and F) are used for each variable.
D) Number of the column that contains the variable (linear or azimuthal). 
E) Alphabetic identifier of the variable (optional).  
F) Check box if this variable is azimuthal.
G) Enter the level of significance, alfa, to be used for the tests of hypothesis. 
Some tests are tabled only for alfa equal to 0.10, 0.05, 0.025, and 0.01 for small
sample sizes.  For larger sample sizes or other tests, any value of alfa may be
used.  If a value for alfa is entered that is not one of these four values, alfa
will be modified to the closest of the four.
H) If resampling methods are used, specify the number of permutation trials to use.
Must be in range (100, 10000), Default is 2500.
I) Are cross plots to be generated?
J) Specify if the calculations are to be written to an output text file.  All
calculations also appear in the Command Window, but this file allows permanent
storage of the results.

3b.  If a Data Array was specified in 1.B, then a window opens that requests
information about the array.  This contains entries as follows:
A) Name of the MATLAB data array that contains the azimuth data.  The default
name is AzData.  If Vector_Corr read an External File when previously executed,
that data will be left in array AzData.
B) Specify if azimuthal variables are in degrees or radians.
Following are 12 sets of three entries.  These sets are labeled X1 - X12 and
consist of controls that specify up to 12 variables to be analyzed.  Entries C), D),
and E) are used for each variable.
C) Number of the column that contains the variable (linear or azimuthal). 
D) Alphabetic identifier of the variable (optional).  
E) Check box if this variable is azimuthal.
F) Enter the level of significance, alfa, to be used for the tests of hypothesis. 
Some tests are tabled only for alfa equal to 0.10, 0.05, 0.025, and 0.01 for small
sample sizes.  For larger sample sizes or other tests, any value of alfa may be
used.  If a value for alfa is entered that is not one of these four values, alfa
will be modified to the closest of the four.
G)  If resampling methods are used, specify the number of permutation trials to use.
Must be in range (100, 10000), Default is 2500.
H) Are cross plots to be generated?
I) Specify if the calculations are to be written to an output text file.  All
calculations also appear in the Command Window, but this file allows permanent
storage of the results.

4.  If calculations are to be written to the optional output file (3a.J or 3b.I), a
window opens that requests the name of the file to be written.  The file should be
of type .txt.

The script then performs the calculations and generates the plots and output text
file if requested.

                             Error return codes

The script returns error codes to the Command Window if problems develop
while filling out the GUI or during execution.

0 = execution OK
1 = the requested data set could not be found or could not be opened
2 = an incorrect value was given for the number of columns (3a.B) in data set, 
  or character-data are in file, or incomplete record is in file;
  file could not be read   
3 = the external output file could not be opened
7 = no name specified for output file, although it had been requested
8 = job was cancelled by user while inputting data-set controls

                           Format of input data set         

The data set read by Vector_Corr consists of (1) header records and (2) an array
of numbers arranged in rows and columns, at least two columns of which 
contain variables of interest to correlate.  

The file begins with header records; any number are allowed, including 0. 
During operation of the GUI, the user specifies the number of records to be
skipped.   

Following the headers are the data records, each representing one observation of
sets of linear and/or azimuth variables. These are the rows in the data array. 
Each record contains two or more columns.  These columns contain the several
variables that are associated with each observation.  These rows and columns
thus make up the data array processed by MATLAB.  

With the exception of the header records, no character or string variables may be
in the file.  If such variables are in the data set, read the file separately into
MATLAB and pick up the resulting array for analysis.

The script may also process a data array that already is in MATLAB's memory. 
If so, the user simply specifies the name of that array, and indicates which
column is to be used for analysis.  If the Vector_Stats or Vector_Corr
script is used to read a text file, the data array that is input will be left in
MATLAB memory, under the name AzData, after execution.       
