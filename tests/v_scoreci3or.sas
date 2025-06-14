**********************************************************
*
* Program Name   : V_SCORECI3or.SAS
* Level / Study  : global reusable macro
* Type           : macro validation
* Description    : Validation of SCORECI macro for contrast = OR with STRATIFY=TRUE
*                  Use in combination with R code (v_scoreci3or.R)
*					(Work in progress)  
* 
* Author         : Pete Laud 
* 
* Date Created   : 08 Jun 2025
* Program Status : IN DEVELOPMENT
*
*                  Discrepancies found:
*
**********************************************************;

*** Validation of SCORECI macro, using a random sample of parameters (alpha, n1, n0, e1, e0);
*** Output results to a text/csv file for comparison against R 'scoreci' function.
***  NOTE SAS/STAT PROC FREQ currently (as of SAS/STAT 15.2 2023) 
***  does *not* provide the stratified M-N method;

*** The validation code below requires the SCORECI macro code to have been run;
*** Some environments (e.g. MS Visual Studio Code) make it difficult to get the path
*** for the current program, so the user may need to run the macro manually.
*** In interactive SAS environment, the following code should submit the macro code
*** from the directory above the location of this test program.;
%macro grabpath; 
  %qsubstr(%sysget(SAS_EXECFILEPATH),
    1,
    %length(%sysget(SAS_EXECFILEPATH)) - %length(%sysget(SAS_EXECFILENAME)) - 6
  )
%mend grabpath;
%let path = %grabpath;
%let filename = "&path.scoreci.sas";
%include &filename.;
* Otherwise, point SAS to the location of the SCORECI macro code 
* (change folder location as appropriate);
* filename prog "C:\Documents\ratesci-sas";
* %inc prog(scoreci);

*** Select number of samples, number of strata, and stratum size;
%let nsamp = 1000;
%let nstrat = 4;
%let maxn = 11;

*Generate a random sample of numerators & denominators;
*data samples longsample; 
data samples(drop=treatment response freq) longsample(drop=n1 n0 e1 e0 p1 p0);
  do sample = 1 to &nsamp.;
  do stratum = 1 to &nstrat.;
   n1 = 9+ceil(rand('uniform')*&maxn.); *select n1 between 10 and &maxn;
   n0 = 9+ceil(rand('uniform')*&maxn.);
   e1 = 2+ceil(rand('uniform')*(n1-1))-1; *select e1 between 2 and n1;
   e0 = 2+ceil(rand('uniform')*(n0-1))-1;
   p1 = e1/n1;
   p0 = e0/n0;
   if mod(sample,3)=0 then alpha=0.01; *select alpha level;
   else if mod(sample,3)=1 then alpha=0.05;
   else alpha=0.1;
   output samples; *this wide version is the correct format for SCORECI macro to use;
   treatment="1test"; response=1; freq=e1; output longsample; *NB trts deliberately labelled in reverse because statxact does 2-1 not 1-2;
   treatment="1test"; response=2; freq=n1-e1; output longsample;
   treatment="2comp"; response=1; freq=e0; output longsample;
   treatment="2comp"; response=2; freq=n0-e0; output longsample;
  end;
 end;
run;

proc freq data=longsample;
  weight freq;
  tables stratum*treatment*response / cmh agree;
run;
%SCORECI(DS = samples, 
         contrast = OR, 
         DELTA = 0.9, 
         LEVEL = 0.95, 
         STRATIFY = TRUE, 
         weight = 3, 
         SKEW = FALSE, 
         ORBIAS = FALSE, 
         DISTRIB = bin, 
         OUTPUT = FALSE);


%macro runstrat(skew=FALSE, orbias=FALSE, distrib=bin, contrast=OR);

data _null_;
  set samples;
  call symputx("nsamples", sample); 
run;

%do sample = 1 %to &nsamples.;

 data onesample;
  set samples;
  if sample = &sample.;
  call symput('alpha', alpha);
 run;

 %SCORECI(DS=onesample, 
          contrast=&contrast., 
          DELTA=0.9,  
          LEVEL=1-&alpha., 
          STRATIFY=TRUE, 
          weight=3, 
          SKEW=&skew., 
          ORBIAS=&orbias., 
          DISTRIB=&distrib., 
          OUTPUT=FALSE);
 
 data oneresult;
  sample = &sample.;
  set result;
 run; 

 proc append base=results data=oneresult force;
 run;
 quit;

 proc datasets lib=work nolist;
  delete homtests weighting result oneresult onesample;
 run;
 quit;

%end;
%mend;

%runstrat(contrast = OR, skew = FALSE, orbias = FALSE);

data output;
 merge samples results;
 by sample;
run;

*Export to csv for validation against R output using program v_scoreci3.R; 
data sasval(keep=sample stratum e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set output;
run;

*Export to csv for validation against R output using program v_scoreci3.R; 
proc export data=sasval 
 outfile = "&path.tests\sasval3or.csv"
 dbms = csv replace;
run;

 proc datasets lib=work nolist;
  delete output results sasval;
 run;
 quit;

* (Repeat with SKEW=TRUE);
%runstrat(contrast = OR, skew = TRUE, orbias = TRUE);

data output;
 merge samples results;
 by sample;
run;

*Export to csv for validation against R output using program v_scoreci3.R; 
data sasval(keep=sample stratum e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set output;
run;

*Export to csv for validation against R output using program v_scoreci.R; 
proc export data=sasval 
 outfile = "&path.tests\sasval3orskew.csv"
 dbms = csv replace;
run;

 proc datasets lib=work nolist;
  delete output results sasval;
 run;
 quit;

