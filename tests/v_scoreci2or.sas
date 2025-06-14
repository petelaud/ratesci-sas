**********************************************************
*
* Program Name   : V_SCORECI2or.SAS
* Level / Study  : global reusable macro
* Type           : macro validation
* Description    : Validation of SCORECI macro for contrast = OR, unstratified
*                  Use in combination with R code (v_scoreci2or.R)
*					(Work in progress)  
* 
* Author         : Pete Laud 
* 
* Date Created   : 30 May 2025
* Program Status : IN DEVELOPMENT
*
*                  Discrepancies found:
*                  SAS PROC FREQ gives no CI for x1=n1 and x2=n2
*                  R: PropCIs::riskscoreci has a bug when x1=n1 and x2/n2 >~0.86
*                     PropCIs::riskscoreci omits 'N-1' correction
*                     PropCIs::riskscoreci gives warnings about acos function
*                     sasLM::RRmn1 omits 'N-1' correction
*                     sasLM::RRmn1 has discrepancy when x1=n1 and x2/n2 is small non-zero
*
**********************************************************;


*Validation of SCORECI macro, using a large random sample of parameters (alpha, n1, n0, e1, e0);
*Output results to a text/csv file for comparison against other sources:
*1: R PropCIs::orscoreci function;
*2: SAS/STAT v9.4 PROC FREQ (NB this will not calculate any interval for 0/n1-0/n0, or n1/n1-n0/n0);
*3: R ratesci::scoreci function for skewness-corrected interval;


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


%let nsamp=1000;

*Generate a random sample of numerators & denominators; 
*(alphas not used in PROC FREQ comparison);
data sample(drop=treatment response freq) longsample(drop=n1 n0 e1 e0 p1 p0);
 do i = 1 to &nsamp.;
  stratum = i;
  n1 = ceil(rand('uniform')*300); *select n1 between 1 and 300;
  n0 = ceil(rand('uniform')*300);
  e1 = ceil(rand('uniform')*(n1+1))-1; *select e1 between 0 and n1;
  e0 = ceil(rand('uniform')*(n0+1))-1;
  p1 = e1/n1;
  p0 = e0/n0;
  if mod(i,3)=0 then alpha=0.01; *select alpha level;
  else if mod(i,3)=1 then alpha=0.05;
  else alpha=0.1;
  output sample; *this wide version is the correct format for NON_INF macro to use;
  *output also in long format for proc freq to use;
   treatment="1test"; response=1; freq=e1; output longsample; *NB trts deliberately labelled in reverse because statxact does 2-1 not 1-2;
   treatment="1test"; response=2; freq=n1-e1; output longsample;
   treatment="2comp"; response=1; freq=e0; output longsample;
   treatment="2comp"; response=2; freq=n0-e0; output longsample;
 end;
run;

**Check random sample;
*proc gplot data=sample;
* plot p1*p2;
*run;
*quit;

* Run the below PROC FREQ comparison with different confidence levels;
%let alpha=0.05;
%let alpha=0.01;
*Run the macro on the sample of data points;
%SCORECI(DS = sample, 
         delta = 0.9, 
         contrast = OR, 
         LEVEL = 1-&alpha.,
         STRATIFY = FALSE, 
         SKEW = FALSE, 
         ORBIAS = FALSE, 
         BCF = TRUE /*FALSE*/, 
         OUTPUT = FALSE, 
         converge = 1E-12);
data sasvalor(keep=i e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
  set result;
  i=row;
run;

* Check against PROC FREQ MN intervals;
* (could also check cl=score(correct=no) vs %SCORECI with BCF=FALSE);
ods output oddsratiocls=orcls;
ods html close;
ods listing close;
proc freq data=longsample; *(where=(stratum=84));
  weight freq;
  by i;
  tables treatment*response / noprint or(cl=score/*(correct=no)*/) alpha=&alpha.;
run; 
ods html;
ods listing;

data check;
 merge sasvalor orcls (keep = i lowercl uppercl);
 by i;
 lcld = lowercl - l_bound;
 if uppercl = .I and u_bound = .I then ucld = 0;
 else  ucld = uppercl - u_bound;
run;

* SAS PROC FREQ gives no result when x1=n1 and x2=n2;
data missing;
 set check;
 where lowercl = .;
run;

* Use proc univariate to summarise differences between macro and PROC FREQ CIs;
proc univariate data=check;
 var lcld ucld;
run;
* max UCL diff is 162/163 vs 1/286;
* Discrepancy can be larger for extremely large UCL;
* because root-finding is conducted on a transformed scale;
proc print data=check;
  where i=591;
run;

*Run the macro again for export to R, this time using the randomly generated alphas;
* Omitting BCF & ORBIAS to check matching against other R packages sasLM and PropCIs;
%SCORECI(DS = sample,
         DELTA = 0.9, 
         contrast = OR,
         LEVEL = 1-alpha,
         STRATIFY = FALSE, 
         SKEW = FALSE, 
         BCF = TRUE, 
         ORBIAS = FALSE, 
         OUTPUT = FALSE);
*%SCORECI(DS=sample1,DELTA=0.9, contrast=OR,LEVEL=1-alpha,STRATIFY=FALSE, SKEW=TRUE, BCF=TRUE, ORBIAS=TRUE, OUTPUT=TRUE);
data sasval(keep=row e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
  set result;
* i=_n_;
run;
*Export to csv for validation against R output using program v_scoreci2.R; 
proc export data=sasval 
  outfile="&path.tests\sasval2or.csv"
  dbms=csv replace;
run;

*Run again, this time with skewness & bias correction;
%SCORECI(DS = sample,
         DELTA = 0.9, 
         contrast = OR,
         LEVEL = 1-alpha,
         STRATIFY = FALSE, 
         SKEW = TRUE, 
         BCF = TRUE, 
         ORBIAS = TRUE, 
         OUTPUT = FALSE);

*Export to csv for validation against R output using program v_scoreci2.R; 
data sasval(keep=e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set result;
* i=_n_;
run;
proc export data=sasval 
outfile="&path.tests\sasval2orskew.csv"
dbms=csv replace;
run;

ods html close;
ods listing;

*Run again, this time with skewness & bias correction;
%SCORECI(DS = sample(where=(i=950)),
         DELTA = 0.9, 
         contrast = OR,
         LEVEL = 1-alpha,
         STRATIFY = FALSE, 
         SKEW = TRUE, 
         BCF = TRUE, 
         ORBIAS = TRUE, 
         OUTPUT = TRUE);
%SCORECI(DS = sample(where=(i=4031)),
         DELTA = 0.9, 
         contrast = OR,
         LEVEL = 1-alpha,
         STRATIFY = FALSE, 
         SKEW = TRUE, 
         BCF = TRUE, 
         ORBIAS = FALSE, 
         OUTPUT = TRUE,
         converge = 1E-10);

