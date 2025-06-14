**********************************************************
*
* Program Name   : V_PAIRBINCI.SAS
* Type           : macro validation
* Description    : Validation of PAIRBINCI macro
*                  Including replication of published examples  
* 
* Author         : Pete Laud 
* 
* Date Created   : 2025-05-14
* Program Status : CREATED
*
**********************************************************;


*** Data from Table II of Fagerland et al, 2014;
*** Tango method for contrast=RD;
*** (-0.517, -0.026) in Fagerland Table V;
%pairbinci(1,1,7,12,level=0.95, bcf=FALSE, skew=FALSE, contrast = RD, dps = 3);

*** Tang method for contrast=RR;
*** (0.065, 0.907) in Fagerland Table VI;
%pairbinci(1,1,7,12,level=0.95, bcf=FALSE, skew=FALSE, contrast = RR, dps = 3);

*** SCAS method for contrast=RD (Laud 2025 Table 5.1);
*** (-0.528, -0.018);
%pairbinci(1,1,7,12,level=0.95, bcf=TRUE, skew=TRUE, contrast = RD, dps = 3);
*** SCAS method for contrast=RR (Laud 2025 Table 5.2);
*** (0.043, 0.928);
%pairbinci(1,1,7,12,level=0.95, bcf=TRUE, skew=TRUE, contrast = RR, dps = 3);

*** Check handling of boundary cases for RR;
*** p2 = 0 -> RR = Inf;
%pairbinci(0,7,0,12,level=0.95, bcf=FALSE, skew=FALSE, contrast = RR);
%pairbinci(0,7,0,12,level=0.95, bcf=TRUE, skew=TRUE, contrast = RR);
*** p1 = 0 -> RR = 0;
%pairbinci(0,0,7,12,level=0.95, bcf=FALSE, skew=FALSE, contrast = RR);
%pairbinci(0,0,7,12,level=0.95, bcf=TRUE, skew=TRUE, contrast = RR);
*** p1 = p2 = 0;
%pairbinci(0,0,0,12,level=0.95, bcf=FALSE, skew=FALSE, contrast = RR);
%pairbinci(0,0,0,12,level=0.95, bcf=TRUE, skew=TRUE, contrast = RR);

data Students;
  input PassPre PassPost Count;
  datalines;
  0 0 12 
  0 1 7
  1 0  1
  1 1 1
  ;
run;
   
proc freq data=Students;
  tables PassPre*PassPost / nopercent norow nocol agree;
  weight Count;
*  ods select CrossTabFreqs McNemarsTest;
run;
  

/*
*** Development code checking the tan transformation used in optimisation subroutine;
data tantest;
 pi = constant('pi'); *3.14159265358979323846264338327950288;
 fx = tan((-0.3609535 + 1)*pi/4);
 do i = -100 to 100;
    x = pi * (i/100 + 1) / 4;
    x_adj = pi * (max(1E-4 - 1, min(i/100, 1 - 1E-4)) + 1) / 4;
	fi = round(tan(x), 1E-10);
	fi_adj = round(tan(x_adj), 1E-10);
	output;
 end;
run;
*/



