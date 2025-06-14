********************************************************************************;
*
* Program Name   : PAIRBINCI.SAS
* Author         : Pete Laud, SSU, University of Sheffield
* Date Created   : 07 Feb 2025
* Version date   : 22 May 2025
* Repository     : Download latest version from https://github.com/PeteLaud/ratesci-sas
*
* Type           : Macro
* Description    : Computes two-sided confidence intervals (CI) for the difference
*                  between two paired binomial proportions 
*
* Macro            a,b,c,d = input values from 2x2 contingency table, where:
* variables:        
*  					 a is the number of pairs with the event (e.g. success) under both
*  					   conditions (e.g. treated/untreated, or case/control)
*  					 b is the count of the number with the event on condition 1 only (= x12) 
* 					 c is the count of the number with the event on condition 2 only (= x21)
*  					 d is the number of pairs with no event under both conditions \cr
*                  
*                  CONTRAST = indicator selecting contrast to be estimated:
*                             RD = paired risk difference
*                             RR = paired risk ratio / relative risk
*                  LEVEL = (2-sided) confidence level required, e.g. 0.95 for 95% CI
*                          NB this corresponds to a NI test at the (1-LEVEL)/2 
*                          significance level*;
*                  SKEW   = indicator for inclusion of skewness correction for confidence interval & tests
*                           FALSE - Tango/Tang asymptotic score method
*                           TRUE (default) - Skewness-corrected asymptotic score method (SCAS)
*                  BCF   = indicator for inclusion of 'N-1' variance bias correction 
*                           FALSE - Tango/Tang asymptotic score method
*                           TRUE (default) - Skewness-corrected asymptotic score method (SCAS)
*                  THETA0 = specified non-inferiority margin
*                          NB THETA=0 for CONTRAST=RD or THETA=1 for CONTRAST=RR
*                             corresponds to a superiority test
*                  DPS   = number of decimal places precision to be displayed
********************************************************************************;

*** Format for readable infinite estimates for RR;
proc format;
 value infty .I = "    Inf"
             other = [best12.];
run;

*** subroutine to evaluate score at a given value of x;
%macro fx(x); 

  %if "&contrast." = "RD" %then %do;
    theta = &x.;
    AA = 2 * n;
    BB = -1*b - c + (2 * n - b + c) * theta;
    CC = -1*c * theta * (1 - theta);
    p21 = (SQRT(BB**2 - 4 * AA * CC) - BB)/(2 * AA);

	denterm = n * (2 * p21 + theta * (1 - theta));
	if denterm = 0 then fxi = 1E8 * sign(b - c - n * theta);
    else fxi = (b - c - n * theta) / SQRT(lambda * denterm);

	* Skewness correction calculations;
    p12 = p21 + theta;
    if a = 0 then p11 = 0;
    else p11 = a / (a + d) * (1 - p12 - p21);
    p22 = 1 - p11 - p12 - p21;
    p2d = min(1, max(0, p21 + p11));
    p1d = p2d + theta;
	k = -1;
  %end;

  %if "&contrast." = "RR" %then %do;
    theta = &x.;
    AA = n * (1 + theta);
    BB = (a + c) * theta**2 - (a + b + 2*c);
    CC = c * (1 - theta) * (a + b + c) / n;
    p21 = divide((SQRT(max(0, BB**2 - 4 * AA * CC)) - BB), (2 * AA));

	denterm = max(0, n * ((1+theta)*p21 + (a+b+c)*(theta - 1)/n));
	if denterm = 0 then fxi = 1E8 * sign(((a+b) - 1 * (a+c)));
    else fxi = ((a+b) - theta * (a+c)) / SQRT(lambda * denterm);

	p12 = divide((p21 + (theta - 1) * (1 - d/n)), theta);
	p11 = divide((1 - d/n - (1 + theta) * p21), theta);
	p22 = (1 - p11 - p12 - p21);
	p1d = p11 + p12;
	p2d = p11 + p21;
    k = -theta;
  %end;

  V = max(0, (p1d * (1 - p1d) + k**2 * p2d * (1 - p2d) +
    k * 2 * (p11 * p22 - p12 * p21)) / n) * lambda;
  mu3 = (p1d * (1 - p1d) * (1 - 2 * p1d) +
      (k**3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
      3 * k * (p11 * (1 - p1d)**2 + p21 * p1d**2 - p1d * p2d * (1 - p1d)) +
      3 * (k**2) * (p11 * (1 - p2d)**2 + p12 * p2d**2 - p1d * p2d * (1 - p2d))
      ) / (n**2);
  if abs(mu3) < 1E-10 or abs(V) < 1E-10 then scterm = 0; * Avoids issues with e.g. x = c(1, 3, 0, 6);
  else scterm = divide(mu3, (6 * V**(3 / 2)));
  As = scterm;
  Bs = 1;
  Cs = -(fxi + scterm);
  num = (-Bs + sqrt(max(0, Bs**2 - 4 * As * Cs)));
  if (skew = "FALSE" | scterm = 0) then score = fxi;
  else score = divide(num, (2 * As));
  fx = score - Z;
  dummy = 0;
%mend fx;

* Subroutine to find lower or upper confidence limit;
%macro Tango (a, b, c, d, X0, X1, Z);
  iteration = 1; * Initializing iteration counter;
  change = 1;
  n = a + b + c + d;
  X0 = &X0;
  X1 = &X1;
  Z = &Z; 
  if bcf = "FALSE" then lambda = 1;
  else lambda = n / (n - 1);

  do until (abs(change) < 10**-(&dps.+1));* or fx1 = 0);* or fx0 <= 0);
    %if "&contrast." = "RR" %then %do;
	  if x0 = -1 then x0t = tan(pi * (1E-12)/4); *%fx(0) fails due to infinite estimate; 
	  else if x0 = 1 then x0t = 1E12; *tan function fails at pi/2;
	  else x0t = tan(pi * (min(1 , max( - 1, x0)) + 1) / 4);
      %fx(x0t);
      fx0 = fx;
	  if x1 = -1 then x1t = tan(pi * (1E-12)/4);
	  if x1 = 1 then x1t = 1E12;
	  else x1t = tan(pi * (min(1 , max( - 1, x1)) + 1) / 4);
      %fx(x1t);
      fx1 = fx;
    %end;
    %if "&contrast." = "RD" %then %do;
      %fx(x0);
      fx0 = fx;
      %fx(x1);
      fx1 = fx;
    %end;

*  +----------------------------------------------------+;
	if iteration = 1 and fx0 < 0 then do;
      x2 = x0;
	  x1 = x0;
	end;
	else do;
      * bisection is slower than secant but works better for RR;
      x2 = max(-1, min(1, x1 + 0.5 * (x1 - x0) * sign(fx0)*sign(fx1) )); 
	end;
    change = x2 - x1;
	
 *+-------------------------------------------------------------------------+
* Set x0 and x1 for next iteration of loop, and increment iteration counter;
*+-------------------------------------------------------------------------+;
    x0 = x1;
    x1 = x2;
*	output;
    iteration = iteration + 1;
  end;
%mend Tango;

%macro pairbinci(a, b, c, d, 
                 CONTRAST = RD, 
                 LEVEL = 0.95, 
                 BCF = FALSE, 
                 SKEW = FALSE, 
				 THETA0 = .,
                 DPS = 6);

  data tangodata(keep = a b c d LEVEL CONTRAST UCL LCL bcf skew 
                        chi2 pval_2sided theta0 Z_STAT pval_l pval_r);
    retain a b c d CONTRAST LEVEL LCL UCL;
	format ucl infty.;
	pi = constant('pi'); 
	length SKEW $8. BCF $8.;
    a = &a;
    b = &b;
    c = &c;
    d = &d;
	CONTRAST = "&contrast.";
    LEVEL = &level.;
    skew = "&skew.";
    bcf = "&bcf.";
    zval = probit(1 - (1 - &level)/2);

    %Tango(a, b, c, d, -1, 1, zval); 
    if contrast = "RD" then LCL = round(x2, 10**-&dps.); 
	else if "&contrast." = "RR" then do;
	  if x2 = -1 then LCL = 0;
      else LCL = round(tan((x2 + 1) * pi / 4), 10**-&dps.);
	end;
	
    %Tango(a, b, c, d, -1, 1, -zval);
    if contrast = "RD" then UCL = round(x2, 10**-&dps.);
	else if "&contrast." = "RR" then do;
	  if x2 = 1 then UCL = .I;
      else UCL = round(tan((x2 + 1) * pi / 4), 10**-&dps.);
	end;


	* Test for association;
	CHI2 = .;
	PVAL_2SIDED = .;
	THETA0 = &theta0.;
    %if "&contrast." = "RD" %then %do;
	  if &theta0. = . then theta0 = 0;
	  %fx(x = 0);
      CHI2 = score**2;
      PVAL_2SIDED = round(1 - PROBCHI(CHI2, 1), 10**-&dps.);
	%end;
    %if "&contrast." = "RR" %then %do;
	  if &theta0. = . then theta0 = 1;
	  %fx(x = 1);
      CHI2 = score**2;
      PVAL_2SIDED = round(1 - PROBCHI(CHI2, 1), 10**-&dps.);
	%end;

     *** Optional noninferiority tests;
      %fx(x = theta0);
      Z_STAT = score;
      PVAL_L = round(PROBNORM(Z_STAT), 10**-&dps.);
      PVAL_R = round(1 - PROBNORM(Z_STAT), 10**-&dps.);
	
    output;
  run;

  proc print;
  run;

%mend pairbinci; 

