/* Nathaniel Simone
   STA5197 Final Project */
  
/* Input the data to SAS from a .txt file. Ignore the first line. */
data obese;
	infile '/home/u59238848/obesity.txt' firstobs=2;
	input ID Sex Baseline_Age Current_Age Occasion Obesity;
run;

/* Obtain estimates of obesity rate by sex*/
proc sort data=obese out=obese_s;
	by obesity occasion;
run;

proc freq data=obese_s;
	tables sex / nocum;
	by obesity;
	
	/* Create bar charts for the distribution of sex by obesity status */
	title 'Plots of Distribution of Sex by Obesity';
	tables sex / plots=Freqplot(scale=Percent);
run;


/* Create bar charts for the distribution of sex by obesity status at each of the 3 occasions */
proc freq data=obese_s;
	tables sex / nocum;
	by obesity occasion;
	
	/* Create bar charts for the distribution of sex by obesity status */
	title 'Plots of Distribution of Sex by Obesity per Occasion';
	tables sex / plots=Freqplot(scale=Percent);
run;

/* Fit a Saturated Model with only sex and occastion in the model */
proc genmod data=obese;
	class ID Sex (ref='0') Occasion (ref='1');
	model obesity = Sex | Occasion / link=logit;
	repeated subject=ID / type=ind covb corrw;
	output out=resids resraw=residuals predicted=pred;
run;

/* Fit spag plot of residuals vs. time */
proc surveyselect data=resids out=resids_samp
	sampsize=100 method=srs seed=1234;
	cluster ID;
run;

data resids_samp;
	set resids_samp;
	tij = current_age-baseline_age;
run;

proc sgplot data=resids_samp;
	series x=tij y=residuals / group = ID;
run;

proc sgplot data=resids_samp;
	series x=baseline_age y=residuals / group = ID;
run;

proc sgplot data=resids_samp;
	scatter x=baseline_age y=residuals;
run;


/* Plot a variogram of the residuals to determine the proper working correlation matrix structure */

proc variogram data=resids outpair=pairs;
	coordinates x=current_age y=id;
	compute robust novariogram;
	var residuals;
run;

data pairs2;
	set pairs;
	if y1 = y2;
	v = ((v1-v2)**2)/2;
run;

data pairs3;
	set pairs;
	if y1<y2;
	v = ((v1-v2)**2)/2;
run;

proc means data=pairs3 mean;
	var v;
run;
*0.1682172;

ods graphics / ANTIALIASMAX=12200 LOESSMAXOBS=15000;
proc sgplot data=pairs2;
	loess y=v x=distance;
	refline 0.1682172 / axis=y;
run;

/* Create a model containing estimate based on the occasion number, sex, and baseline age. */
proc genmod data=obese;
	class ID Sex (ref='0') Occasion (ref='1');
	model obesity = Occasion Sex*Occasion Baseline_age / link=logit;
	repeated subject=ID / type=cs covb corrw;
run;

proc genmod data=obese;
	class ID Sex (ref='0') Occasion (ref='1');
	model obesity = Occasion Baseline_age / link=logit;
	repeated subject=ID / type=cs covb corrw;
run;

/* Add a random intercept into a model to account for within subject variablility */
proc nlmixed data=obese method=gauss qpoints=10;
	dummy1=0; 
	dummy2=0; 
	dummy3=0; 

	if occasion = 1 then dummy1 = 1;
	if occasion = 2 then dummy2 = 1;
	if occasion = 3 then dummy3 = 1;

	parms b01=0 b02=0 b03=0 b11=0 b12=0 b13=0 b2=0 sb2=1;
	eta = b01*dummy1+b02*dummy2+b03*dummy3+b11*sex*dummy1+b12*sex*dummy2+b13*sex*dummy3+b2*baseline_age+bi;
	m=exp(eta)/(1+exp(eta));
	random bi ~ normal(0,sb2) subject=ID;
	model obesity ~ binary(m);
run;

/* Use Monte Carlo Method to generate data from nlmixed model */
/* Use this model to estimate odds ratio */
data generate;
	set obese;
	call streaminit(1234);
		U_ij = rand("normal", 0, sqrt(12.4242));
		B_01 = -4.1446;
		B_02 = -3.9200;
		B_03 = -3.7166;
		B_11 = 0.2492;
		B_12 = 0.4295;
		B_13 = 0.3036;
		B2 = 0.07121;
		
		if (sex=1 and occasion=2) then odds_1 = exp(B_02 + B_12 + B2*8 + U_ij);
		if (sex=0 and occasion=2)then odds_0 = exp(B_02 + B2*8 + U_ij);
run;

proc means data=generate mean;
	var odds_1 odds_0;
run;
