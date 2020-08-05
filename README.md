# rac-intensity
Title:
------
Spatial modeling for the probability of accidental mark locations on a shoe sole

Researchers: 
------------
Naomi Kaplan-Damary, Department of Statistics, The Hebrew University of Jerusalem, 91905 Jerusalem, Israel. naomi.kaplan@mail.huji.ac.il.

Micha Mandel, Department of Statistics, The Hebrew University of Jerusalem, 91905 Jerusalem, Israel. micha.mandel@mail.huji.ac.il.

Yoram Yekutieli, Hadassah Academic College Dept. of Computer Science 37 Hanevi'im St. P.O.Box 1114, 9101001 Jerusalem, Israel. yoramye@hac.ac.il.

Sarena Wiesner, Israel National Police Division of Identification and Forensic Science (DIFS) 1 Bar-Lev Road 91906 Jerusalem, Israel. sarenawiz@gmail.com.

Yaron Shor, Israel National Police Division of Identification and Forensic Science (DIFS) 1 Bar-Lev Road 91906 Jerusalem, Israel. yaronshor@gmail.com.

Description:
------------
The location of randomly acquired characteristics (RACs) is modeled here as a point process over the shoe sole.
A database of RACs collected by the Israeli Police Division of Identification and Forensic Science, which includes 13,000 RACs from 386 lab shoeprints (Wiesner et al., 2019) is used  to estimate its intensity function. 
The analysis is somewhat complicated as the shoes are differentiated by shape, level of wear and tear and contact surface. 
We present methods that take into account these challenges, either by using natural cubic splines on high resolution data, or by using a piecewise-constant model on larger regions defined by experts' knowledge.

Files
----- 
1.Programs:

	A. pixel analysis.R
		 Estimation using maximum resolution (pixel)
		 The code  produces Figures 2,3 and the global test of all coefficients equal to zero using the random effects model conducted in 4.5.1
		 In addition,the code produces Figures 1-4 in the Web appendix 

	B. Organizing data to include subsets (larger areas).R
		 In this code we adjust the data files to include 14 subsets and 36 subsets

	C. larger areas analysis.R
		 This is the code for estimation using larger areas
		 The code here provides Figures 5, Table 1, it is testing the hypothesis that the lambda parameters are equal for all j using the random effects model conducted in Section 5.2.
		 In addition it produces the histograms of contact surface in 14 areas in Section 6 of the web appendix.
		 In addition, using contacts_data.txt we create the list cont_use  and allcount matrix which are used to determine the contour of the shoe - this is only used to improve the appearance of figure 5. 
		 See the "Organizing data to include subsets (larger areas).R" for details.

	D. functions4simulation sub-sampling techniques.R
		 Contains the functions to be used in "sim_logistic_reg sub-sampling techniques" to create a simulation that compares the random effects and CML estimators
  	 using different types of within-cluster case-control sub-sampling and sub-sampling across the whole data frame.

	E. sim_logistic_reg sub-sampling techniques.R
		 Creates a simulation that compares the random effects and CML estimators.
  	 using different types of within-cluster case-control sub-sampling and sub-sampling across the whole data frame
		 using the functions in "functions4simulation sub-sampling techniques".
		 The code  produces Figures 6.


2. Data sets:


	A. locations_data.CSV: A data set of RAC locations. 
     The first three columns are used. 
     The first column indicates the shoe number.
     The second indicates the x axis of the RAC location .
     The third indicates the Y axis of the RAC location.
 
	B. contacts_data.txt: A data set of the contact surface 
     This is a pixel data where 1 indicates there is a contact surface and 0 otherwise.
     There are 307 columns in each shoe and 395 is the number of rows.
     The number of shoes is 387 but 386 is the number of shoes with RACs - shoe 127 has no RACS.

Last updated: 8/5/2020


These analyses are created as part of the manuscript "Spatial modeling for the probability of accidental mark locations on a shoe sole": https://arxiv.org/abs/1912.08272.  
