
******************************************************************************
****  Java code for multipitch tracking in noisy and reverberant speech	  ****
****                                                              	  ****
****  Written by Zhaozhang Jin, 2010	                                  ****
****  Ohio State University Perception and Neurodynamics Laboratory (PNL) ****
******************************************************************************

Important Notice: These programs are not for public distribution. Prior approval must be obtained
	from the PNL for any distribution not authorized by the PNL. 

Reference:
	Jin Z. and Wang D.L. (2011): "HMM-based multipitch tracking for noisy and reverberant speech,"
		IEEE Transactions on Audio, Speech, and Language Processing, vol. 19, pp. 1091-1102.


Usage: java -jar RevMultiPitch.jar mixture outpitch1 outpitch2

Contents:

- JMultiPitch2.java		main class
- casa/Gammatone.java		class of Gammatone Filterbank
- casa/FrontEnd.java		class of Preprocessing (e.g., envelope extraction)
- casa/Correlogram.java		class of Correlogram
- PitchProbModel.java		class of Pitch Probability Modeling and HMM

* JAR is a packaging solution for JAVA. It can include both source and compiled files, 
making it directly runnable on any platforms. To extract, edit, or update the package, 
please refer to manual by simply typing "jar" in command line.
