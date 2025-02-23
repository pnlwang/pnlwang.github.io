******************************************************************************
****  Matlab code for generating Gammatone Features (GF) and	   	  ****
****  Gammatone Frequency Cepstral Coefficients (GFCC, 30 dimensions)	  ****
****                                                              	  ****
****  Written by Yang Shao, and adapted by Xiaojia Zhao, Nov. 2011	  ****
****  Ohio State University Perception and Neurodynamics Laboratory (PNL) ****
******************************************************************************

Important Notice: These programs are not for public distribution. Prior approval must be obtained
	from the PNL for any distribution not authorized by the PNL. 

References:
	Shao Y., Srinivasan S., and Wang D.L. (2007): "Incorporating auditory feature uncertainties 
		in robust speaker identification," Proceedings of ICASSP-07, pp. IV.277-280.
	Shao Y. (2007): Sequential organization in computational auditory scene analysis. 
		Ph.D. Dissertation, OSU Department of Computer Science and Engineering 
		(http://www.cse.ohio-state.edu/pnl/theses.html)


****** This package contains the following routines ******

gtfeatuers.m - This is the main entrance function, reading a list of wav files, calculating corresponding GF and GFCC features, and writing into HTK format.

gen_gammatone.m - Produce an array of impulse responses from a Gammatone filterbank (by default it has 64 channels).

fgammatone.m - Apply Gammatone filter impulse response on the input signal, and output GF features.

calDelta.m - Given GF or GFCC features, calculate delta features.

gtf2gtfcc.m - Convert GF features into GFCC features.

loadData.m - Load HTK format files.

writeHTK.m - Write HTK format files.

hz2erb.m - Convert frequency scale in Hz to ERB-rate scale.

erb2hz.m - Convert ERB-rate scale to frequency scale in Hz.

loudness.m - Compute loudness level in Phons on the basis of equal-loudness functions.


In addition:

f_af_bf_cf.mat - Contain parameter values of equal-loudness functions from BS3383,
"Normal equal-loudness level contours for pure tones under free-field listening conditions", table 1. 


An example:

There is a control file called "list.txt" in this folder. It has the paths of two ASCII format speech files sampled at 16 kHz. To generate GF and GFCC features, call the routine gtfeatures.m as follows:

>> gtfeatures('list.txt', 16000);

** It will generate HTK format feature files in the same folder as speech files for you. Here the samping frequency is explicitly set to 16000, and the default value is 8000.

** The code can be modified to generate different dimensions of GFCC features.
