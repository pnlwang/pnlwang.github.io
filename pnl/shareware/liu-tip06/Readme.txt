Jiangye,

Attached is a gzipped tar file ('seg.source.osu.tar.gz') that includes the source code (*.h, *.c, makefile and just type in "make" to generate executable seg). Then you can run the program by doing
 ./seg data.20.par

Please let me know if you encounter any problems. Note that it runs on a linux machine or a cygwin environment on a Windows machine.

Xiuwen 

------------------------------------------------

Jiangye,

Sorry for that. After I sent you the files, I realized that I did not the filter list file and I thought you have those files. Anyway, attached is a tar file (othe_files.tar) that includes feature_filter (which specifies the filters to be used), data_20_feature_filter_mark.dat (histogram bin parameters, which can be generated automatically), and kernels.tar.gz. "kernels.tar.gz" includes all the filters I used. To do this, first you need to create a directory and gunzip and untar the filters to the directory and then you need to change the filters in feature_filter to the new location where these filters can be found. For the example, you need to change the first filter "/cavis/group/liux/percept2/kernels/intensity_1x1.ker" to the new location where intensity_1x1.ker can be found.

By the way, you will find another tar file (refineBoundary.tar), which implements the boundary refinement.

Let me know if you have any questions.


Xiuwen 