

** This directory contains some simple plot style scripts for CEPC analysers or someone who needs a simple plot utility.

** The origin code/script is invented by Nik (BESIII), here I made some improvement.

** two c++ files/scripts are provided: cepcplotstyle.C and cepcplotstyle.h

** The simplest usage in root interactively:
   ---load the files first
   .L cepcplotstyle.C
   ---then execute your script
   .x script.C
   --- examples in this directory 
   .x DataMCPlots.C // Histogram comparison 
   .x DataFitPlot.C // Fit data/models comparison 
   .x Data2DPlot.C  // Scattering plot 

**  An example of advanced usage 

   3 code files:
     plots.C cepcplotstyle.C  cepcplotstyle.h  
    
   1 makefile
   
   3 root files with SAME stuctures 
     h0bq-91.2GeV.root  h0cq-91.2GeV.root  h0lq-91.2GeV.root

   * In the command line, you only need to type 
   ~$ make 
   ~$ ./plots
   
    It will plot all the variables in the root file and put them together for comparision, the figures will be put in the sub-dir ./figs


Any problem, please contact Gang LI   ---- li.gang@ihep.ac.cn or discussed in the CEPC QQ group: 392102847 ( encouraged ) 

