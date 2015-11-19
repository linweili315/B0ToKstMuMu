##########################
# B0 To Kst Mu Mu Analysis 
##########################

(a) You need to setup some environment variables with the script
“InitAnalysis.sh”
which is under the directory “ plugins/Scripts”
For example you can call it with the calling sequence:

export B0=/Users/dinardo/Documents/MyDocuments/CMS/PhysicsAnalysis/B0Analysis/B0KstMuMu/
source $B0/plugins/Scripts/InitAnalysis.sh $B0 0

(b) You need to change the directory of reading efficiency, which is 
under the directory "src/Utils.cc/ReadRTEffPDF & ReadWTEffPDF"
For example:

TFile* file=TFile::Open("/afs/cern.ch/user/l/llinwei/work/B0KstMuMu/efficiency/effKEpdf_out_RT.root","READ");
just change this to your efficiency's directory.

(c) You can now compile the fitter program:
make ExtractYield

(d) Every program that I wrote as a synopsis if you call it without
any argument.
Anyway, to run the fit you can for example type:

./ExtractYield 6 singleCand_B0ToKstMuMu_MC_NTuple.root yesEffCorr 0
