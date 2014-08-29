import os
import glob
import math
import time

import ROOT
from ROOT import *
from ROOT import gROOT, gStyle, gSystem, TLatex
import subprocess

from subprocess import Popen
from optparse import OptionParser

############################################
#            Job steering                  #
############################################
parser = OptionParser()
parser.add_option("-f","--inputFilePath"     , dest="inputFilePath"    , type="string", help="list of the file to  be processed -> path is enough")
parser.add_option("-w","--workdir"           , dest="workdir"          , type="string", help="Name of the directory for jobs")
parser.add_option("-c","--config"            , dest="config"           , type="string", default="TMVAApplyWeights_cfg.py",help="Ntuplizer config file")
parser.add_option("-e","--executable" , dest="executable"  , type="string", default="TMVAApplyWeights",help="Name of the executable. Default is: OptimizeSelections")
parser.add_option("-q","--queue"      , dest="queue"       , type="string", default="1nh",help="Name of the queue on lxbatch")
parser.add_option("-s","--submit"     , dest="submit"      , type="int"   , default=1,help="just create or create + submit")
parser.add_option("-a","--storeOneos" , dest="storeOneos"  , type="int"   , default=1,help="to understand if files have open from eos and read form eos or not")
parser.add_option("-n","--nfileperJob" ,dest="nfileperJob" , type="int"   , default=1,help="number of input files per each job")

(options, args) = parser.parse_args()
############################################

#-----------------------------------------                                                                                                                                                
#--- MAIN                                                                                                                                                                              
#-----------------------------------------                                                                                                                                              
                                                                                                                                                                        
path   = os.getcwd() ## path where create working directory
config = path+'/'+options.config  ## run the job where the template config is located                                                                                              
workingdir = path+'/'+options.workdir ## working directory

if os.path.isdir(workingdir): ## if exist, remove working dir
   os.system("rm -r "+workingdir);

os.system("mkdir -p "+workingdir); ## create working dir

lines = [];
njobs = 0 ;

if os.path.isfile("temp_list.txt"):
 os.system("rm -r temp_list.txt"); ## create a temp list of files

if options.storeOneos == 1 : ## if input files are stored  on eos use cmsLs command to read the input directory and make the file list
  os.system("cmsLs "+options.inputFilePath+" | grep root | awk '{print $5}' >> temp_list.txt");
else:
  os.system("ls "+options.inputFilePath+" | grep root >> temp_list.txt");

inputFileList = open('temp_list.txt', 'r') ## open the input list and read the lines
lines = inputFileList.readlines();

if len(lines)/options.nfileperJob-int(len(lines)/options.nfileperJob) > 0.5:
  njobs = len(lines)/options.nfileperJob+1;
else:
  njobs = len(lines)/options.nfileperJob;

residualFiles = len(lines)-options.nfileperJob*njobs;
if residualFiles > 0 : njobs = njobs +1;

iLine = 0; ## count the line of file list

for ijob in range(njobs): ## loop on the lines
 jobdir = workingdir+"/JOB_%d"%(ijob); ## for each job create a directory
 os.system("mkdir -p "+jobdir);
 os.system("cp "+options.config+" "+jobdir); ## copy the cfg in the jobdir
 
 inputListLine = ""; ## create the vector of string with the input file for each job
 for ifile in range(options.nfileperJob): ## for each job loop
  if iLine >= len(lines) : continue ;
  line = lines[iLine];
  line = line[:-1];
  if options.storeOneos == 1 :
        inputListLine = inputListLine+"root://eoscms.cern.ch//"+line+"\\\",\\\"" ;
  else:
        inputListLine = inputListLine+line+"\\\," ;
 
  iLine = iLine+1;

 command = "cat "+jobdir+"/"+options.config+" | sed -e s%INPUTFILELIST%\\\""+str(inputListLine)+"\\\"%g > "+jobdir+"/tmp.txt"; ## make a tmp config file
 os.system(command);
 command = "mv "+jobdir+"/tmp.txt"+" "+jobdir+"/"+options.config; ## move in the new one
 os.system(command);

 #--- prepare the jobs scripts                    
 jobscript = open('%s/sub_%d.sh'%(jobdir,ijob),'w')
 jobscript.write('cd %s \n'%jobdir)
 jobscript.write('eval ` scramv1 runtime -sh ` \n')
 jobscript.write('cd - \n')
 jobscript.write('if ( \n')
 jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,ijob));
 jobscript.write('\t %s %s'%(options.executable,jobdir+"/"+options.config));
 jobscript.write(') then \n');
 jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,ijob));
 jobscript.write('else \n');
 jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,ijob));
 jobscript.write('fi \n');
 os.system('chmod a+x %s/sub_%d.sh'%(jobdir,ijob));

if options.submit :
    print "for job in range(njobs): ",njobs
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(workingdir,job)
        jobname = '%s/sub_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/sub_%d.log %s'%(options.queue,jobdir,job,jobname)
        os.system('bsub -q %s -o %s/sub_%d.log %s'%(options.queue,jobdir,job,jobname));
