import os
import sys
import numpy
import glob
import ROOT
from ROOT import TFile, TTree

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inputTrainingList" , dest="inputTrainingList"  , type="string", help="Path where list of training setup is located. Example: ./fileSingleVariables.txt")
parser.add_option("-w","--workdir"    , dest="workdir"    , type="string", default="mydir",help="Name of the directory for jobs")
parser.add_option(""  ,"--outputdir"  , dest="outputdir"  , type="string", default="",help="Name of the local afs output directory for jobs")
parser.add_option("-c","--config"     , dest="config"     , type="string", default="OptimizeSelections_cfg_lowPT_lowPU.py",help="Ntuplizer config file. Default is: OptimizeSelections_cfg_lowPT_lowPU.py")
parser.add_option("-m","--methodName" , dest="methodName" , type="string", default="BDT",help="MVA method to be used")
parser.add_option("-e","--executable" , dest="executable" , type="string", default="OptimizeSelections",help="Name of the executable. Default is: OptimizeSelections")
parser.add_option("-q","--queue"      , dest="queue"      , type="string", default="1nh",help="Name of the queue on lxbatch")
parser.add_option("-s","--submit"     , dest="submit"     , type="int"   , default=1,help="just create or create + submit")

(options,args)=parser.parse_args()

#-----------------------------------------                                                                                                                                                
#--- MAIN                                                                                                                                                                                
#-----------------------------------------                                                                                                                                                 
path = os.getcwd()
conf = path+'/'+options.config  ## run the job where the template config is located
workingdir = path+'/'+options.workdir
outputdir = '' 

if os.path.isdir(workingdir):
   os.system("rm -r "+workingdir);

os.system("mkdir -p "+workingdir);

if (options.outputdir !=''):        
   outputdir = workingdir+"/"+options.outputdir ;
   mkdir = 'mkdir -p '+outputdir ;
   os.system(mkdir)

lines = [];
listofVariables = [];
njobs = 0 ;

file = open('%s'%options.inputTrainingList, 'r')
lines = file.readlines();

for iLine in range(len(lines)):
 listofVariables = lines[iLine].split(';');
 if listofVariables[0][0] == "#" : continue ; 
 jobdir = workingdir+"/JOB_%d"%(njobs); 
 os.system("mkdir -p "+jobdir);
 os.system("cp "+options.config+" "+jobdir);
 command = "cat "+jobdir+"/"+options.config+" | sed -e s%INPUFILELIST%"+str(listofVariables[0])+"%g > "+jobdir+"/temp.txt";  
 os.system(command);
 command = "cat "+jobdir+"/temp.txt | sed -e s%LABELNAME%"+str(listofVariables[1])+"%g > "+jobdir+"/temp2.txt";  
 os.system(command);
 command = "cat "+jobdir+"/temp2.txt | sed -e s%METHOD%"+options.methodName+"%g > "+jobdir+"/temp.txt";  
 os.system(command);
 command = "cat "+jobdir+"/temp.txt | sed -e s%OUTPUTDIR%"+outputdir+"%g > "+jobdir+"/"+options.config;  
 os.system(command);
 os.system(" rm "+jobdir+"/*temp*");
 
 #--- prepare the jobs scripts                                                                                                                                                      
 jobscript = open('%s/sub_%d.sh'%(jobdir,njobs),'w')
 jobscript.write('cd %s \n'%jobdir)
 jobscript.write('eval ` scramv1 runtime -sh ` \n')
 jobscript.write('cd - \n')
 jobscript.write('if ( \n')
 jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,njobs));
 jobscript.write('\t %s %s'%(options.executable,jobdir+"/"+options.config));
 jobscript.write(') then \n');
 jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,njobs));
 jobscript.write('else \n');
 jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,njobs));
 jobscript.write('fi \n');
 os.system('chmod a+x %s/sub_%d.sh'%(jobdir,njobs));
 njobs = njobs +1
 

if options.submit :
    print "for job in range(njobs): ",njobs
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(workingdir,job)
        jobname = '%s/sub_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/sub_%d.log %s'%(options.queue,jobdir,job,jobname)
        os.system('bsub -q %s -o %s/sub_%d.log %s'%(options.queue,jobdir,job,jobname));

