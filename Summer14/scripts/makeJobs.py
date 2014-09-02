import os
import sys
import numpy
import glob
import ROOT
from ROOT import TFile, TTree

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inputdir"   , dest="inputdir"   , type="string", help="Path where ntuples are located. Example: /store/group/phys_jetmet/ntran/PUPPI/miniSamples/62x/qcd300-470_62x_PU40BX50")
parser.add_option("-w","--workdir"    , dest="workdir"    , type="string", default="mydir",help="Name of the directory for jobs")
parser.add_option("-o","--outputname" , dest="outputname" , type="string", default="outtree",help="Name of the output file. Default is: outtree")
parser.add_option(""  ,"--eosdir"     , dest="eosdir"     , type="string", default="",help="Name of the eos output directory for jobs")
parser.add_option("-c","--config"     , dest="config"     , type="string", default="minintuplizer_cfg.py",help="Ntuplizer config file. Default is: minintuplizer_cfg.py")
parser.add_option("",  "--puppiConfig", dest="puppiConfig", type="string", default="Puppi_cff.py",help="Puppi config file. Default is: Puppi_cff.py")
parser.add_option("-n","--njobs"      , dest="njobs"      , type="int"   , help="Number of jobs")
parser.add_option("-e","--executable" , dest="executable" , type="string", default="MiniNtuplizer",help="Name of the executable. Default is: MiniNtuplizer")
parser.add_option("-q","--queue"      , dest="queue"      , type="string", default="1nh",help="Name of the queue on lxbatch")
parser.add_option(""  ,"--checkJobs"  , dest="checkJobs"  , action="store_true", default=False,help="Checks job status")
parser.add_option(""  ,"--resubmit"   , dest="resubmit"   , action="store_true", default=False,help="Resubmit job ")
parser.add_option("-j","--jobmin"     , dest="jobmin"     , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option("-x","--jobmax"     , dest="jobmax"     , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option(""  ,"--dryRun"     , dest="dryRun"     , default=False, action="store_true",help="Do not submit jobs")
parser.add_option(""  ,"--eventsPerJob" , dest="eventsPerJob", default=3500, type="int",help="Number of events in each job .. don't look at the number of files")
parser.add_option(""  ,"--jobtype" , dest="jobtype", default=0, type="int",help="1 means specify number of jobs, 0 means number of jobs per event")


(options,args)=parser.parse_args()

eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'

def makeFilesList(indir,wdir):
    list = []
    command = ('%s find -f %s | grep root > %s/list.txt' % (eos,indir,wdir))
    #print command
    os.system(command)
    file = open('%s/list.txt'%wdir, 'r')
    list =[line.replace('/eos/cms/','root://eoscms.cern.ch//').replace('\n','') for line in file ]
    print 'Found %d files' %len(list)
    #print list             
    return list


def writeJobs(wdir, analysis, config, puppiconfig, indir, output, eosoutdir, njobs):
    #---------------------------------------------
    # --- prepare the list of files to be analyzed
    #---------------------------------------------
    listoffiles = []
    listoffiles = makeFilesList(indir,wdir) ## make the file list for the input directory on eos

    #---------------------------------------------
    # --- now split the jobs
    #---------------------------------------------
    if options.jobtype == 0:
      ## loop on the file list
      jobid = 0 ;
      for ifile in listoffiles:
       tfile    = TFile.Open(ifile,"READ");
       if not tfile : continue 
       ttree    = tfile.Get("Events");
       if not ttree : continue 
       nentries = ttree.GetEntries();
       
       if nentries/options.eventsPerJob - int(nentries/options.eventsPerJob) > 0.5 :  
         njobs    = int(nentries/options.eventsPerJob)+1;  
       else:
         njobs    = int(nentries/options.eventsPerJob);  

       residualEvents = nentries-options.eventsPerJob*njobs ;
       if residualEvents > 0 : njobs = njobs+1;
       for i in range(njobs):         
        jobdir = '%s/JOB_%d'%(wdir,jobid)
        os.system("mkdir -p "+jobdir)        
        f = open('%s/input_%d.txt'%(jobdir,jobid), 'w')
        f.write('%s \n'%ifile)
        os.system("cp "+config+" "+jobdir);
        configName = config.split("/");

        os.system("cat "+jobdir+"/"+configName[len(configName)-1]+" | sed -e s%MAXEVENTS%"+str(options.eventsPerJob*(i+1))+"%g > "+jobdir+"/temp.txt");
        if i == 0: 
         os.system("cat "+jobdir+"/temp.txt | sed -e s%MINEVENTS%"+str(0)+"%g > "+jobdir+"/temp2.txt");
        else:
         os.system("cat "+jobdir+"/temp.txt | sed -e s%MINEVENTS%"+str(options.eventsPerJob*i+1)+"%g > "+jobdir+"/temp2.txt");
        os.system("mv "+jobdir+"/temp2.txt "+jobdir+"/"+configName[len(configName)-1]) ;
        os.system("rm "+jobdir+"/temp.txt");

        os.system("cat "+jobdir+"/"+configName[len(configName)-1]+" | sed -e s%PUPPICONFIG%"+(puppiconfig.split('/')[-1])+"%g > "+jobdir+"/temp3.txt");
        os.system("mv "+jobdir+"/temp3.txt "+jobdir+"/"+configName[len(configName)-1]) ;
        
        #--- prepare the jobs scripts
        jobscript = open('%s/sub_%d.sh'%(jobdir,jobid),'w')
        jobscript.write('cd %s \n'%jobdir)
        #jobscript.write('export SCRAM_ARCH=slc5_amd64_gcc472 \n')
        jobscript.write('export SCRAM_ARCH=slc6_amd64_gcc472 \n')
        jobscript.write('eval ` scramv1 runtime -sh ` \n')
        jobscript.write('cd - \n')
        jobscript.write('cp %s ./ \n'%(puppiconfig)) 
        jobscript.write('if ( \n')
        jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,jobid))
        jobscript.write('\t %s %s %s/input_%d.txt %s_%d.root'%(analysis, jobdir+"/"+configName[len(configName)-1], jobdir, jobid, output, jobid))
        jobscript.write(') then \n')
        if (eosoutdir == ''):
            jobscript.write('\t cp ./%s_%d.root %s \n'%(output,jobid,jobdir))
        else:
            jobscript.write('\t cmsStage ./%s_%d.root %s/ \n'%(output,jobid,eosoutdir))
        jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,jobid))
        jobscript.write('else \n')
        jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,jobid))
        jobscript.write('fi \n')
        os.system('chmod a+x %s/sub_%d.sh'%(jobdir,jobid))
        jobid += 1; 

      options.njobs = jobid;
        
    elif options.jobtype == 1:
     #---------------------------------------------
     # --- now split the jobs
     #---------------------------------------------
     for job in range(njobs):
        jobdir = '%s/JOB_%d'%(wdir,job)
        os.system("mkdir -p "+jobdir)
        
        #--- prepare the list of files for each job
        f = open('%s/input_%d.txt'%(jobdir,job), 'w')
        sublist = [file for i,file in enumerate(listoffiles) if (i%njobs==job)]
        for fname in sublist:
            f.write('%s \n'%fname)

        #--- prepare the jobs scripts
        jobscript = open('%s/sub_%d.sh'%(jobdir,job),'w')
        jobscript.write('cd %s \n'%jobdir)
        #jobscript.write('export SCRAM_ARCH=slc5_amd64_gcc472 \n')
        jobscript.write('export SCRAM_ARCH=slc6_amd64_gcc472 \n')
        jobscript.write('eval ` scramv1 runtime -sh ` \n')
        jobscript.write('cd - \n')
        jobscript.write('cp %s ./ \n'%(puppiconfig))
        jobscript.write('if ( \n')
        jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,job))
        jobscript.write('\t %s %s %s/input_%d.txt %s_%d.root'%(analysis, config, jobdir, job, output, job))
        jobscript.write(') then \n')
        if (eosoutdir == ''):
            jobscript.write('\t cp ./%s_%d.root %s \n'%(output,job,jobdir))
        else:
            jobscript.write('\t cmsStage ./%s_%d.root %s/ \n'%(output,job,eosoutdir))
        jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,job))
        jobscript.write('else \n')
        jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,job))
        jobscript.write('fi \n')
        os.system('chmod a+x %s/sub_%d.sh'%(jobdir,job))

    else: 
      print" Error in defining jobtype .. options not known " ;
      return ;   

def submitJobs(wdir, njobs, queue):
    print "for job in range(njobs): ",njobs
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(wdir,job)
        jobname = '%s/sub_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/sub_%d.log -e %s/sub_%d.err %s'%(queue,jobdir,job,jobdir,job,jobname )
        os.system('bsub -q %s -o %s/sub_%d.log -e %s/sub_%d.err %s'%(queue,jobdir,job,jobdir,job,jobname ))
        


def checkJobs(wdir, output, queue, eosoutdir):
    jobs = glob.glob( '%s/JOB_*/sub*.sh'% (wdir) )
    print 'Total number of jobs: %d' %len(jobs)
    
    listdone = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.done' % (wdir,j,j))]
    print 'Total number of DONE jobs: %s ' % len(listdone)
    print '  %s' %listdone
    for j in listdone:
        f = '%s/JOB_%d/sub_%d.run'%(wdir,j,j)
        #print 'rm %s/JOB_%d/sub_%d.run'%(wdir,j,j)
        if (os.path.isfile(f)):
            os.system('rm %s'%f)
            
    listrun = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.run' % (wdir,j,j))]
    print 'Total number of RUNNING jobs: %d ' %len(listrun)
    print '   %s' %listrun
            
    listfailed = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.fail' % (wdir,j,j))]
    print 'Failed jobs: %s ' % listfailed
    print '   %s' %listfailed
                        
    if (len(listdone) == len(jobs)):
        print "All jobs successful! Merging output files..."
        os.chdir(wdir)
        #        os.system('ls JOB*/*.root')
        # add automathically only if saved locally
        if (eosoutdir == ''):
            os.system('hadd -f %s.root JOB_*/*.root'%output)
    else:
        for j in listfailed:
            print 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(queue,wdir,j,j,wdir,j,j )
                    


#-----------------------------------------
#--- MAIN
#-----------------------------------------

path = os.getcwd()
conf = path+'/'+options.config
puppiconf = path+'/'+options.puppiConfig
workingdir = path+'/'+options.workdir

eosoutdir = ''

if not options.checkJobs and not options.resubmit:

    # -- write jobs scripts
    if os.path.isdir(workingdir):
     os.system("rm -r "+workingdir)

    os.system("mkdir -p "+workingdir)
    if (options.eosdir !=''):
        eosoutdir = options.eosdir+'/'+options.workdir
        mkdir = 'cmsMkdir '
        command = '%s %s'%(mkdir,eosoutdir)
        print command
        os.system(command)
    writeJobs(workingdir, options.executable, conf, puppiconf, options.inputdir, options.outputname, eosoutdir, options.njobs)
    
    # -- submit jobs 
    if not options.dryRun:
         submitJobs(workingdir, options.njobs, options.queue)
elif options.resubmit and options.jobmin >-1 and options.jobmax >-1:
    for ijob in range(options.jobmax-options.jobmin):
     print 'Resubmitting job %d ' %(ijob+options.jobmin)
     resubcmd = 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(options.queue,workingdir,options.jobmin+ijob,options.jobmin+ijob,workingdir,options.jobmin+ijob,options.jobmin+ijob)
     #print resubcmd
     os.system(resubcmd)

elif options.checkJobs:
    checkJobs(workingdir,options.outputname, options.queue, eosoutdir)
    
