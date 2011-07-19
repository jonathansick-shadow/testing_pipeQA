import sys, os

homedir = "/lsst/home/becker/lsst_devel/pipeQA_trunk/examples/testPbs"
wwwdir  = "/lsst/home/becker/public_html"
db      = "rplante_PT1_2_u_pt12prod_im2000"

nJobs   = int(sys.argv[1])
for i in range(nJobs):
    jobId   = "Job_%03d" % (i)
    jobDir  = os.path.join(homedir, jobId)
    if not os.path.isdir(jobDir):
        os.makedirs(jobDir)

    outfile = "%s.pbs" % (jobId)
    buff    = open(os.path.join(jobDir, outfile), "w")
    buff.write("#PBS -S /bin/tcsh \n")
    buff.write("#PBS -V \n")
    buff.write("#PBS -d %s \n" % (jobDir))
    buff.write("#PBS -o %s/%s.out \n" % (jobDir, jobId))
    buff.write("#PBS -e %s/%s.err \n" % (jobDir, jobId))
    buff.write(" \n")
    buff.write("source /share/apps/lsst_gcc440/loadLSST.csh \n")
    buff.write("setup testing_pipeQA \n")
    buff.write("setup testing_displayQA \n")
    buff.write("setenv WWW_ROOT %s \n" % (wwwdir))
    buff.write("$TESTING_DISPLAYQA_DIR/bin/newQa.py %s \n" % (db))
    buff.write("setenv WWW_RERUN %s \n" % (db))
    buff.write("$TESTING_PIPEQA_DIR/bin/pipeQa.py --noWwwCache -k -b raft -g 1:%d -v '.*' %s \n" % (i, db))
    buff.close()
