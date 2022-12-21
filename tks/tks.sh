#!/bin/bash

# TODO:
# remove -prefix option; replace with -analyse VT2005_Inf_bernoulli.txt

version=
#version=2010.09 # most recent
#export R_HOME=/usr/contrib/lib64/R
export TKS=${HOME}/apps/tks/$version
script=$TKS/R/tks.R
RSCRIPT=/usr/contrib/bin/Rscript

function help_help () {
  echo "--h | --help    : for list of all options"
  echo " -w | -workflow : basic workflow description"
  echo " -e | -example  : step by step example"
  echo " -u | -unix     : basic unix commands"
  echo "      -version  : how to change version"
}

function help_workflow () {
  echo usage: tks.sh -i input.xml [-model VT2005 -nrep 4 -njobs 32 -script $TKS/R/tks.R -queue blade.q -report 0 -template $TKS/R/tks.Rnw
  echo                    -seed 314159 -keep]
  echo
  echo Basic workflow:
  echo
  echo "1. optimize parameters: tks -input input.xml -model VT2005 -nrep 8 -njobs 32 -log2sde -Inf -log2sdo -1"
  echo "   as a first step we need to find reasonable values of parameters so we start without adding any noise to experimental data (-log2sde -Inf)"
  echo "   you may need to run more extensive search (e.g. -nrep 16 -njobs 64) if initial runs do not produce RMSE values < 1"
  echo "   take a look at produced VT2005_Inf_bernoulli*.txt file, e.g. head VT2005_Inf_bernoulli*.txt, tail VT2005_Inf_bernoulli*.txt"
  echo
  echo "2. analyse results: tks -input input.xml -model VT2005 -nrep 8 -njobs 32 -analyse -prefix FILE (this will generate .pdf for FILE.txt)"
  echo "   this command will generate report pdf file for FILE.txt using ALL rows in the file"
  echo
  echo "3. analyse results: tks -input input.xml -model VT2005 -nrep 8 -njobs 32 -analyse -prefix FILE -cutRMSE 1.0 (this will use ony entries from FILE.txt with RMSE < 1.0)"
  echo "   as a next step we want to use only rows with the lowest RMSE values, you have to decide what cutRMSE to use after looking at RMSE distribution in the report file"
  echo "   generated in the previous step. At this stage it might be useful to rerun calculations with optimal paramater values as a starting point"
  echo "   and with 'ode' solver instead of default 'bernoulli' in order to verify that obtained solution is the same"
  echo "   bernoulli solver is ~x10 faster than ode, but might be unstable and might require time step (dt) adjustment"
  echo "   e.g. parms='c(N0=7.46e+06,Nmax=4.05e+08,Kg=0.593,Kk=13.97,C50k=1.192,H=4.208,beta=37.46,tau=0.0219)' for AB.Amik"
  echo "        tks -input input.xml -model VT2005 -nrep 2 -njobs 8 -parms '\$parms' -log2sde -Inf -log2sdo -8 -solver ode"
  echo
  echo '4. add noise: tks -input input.xml -model VT2005 -nrep 8 -njobs 32 -parms $parms -log2sde -1 -log2sdo -8'
  echo "   with e.g. parms='c(N0=7.46e+06,Nmax=4.05e+08,Kg=0.59,Kk=14.0,C50k=1.19,H=4.21,beta=37.5,tau=0.022)'"
  echo "   at this stage we are adding noise to experimental data to generate distribution of parameter values for the model around optimal values"
  echo "   selected by fitting without noise"
  echo "   you need to select appropriate noise level corresponding to standard deviations in measured log10(CFU) values"
  echo "   in the example input line above we are assuming that StDev(log10(CFU)) is 0.5 = 2^(-1) e.g. on command line: -log2sde -1"
  echo
  echo -version : for help with version changes
  echo -unix    : for help with basic unix commands
  echo --help   : for list of all available options
  echo
  #Rscript $script -help
}

function help_R () {
  $RCRIPT $script -help
}

function help_example () {
  echo
  echo "Step by step example:"
  echo
  echo "# create input xml file with experimental data"
  echo "cp $HOME/apps/tks/data/VT/PA.Levo/PA.Levo.xml . # for the purpose of this example we copy existing file to the current working directory"
  echo "tks -input PA.Levo.xml -m VT2005                         # run parameter fitting for selected model (VT2005 or VT2007); this should take about 3-5 minutes (if cluster is not to busy)"
  echo "head VT2005_Inf_bernoulli.txt                            # take a look at top few lines of file with optimized parameters"
  echo "# the lowest RMSE values should be ~0.67" 
  echo "# if values are much different try again with larger ntry option, e.g. -ntry 256"
  echo "xpdf VT2005_Inf_bernoulli.pdf                            # inspect generated PDF file"
  echo
  echo "parms='c(N0=5.96e+07,Nmax=1.504e+09,Kg=0.49,Kk=28.7,C50k=0.91,H=1.71,beta=148.2,tau=0.018)' # optimized parameter values"
  echo "tks -input PA.Levo.xml -m VT2005 -parms \$parms -solver ode -nrep 1 -njobs 4 -log2sdo -8     # restart run with ode solver to verify convergence"
  echo "head VT2005_Inf_ode.txt                                                                     # params should be close to the starting values"
  echo "tks -input PA.Levo.xml -m VT2005 -parms \$parms -log2sdo -8 -log2sde -2                      # add noise to experimental data"
  echo "xpdf VT2005_2_bernoulli.pdf"
}

function help_unix () {
  echo
  echo Few basic Unix commands:
  echo
  echo pwd : print working directory
  echo "cd DIR : change directory to directory DIR, e.g. cd ../ (one level up) or cd ~/ (home directory)"
  echo ls, ls -lt : to list all files in current directory, to list files in long format sorted by date
  echo cp FILE1 FILE2 : copy FILE1 to FILE2
  echo 'mv FILE1 FILE2 : move (i.e. rename) FILE1 to FILE2'
  echo 'head -4 FILE; more FILE, less FILE : head -4 command shows first 4 lines of the FILE; more/less commands show the whole FILE'
  echo "ls -lt *.txt | head -4" : this will show only 4 most recent files with txt extension
  echo rm FILE1 FILE2, rm *.LOG : remove FILE1 and FILE2, remove all files with LOG extension: WARNING!!! 'rm *' will remove all files in current directory
  echo chmod -w FILE : change FILE status to read-only
  echo 'echo $PATH : shows the current search path for executable programs'
  echo which tks : shows which tks will be used
  echo "CTRL Z; bg : to put command line task in the background"
  echo 'qstat : command to monitor jobs on the cluster (you may need to run 'module load sge' first if qstat is not in your path)'
  echo "qstat -u '*' : will show all currently running jobs on the cluster"
}

function help_version () {
  echo current version=$version
  echo to change version \(e.g.\) : version=2009.11\; export TKS=$HOME/apps/tks/\$version\; export PATH=\$TKS/bin
}

function update () {
  if [ ! -d $TKS ]; then
    mkdir -p $TKS/data/models
    mkdir -p $TKS/{bin,R}
  fi
  cp $HOME/apps/tks/R/tks.{R,Rnw} $TKS/R/
  cp $HOME/apps/tks/bin/tks.sh $TKS/bin/tks
  cp $HOME/apps/tks/data/models/*.opt $TKS/data/models
}

function config ()
{
echo '[TKS Options]'
echo TKS:$TKS
echo RSCRIPT:$RSCRIPT
echo version:$version
grep ':=' $0 | grep '=\$' | sed -e '/^#/d' -e 's/#.*//' -e '/:=}/d' -e 's/^.*${//' -e 's/=//' -e 's/}//'
}

function count_free_cpus () {
  lst=`qstat -f | grep blade | sed -e '/amd64    au/d' -e '/amd64    E/d' -e 's/blade.*BIP//' -e 's/ //g' -e 's/\/.*//'`
  declare -i nfree
  nfree=0
  for x in $lst; do
    nfree=$nfree+$x
  done
  echo $nfree
}

while [ $# -gt 0 ]; do
   case $1 in
       -h|-help)        help_help; exit;;
       --h|--help)      help_R; exit;;
       -u|-unix)        help_unix; exit;;
       -e|-example)     help_example; exit;;
       -w|-workflow)    help_workflow; exit;;
       -version)        help_version; exit;;
       -update)         update; exit;;        # copy current to named version
       -config)         config; exit;;        # create config file
       -m|-model)       model=$2; shift;;
       -p|-prefix)      prefix=$2; shift;;
       -parms)          parms=$2; shift;;     # input values for parameters
       -lower)          lower=$2; shift;;     # lower bound for parameters
       -upper)          upper=$2; shift;;     # upper bound for parameters
       -fixed)          fixed=$2; shift;;     # fixed parameters
       -o|-optimized)   optimized=$2; shift;;
       -i|-input)       input=$2; shift;;
       -j|-jobname)     jobname=$2; shift;;
       -njobs)          njobs=$2; shift;;
       -n|-nrep)        nrep=$2; shift;;
       -ntry)           ntry=$2; shift;;
       -s|-script)      script=$2; shift;;
       -seed)           seed=$2; shift;;      # RNG seed (default=$RANDOM)
       -sleep)          sleep=$2; shift;;
       -solver)         solver=$2; shift;;
       -dt)             dt=$2; shift;;
       -maxit)          maxit=$2; shift;;
       -q|-queue)       queue=$2; shift;;
       -r|-report)      report=$2; shift;;   # 0 - dont generate report, 1 - generate report only (do not optimize)
       -sde)            sde=$2; shift;;
       -sdo)            sdo=$2; shift;;
       -log2sde)        log2sde=$2; shift;;
       -log2sdo)        log2sdo=$2; shift;;
       -t|-template)    template=$2; shift;;
       -c|-cutRMSE)     cutRMSE=$2; shift;;
       -d|-deltaRMSE)   deltaRMSE=$2; shift;;
       -user)           user=$2; shift;;
       -firstName)      firstName=$2; shift;;
       -lastName)       lastName=$2; shift;;
       -a|-analyse)     analyse=1;;
       -clean)          clean=1;;
       -keep)           keep=1;;
       -debug)          debug=1;;
       -v|-verbose)     verbose=$2; shift;;
        * )
             echo "WARNING: unknown command line parameter: $1"; help_help; exit;
   esac
   shift
done

abs () {                               # Function to evaluate the absolute value of a number
if expr $1 + 0 2>/dev/null 1>&2 ; then # Check for numeric input
  if [ $1 -lt 0 ] ; then               # Is the number negative?
    echo `expr 0 - $1`
  else
    echo $1
  fi
  return 0                             # OK
else
  if [ $1 == '-Inf' ]; then
    echo Inf
  else
    return 1                             # Not a number
  fi
fi
}

function pow2 () {             # calculates 2^x
  if [ $1 == '-Inf' ]; then
    echo 0
  else
    echo "e($1*l(2))" | bc -l
  fi
}

model=${model:=VT2005}
prefix=${prefix:=}
nrep=${nrep:=4}
njobs=${njobs:=64}
script=${script:=$TKS/R/tks.R}
template=${template:=$TKS/R/tks.Rnw}
optimized=${optimized:=$TKS/data/models/$model.opt}
queue=${queue:=blade.q}
solver=${solver:=bernoulli} # ode|bernoulli
dt=${dt:=0.015625} # 2^(-6)
maxit=${maxit:=128}
jobname=${jobname:=tks}
log2sde=${log2sde:=-Inf}
sde=${sde:=`pow2 $log2sde`}
log2sdo=${log2sdo:=-1}
sdo=${sdo:=`pow2 $log2sdo`}
deltaRMSE=${deltaRMSE:=0.2}
#seed=${seed:=$RANDOM}
sleep=${sleep:=3}
clean=${clean:=1}
keep=${keep:=}

if [ $parms ]; then optimized=NONE ; fi
#echo parms=$parms
#echo optimized=$optimized

if [ $keep ]; then unset clean; fi
if [ $debug ]; then unset clean; fi

source /usr/contrib/etc/bash.modules
module load modules
module load sge


function write_bash_file ()
{
cat > $1 <<EOF
#!/bin/bash
cd `pwd`
$RSCRIPT $2
EOF
}

function optimize () {
  if [ $ntry ]; then
    nfree=`count_free_cpus`
    #echo optimize: after count_free_cpus: nfree = $nfree
    x=`echo $ntry/$nfree | bc -l`
    #echo after count_free_cpus: x, nfree = $x $nfree
    python_cmd="x=$x; from math import *; x=int(ceil(x)); print x"
    nrep=`python -c "$python_cmd"`
    njobs=$nfree
  fi
  #echo ntry=$ntry
  #echo nrep=$nrep
  #echo njobs=$njobs
  #exit
  declare -i i
  i=0
  while [ $i -lt $njobs ]; do
    i=$i+1
    job_prefix=${process_prefix}_$i
    #echo in optimize: prefix = $prefix
    #echo in optimize: i, job_prefix = $i  $job_prefix
    qsub_sh=q.$job_prefix.sh
    if [ ! $seed ]; then
      qsub_cmd=`echo $cmd -output $job_prefix.txt -prefix $job_prefix -seed $RANDOM`
    else
      qsub_cmd=`echo $cmd -output $job_prefix.txt -prefix $job_prefix -seed $seed`
    fi
    #echo qsub_cmd = $qsub_cmd
    write_bash_file $qsub_sh "$qsub_cmd"
    chmod +x $qsub_sh
    qsub -r yes -p 0 -S /bin/bash -V -q $queue -e `pwd` -o `pwd` -N $jobname.$pid $qsub_sh 1> /dev/null &
  done
  #exit
  echo ...submitted $njobs jobs: jobname=$jobname, pid=$pid
}

function wait_for_jobs () {
  #declare -i nfinished nexpected
  echo ...waiting for jobs to finish: process_prefix=${process_prefix}
  njobs=$1
  sleep $sleep
  while [ 1 ]; do
    nrunning=`qstat -xml | grep $jobname.$$ | wc -l`
    nfinished=`ls ${process_prefix}*.txt 2> /dev/null | wc -l`
    if [ $nrunning -gt 0 ]; then
      #echo njobs=$njobs, nfinished=$nfinished, nrunning=$nrunning
      #echo sleep
      sleep $sleep
    else
      #echo break
      break
    fi
  done
  sleep $sleep
  if [ $njobs != $nfinished ]; then
    echo ERROR: not all jobs finished properly: njobs, nfinished = $njobs, $nfinished
    err=1
  else
    echo ...all jobs finished: njobs = $njobs
  fi
}

function gather () {
  p=$1; x=$2
  #echo p, x = $p $x
  head -1 ${p}_1.$x > $p.$x
  cat ${p}_*.$x | sed -e '/RMSE/d' | sort -n -k 1  >> $p.$x
  echo ...gather: file $p.$x is ready
}

# check for errors

function check_for_errors () {
  declare -i err warn; err=0; warn=0
  if [ ! -s $input ]; then
    echo ERROR: input file: $input : does not exist
    err=$err+1
  fi
  if [ $keep ] && [ $clean ]; then
    echo ERROR: keep and clean are mutually exclusive
    err=$err+1
  fi
  if [ ! -s $optimized ]; then
    warn=$warn+1
    echo WARNING: optimized=$optimized file does not exists!
    optimized=NONE
  fi
  if [ $err -gt 0 ]; then
    exit 1
  fi
}

check_for_errors

#count_free_cpus; echo after check_for_errors: nfree = $nfree

#echo BEFORE: prefix = $prefix
if [ ! $prefix ]; then # if prefix is undefined make it e.g. VT2005_2_ode
  d=`abs $log2sde` # d - for DATA
  p=`abs $log2sdo` # p - for PARAMETERS
  #prefix=${model}_d${d}_p${p}_${solver}
  prefix=${model}_${d}_${solver}
fi

pid=$$
process_prefix=${prefix}_${pid}

#cmd="$script -input $input -output $job_prefix.txt -model $model -prefix $job_prefix -nrep $nrep -sde $sde -sdo $sdo -solver $solver -dt $dt -maxit $maxit -deltaRMSE $deltaRMSE"
cmd="$script -input $input -model $model -nrep $nrep -sde $sde -sdo $sdo -solver $solver -dt $dt -maxit $maxit -deltaRMSE $deltaRMSE -optimized $optimized"
if [ $cutRMSE ]; then cmd=`echo $cmd -cutRMSE $cutRMSE`; fi
if [ $template ]; then cmd=`echo $cmd -template $template`; fi
if [ $parms ]; then cmd=`echo $cmd -parms \"$parms\"`; fi
if [ $lower ]; then cmd=`echo $cmd -lower \"$lower\"`; fi
if [ $upper ]; then cmd=`echo $cmd -upper \"$upper\"`; fi
if [ $fixed ]; then cmd=`echo $cmd -fixed \"$fixed\"`; fi
if [ "$user" ]; then cmd=`echo $cmd -user \"$user\"`; fi
if [ $firstName ]; then cmd=`echo $cmd -firstName \"$firstName\"`; fi
if [ $lastName ]; then cmd=`echo $cmd -lastName \"$lastName\"`; fi

echo version = $version

if [ $analyse ]; then
  cmd=`echo $cmd -output $prefix.txt -prefix $prefix -report 1`
  echo ...reading optimized parameters and generating report pdf
elif [ $report ]; then
  echo ...generating report pdf for optimized paramaters
  cmd=`echo $cmd -report 1`
else
  #echo before optimize: prefix = $prefix
  #count_free_cpus; echo before optimize: nfree = $nfree
  echo RSCRIPT=$RSCRIPT
  #exit
  optimize
  wait_for_jobs $njobs
  if [ $err ]; then
     cat $jobname.$$.e* > error.log
     exit 1
  fi
  gather $process_prefix txt
  mv $process_prefix.txt $prefix.txt
  cmd=`echo $cmd -output $prefix.txt -prefix $prefix -report 1`
  #$RSCRIPT $cmd >& after_optimize.$pid.log
fi
report_sh=report.$pid.sh
report_sh=report.sh
echo ...generating report
echo cmd=$cmd
echo report_sh=$report_sh
write_bash_file $report_sh "$cmd"
chmod +x $report_sh
./report.sh >& report.$pid.log

if [ $clean ]; then
  rm -f ${process_prefix}_*.txt
  rm -f $jobname.$$.[eo]*
  rm -f q.${process_prefix}_*.sh
  rm -f *.$pid.log
  rm -f Rplots.pdf Makefile.$prefix.log
fi

exit

