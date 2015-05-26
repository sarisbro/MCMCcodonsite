#!/usr/bin/perl -w

##################################################################################
# Oct. 2004 - Bioinformatics Research Center - North Carolina State University   #
#                                                                                #
# Comments to sarisbro at uottawa dot ca - Thanks.                               #
##################################################################################
# requires: codeml (3.14) in your path and PsTools from www.sysinternals.com (free)
# perl scripts mcmc_041105_Mx.pl && mcmc_041105_M2a.pl in same dir as this one is run
##################################################################################
# to run on SunGridEngine:
# qsub -cwd -S /usr/bin/perl -j y /Volumes/Home/stephane/work/mcmc_codon/M8/1/mcmc_041105_Mx.pl
##################################################################################

use strict;
use Math::Random;            # for random number generator
use File::Copy;              # used to move files

# don't change too much above this line... #######################################
my $NREPS = 100; # 100;
my $MAX_PROCESSES = 20; # for the fork or the "pseudo-fork"
# don't change too much below this line... #######################################

my($OStype,$pwd) = get_pwd();
my $curdir = ''; if($curdir eq ''){$curdir = $pwd;}else{$curdir = $ARGV[0];}
my $PI = 3.141592653; random_set_seed(time(),int(time()/$PI));

my $SEQLEN = 91; my $NBSEQS = 13;
my $M0_tree = "(1: 0.023747, 2: 0.079178, ((3: 0.050407, (4: 0.028430, 5: 0.099890): 0.083461): 0.023648, (6: 0.131681, (7: 0.191149, (8: 0.193336, ((((9: 0.018985, 11: 0.070370): 0.055604, 10: 0.030302): 0.036688, 12: 0.061120): 0.034773, 13: 0.198347): 0.050591): 0.064014): 0.033695): 0.048102): 0.152814);";
my $M0_kappa = 2.47175;
my $M0_rateclass = "2 \n 1.0 0.0 \n 0.90 1.00\n";
my $M3_tree = "(1: 0.025745, 2: 0.088718, ((3: 0.046989, (4: 0.027474, 5: 0.111751): 0.103917): 0.026244, (6: 0.148681, (7: 0.243290, (8: 0.248774, ((((9: 0.020526, 11: 0.073378): 0.058789, 10: 0.033768): 0.039741, 12: 0.065310): 0.003991, 13: 0.270271): 0.069933): 0.068510): 0.026430): 0.054465): 0.175924);";
my $M3_kappa = 2.78820;
my $M3_rateclass = "3 \n 0.28   0.48   0.24 \n 0.00   0.73  3.26";
my $M8_tree = "(1: 0.025625, 2: 0.088759, ((3: 0.046862, (4: 0.026949, 5: 0.112106): 0.104733): 0.026482, (6: 0.148407, (7: 0.244638, (8: 0.249447, ((((9: 0.020551, 11: 0.073273): 0.058703, 10: 0.033842): 0.039294, 12: 0.066265): 0.001793, 13: 0.272319): 0.070910): 0.070110): 0.025919): 0.054930): 0.175460);";
my $M8_kappa = 2.78831;
my $M8_rateclass = "6 \n 0.16  0.16  0.16  0.16  0.16  0.20 \n 0.00  0.05  0.57  0.97  1.00  3.43 \n";
my $M2a_tree = "(1: 0.025591, 2: 0.088874, ((3: 0.046978, (4: 0.026405, 5: 0.112618): 0.105490): 0.026712, (6: 0.148327, (7: 0.246062, (8: 0.250033, ((((9: 0.020605, 11: 0.073256): 0.058693, 10: 0.033955): 0.038819, 12: 0.067259): 0.000645, 13: 0.273579): 0.071600): 0.071541): 0.025435): 0.055482): 0.175161);";
my $M2a_kappa = 2.78555;
my $M2a_rateclass = "3 \n 0.38   0.44   0.18 \n 0.06   1.00   3.63\n";

my $curdir_rep = "$curdir"; chdir "curdir";
my($simul_model,$anal_model,$curmodel);

# list of things to do...
#$simul_model = "M0";  $anal_model = "M2a"; do_anal();
$simul_model = "M0";  $anal_model = "M8";  do_anal();

#$simul_model = "M2a"; $anal_model = "M2a"; do_anal();
#$simul_model = "M2a"; $anal_model = "M8";  do_anal();

#$simul_model = "M3";  $anal_model = "M2a"; do_anal();
#$simul_model = "M3";  $anal_model = "M8";  do_anal();

#$simul_model = "M8";  $anal_model = "M2a"; do_anal();
#$simul_model = "M8";  $anal_model = "M8";  do_anal();


##################################################################################
######################### SUBS ONLY BELOW THIS POINT ############################
##################################################################################
sub CheckNbProcsRunning{
	my $nbrunningjobs = 0;

	if($OStype eq "Darwin"){
		system("qstat -f > $curdir/qstat.out");
		open(PSOUT,"<$curdir/qstat.out") || die "\ncould not read $curdir/qstat.out... $!\n";
		while(<PSOUT>){
			chomp(my $curline = $_);
			if($curline =~ /stephane/){
				$nbrunningjobs++;
			}
		}
		close PSOUT;
		unlink "$curdir/qstat.out";
	}else{
		if($OStype eq "Windows"){
			system "pslist|findstr codeml_old > $curdir_rep/ps.out";
		}else{
			system "ps -u stephane|grep mcmc_041105 > $curdir_rep/ps.out";
		}
		open(PSOUT,"<$curdir_rep/ps.out") || die "\ncould not read $curdir_rep/ps.out... $!\n";
		$nbrunningjobs++ while <PSOUT>;
		close PSOUT;
		unlink "$curdir_rep/ps.out";
	}

	return $nbrunningjobs;
}
##################################################################################
sub CompletedRun{
	my $simulation_path = shift @_;
	my $completed = 0;

	if(-e "$simulation_path/mcmc_time.out"){
		open(TIMEFILE,"<$simulation_path/mcmc_time.out") || warn "\ncould not read $simulation_path/mcmc_time.out... $!\n";
		while(<TIMEFILE>){
			chomp(my $curline = $_);
			if($curline =~ /^MCMC stopped/){
				$completed = 1;
			}
		}
		close TIMEFILE;
	}
	
	return $completed;
}
##################################################################################
sub do_anal{

	my $StartedProcsCount = 0; my($i,$testPassed);
	if($OStype eq "Darwin"){ $MAX_PROCESSES = 9; }

	while($StartedProcsCount < $NREPS+1){
		if(CheckNbProcsRunning() < $MAX_PROCESSES){
			$i = $StartedProcsCount;
			chdir "$curdir_rep/$curmodel/$i";
			$testPassed = CompletedRun("$curdir_rep/$curmodel/$i");
			if($testPassed){
				print "\n$curdir_rep/$curmodel/$i -- all done!";
			}else{
				if($anal_model eq "M2a"){
					if($OStype eq "Windows"){
						system "$curdir_rep/$curmodel/$i/mcmc_041105_M2a.pl &";
					}
					elsif($OStype eq "Darwin"){
						system "qsub -cwd -S /usr/bin/perl -j y $curdir_rep/$curmodel/$i/mcmc_041105_M2a.pl";
					}else{
						system "nice -n10 $curdir_rep/$curmodel/$i/mcmc_041105_M2a.pl &";
					}
				}else{
					if($OStype eq "Windows"){
						system "$curdir_rep/$curmodel/$i/mcmc_041105_Mx.pl &";
					}
					elsif($OStype eq "Darwin"){
						system "qsub -cwd -S /usr/bin/perl -j y $curdir_rep/$curmodel/$i/mcmc_041105_Mx.pl";
					}else{
						system "nice -n10 $curdir_rep/$curmodel/$i/mcmc_041105_Mx.pl &";
					}
				}
			}
			$StartedProcsCount++;
		}
		sleep 3;
	}

}
##################################################################################
sub simul_data{

	my $RNseed = 2*int(random_uniform(1,0,1000000))+1;
	my($my_tree,$my_rateclass,$my_kappa);

	if($simul_model eq "M0"){
		$my_tree = $M0_tree; $my_rateclass = $M0_rateclass; $my_kappa = $M0_kappa;
	}
	elsif($simul_model eq "M2a"){
		$my_tree = $M2a_tree; $my_rateclass = $M2a_rateclass; $my_kappa = $M2a_kappa;
	}
	elsif($simul_model eq "M3"){
		$my_tree = $M3_tree; $my_rateclass = $M3_rateclass; $my_kappa = $M3_kappa;
	}else{
		$my_tree = $M8_tree; $my_rateclass = $M8_rateclass; $my_kappa = $M8_kappa;
	}

	$curmodel = "sim$simul_model-ana$anal_model";
	mkdir("$curmodel",0777);

	open(EVOLVERFILE,">$curdir_rep/MCcodonNSsites.dat") || die "\ncould not create evolver... $!\n";
	print EVOLVERFILE "0            * 0:paml format (mc.paml); 1:paup format (mc.paup)
$RNseed       * random number seed (odd number)
$NBSEQS $SEQLEN $NREPS   * <# seqs>  <# nucleotide sites>  <# replicates>

-1           * <tree length, use -1 if tree has absolute branch lengths>

$my_tree

$my_rateclass

$my_kappa     * kappa

0.01143173  0.00500138  0.02153657  0.00227954
0.00869404  0.00380364  0.01637894  0.00173363
0.01564926  0.00684655  0.00000000  0.00000000
0.00799111  0.00349611  0.00000000  0.00159347
0.00939504  0.00411053  0.01769959  0.00187342
0.00714510  0.00312598  0.01346085  0.00142477
0.01286118  0.00562676  0.02422954  0.00256458
0.00656741  0.00287324  0.01237253  0.00130957
0.03804007  0.01664253  0.07166478  0.00758537
0.02893015  0.01265694  0.05450234  0.00576881
0.05207427  0.02278249  0.09810421  0.01038386
0.02659112  0.01163361  0.05009577  0.00530240
0.01885579  0.00824941  0.03552296  0.00375993
0.01434016  0.00627382  0.02701584  0.00285950
0.02581229  0.01129288  0.04862851  0.00514709
0.01318074  0.00576658  0.02483158  0.00262830

//end of file
";
	close EVOLVERFILE;
	system "evolverNSsites";
	unlink "$curdir_rep/evolver.out","$curdir_rep/ancestral.seq";
	move("$curdir_rep/MCcodonNSsites.dat","$curdir_rep/$curmodel/MCcodonNSsites.dat");

	# create tree file, codeml.ctl and copy scripts for MCMCs
	my $NCATG = 3; if($anal_model eq "M8"){$NCATG = 5;}
	my $curmodelcode = $anal_model; 
	$curmodelcode =~ s/M//;	if($anal_model eq "M2a"){$curmodelcode = 3;}
	for(my $i=1; $i < $NREPS+1; $i++){
		mkdir("$curdir_rep/$curmodel/$i",0777);
		if($anal_model eq "M2a"){
			copy("$curdir/mcmc_041105_M2a.pl","$curdir_rep/$curmodel/$i/mcmc_041105_M2a.pl");
			system "chmod 700 $curdir_rep/$curmodel/$i/mcmc_041105_M2a.pl" if !($OStype eq "Windows");
		}else{
			copy("$curdir/mcmc_041105_Mx.pl","$curdir_rep/$curmodel/$i/mcmc_041105_Mx.pl");
			system "chmod 700 $curdir_rep/$curmodel/$i/mcmc_041105_Mx.pl" if !($OStype eq "Windows");
		}
		open(TREEFILE,">$curdir_rep/$curmodel/$i/data.tre") || die "\ncould not create tree @ $i... $!\n";
		print TREEFILE " $NBSEQS  1\n$my_tree";
		close TREEFILE;
		open(CTLFILE,">$curdir_rep/$curmodel/$i/codeml.ctl") || die "\ncould not create ctl file @ $i... $!\n";
		print CTLFILE "      seqfile = data.nuc 
     treefile = data.tre
      outfile = $anal_model.txt

        noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 2   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        model = 0
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = $curmodelcode  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = $NCATG   * # of categories in the dG or AdG models of rates

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = -1   * 0: simultaneous; 1: one branch at a time
";
		close CTLFILE;
	}

	# parse evolver's outputs
	open(MCFILE,"$curdir_rep/mc.paml") || die "\ncould not open mc.paml... $!\n";
	my $datasetscounter = 0;
	while(<MCFILE>){
		chomp(my $curline = $_);
		if($curline =~ /^ +13/){
			$datasetscounter++;
			open(TMPFILE,">$curdir_rep/$curmodel/$datasetscounter/data.nuc") || die "\ncould not create data @ $datasetscounter... $!\n";
			print TMPFILE "$curline\n";
		}
		elsif(!($curline =~ /^ +13/) && !($curline =~ /^seq/)){
			close TMPFILE;
		}else{
			print TMPFILE "$curline\n";
		}
	}
	close MCFILE;
	unlink "$curdir_rep/mc.paml";

	# move siterates file
	move("$curdir_rep/siterates", "$curdir_rep/$curmodel/siterates");
}
##################################################################################
sub simul_data_tmp{

	my $RNseed = 2*int(random_uniform(1,0,1000000))+1;
	my($my_tree,$my_rateclass,$my_kappa);

	if($simul_model eq "M0"){
		$my_tree = $M0_tree; $my_rateclass = $M0_rateclass; $my_kappa = $M0_kappa;
	}
	elsif($simul_model eq "M2a"){
		$my_tree = $M2a_tree; $my_rateclass = $M2a_rateclass; $my_kappa = $M2a_kappa;
	}
	elsif($simul_model eq "M3"){
		$my_tree = $M3_tree; $my_rateclass = $M3_rateclass; $my_kappa = $M3_kappa;
	}else{
		$my_tree = $M8_tree; $my_rateclass = $M8_rateclass; $my_kappa = $M8_kappa;
	}

	$curmodel = "sim$simul_model-ana$anal_model";
	if(!(-e "$curmodel/siterates")){ # to be able to restart script after crash without doing everything again
		mkdir("$curmodel",0777);

		open(EVOLVERFILE,">$curdir_rep/MCcodonNSsites.dat") || die "\ncould not create evolver... $!\n";
		print EVOLVERFILE "0            * 0:paml format (mc.paml); 1:paup format (mc.paup)
$RNseed       * random number seed (odd number)
$NBSEQS $SEQLEN $NREPS   * <# seqs>  <# nucleotide sites>  <# replicates>

-1           * <tree length, use -1 if tree has absolute branch lengths>

$my_tree

$my_rateclass

$my_kappa     * kappa

0.01143173  0.00500138  0.02153657  0.00227954
0.00869404  0.00380364  0.01637894  0.00173363
0.01564926  0.00684655  0.00000000  0.00000000
0.00799111  0.00349611  0.00000000  0.00159347
0.00939504  0.00411053  0.01769959  0.00187342
0.00714510  0.00312598  0.01346085  0.00142477
0.01286118  0.00562676  0.02422954  0.00256458
0.00656741  0.00287324  0.01237253  0.00130957
0.03804007  0.01664253  0.07166478  0.00758537
0.02893015  0.01265694  0.05450234  0.00576881
0.05207427  0.02278249  0.09810421  0.01038386
0.02659112  0.01163361  0.05009577  0.00530240
0.01885579  0.00824941  0.03552296  0.00375993
0.01434016  0.00627382  0.02701584  0.00285950
0.02581229  0.01129288  0.04862851  0.00514709
0.01318074  0.00576658  0.02483158  0.00262830

//end of file on $OStype
";
		close EVOLVERFILE;
		system "evolverNSsites";
		unlink "$curdir_rep/evolver.out","$curdir_rep/ancestral.seq";
		move("$curdir_rep/MCcodonNSsites.dat","$curdir_rep/$curmodel/MCcodonNSsites.dat");
		
		# parse evolver's outputs
		open(MCFILE,"$curdir_rep/mc.paml") || die "\ncould not open mc.paml... $!\n";
		my $datasetscounter = 0;
		while(<MCFILE>){
			chomp(my $curline = $_);
			if($curline =~ /^ +13/){
				$datasetscounter++;
				open(TMPFILE,">$curdir_rep/$curmodel/$datasetscounter/data.nuc") || die "\ncould not create data @ $datasetscounter... $!\n";
				print TMPFILE "$curline\n";
			}
			elsif(!($curline =~ /^ +13/) && !($curline =~ /^seq/)){
				close TMPFILE;
			}else{
				print TMPFILE "$curline\n";
			}
		}
		close MCFILE;
		unlink "$curdir_rep/mc.paml";

		# move siterates file
		move("$curdir_rep/siterates", "$curdir_rep/$curmodel/siterates");
	}
	# create tree file, codeml.ctl and copy scripts for MCMCs
	my $NCATG = 3; if($anal_model eq "M8"){$NCATG = 5;}
	my $curmodelcode = $anal_model; 
	$curmodelcode =~ s/M//;	if($anal_model eq "M2a"){$curmodelcode = 3;}
	for(my $i=1; $i < $NREPS+1; $i++){
		mkdir("$curdir_rep/$curmodel/$i",0777);
		if($anal_model eq "M2a"){
			copy("$curdir/mcmc_041105_M2a.pl","$curdir_rep/$curmodel/$i/mcmc_041105_M2a.pl");
			system "chmod 700 $curdir_rep/$curmodel/$i/mcmc_041105_M2a.pl" if !($OStype eq "Windows");
		}else{
			copy("$curdir/mcmc_041105_Mx.pl","$curdir_rep/$curmodel/$i/mcmc_041105_Mx.pl");
			system "chmod 700 $curdir_rep/$curmodel/$i/mcmc_041105_Mx.pl" if !($OStype eq "Windows");
		}
		open(TREEFILE,">$curdir_rep/$curmodel/$i/data.tre") || die "\ncould not create tree @ $i... $!\n";
		print TREEFILE " $NBSEQS  1\n$my_tree";
		close TREEFILE;
		open(CTLFILE,">$curdir_rep/$curmodel/$i/codeml.ctl") || die "\ncould not create ctl file @ $i... $!\n";
		print CTLFILE "      seqfile = data.nuc 
     treefile = data.tre
      outfile = $anal_model.txt

        noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 2   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        model = 0
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = $curmodelcode  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = $NCATG   * # of categories in the dG or AdG models of rates

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = -1   * 0: simultaneous; 1: one branch at a time
";
		close CTLFILE;
	}
}
##################################################################################
sub copy_script{

	open(SCRIPTFILE,">$curdir_rep/mcmc_.pl") || die "\ncould not create script file ... $!\n";
	print SCRIPTFILE "";

	close SCRIPTFILE;
	system "chmod 700 $curdir_rep/mcmc_.pl";
}
##################################################################################
sub get_pwd{
	my $OStype = ''; my($this_pwd,@dum);
	# test if OS is Windows
	system "ver > ostype.txt";
	open(TESTFILE,"ostype.txt") || die "\ncould not open ostype.txt(1)... $!\n";
	while(<TESTFILE>){
		chomp(my $curline = $_);
		if($curline =~ /Windows/){
			$OStype = "Windows";
			next;
		}
	}
	close TESTFILE;

	# test if OS is linux/unix, & pretty much assumes this type if you get here
	if($OStype eq ''){
		system "uname > ostype.txt";
		open(TESTFILE,"ostype.txt") || die "\ncould not open ostype.txt(L)... $!\n";
		while(<TESTFILE>){
			chomp($OStype = $_);
		}
		close TESTFILE;
	}
	unlink "ostype.txt";

	if($OStype eq "Windows"){
		system "dir > pwd.txt";
		open(PWDFILE,"pwd.txt") || die "\ncould not open pwd.txt(W)... $!\n";
		while(<PWDFILE>){
			chomp(my $curline = $_);
			if($curline =~ /Directory of/){
				@dum = split(/ of /,$curline);
				$this_pwd = $dum[1];
				next;
			}
		}
		system "cls";
	}else{
		system "pwd > pwd.txt";
		open(PWDFILE,"pwd.txt") || die "\ncould not open pwd.txt(L)... $!\n";
		while(<PWDFILE>){
			chomp($this_pwd = $_);
		}
		system "clear";
	}
	close PWDFILE;
	unlink "pwd.txt";
	$this_pwd =~ s/\\/\//g if $OStype eq "Windows";

	return($OStype,$this_pwd);
}
##################################################################################
##############################   DA *real* END   #################################
##################################################################################
