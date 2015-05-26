#!/usr/bin/perl -w

##################################################################################
# Sep. 2004 - Bioinformatics Research Center - North Carolina State University   #
#                                                                                #
# Comments to sarisbro at uottawa dot ca - Thanks.                               #
##################################################################################
# requires: Math-Random module (0.67 or >), Statistics::GammaDistribution (0.01 or >)
#           and codeml (3.14) in your path

use POSIX;
use Math::Random;
use Statistics::GammaDistribution;
use strict;
use Term::Cap;

my($OStype,$pwd) = get_pwd();
my $curdir = ''; if($curdir eq ''){$curdir = $pwd;}else{$curdir = $ARGV[0];}

my $ctlfile = "$curdir/codeml.ctl";
# my $N_SAMPLES = 100000; my $THINNING = 10; # 100;
my $N_SAMPLES = 10000; my $THINNING = 10; # 100;
my $start_from_mle = 1;
my $set_brlens_2_mles = 1;
my $model_M2a = 0;

##################################################################################
##############################    MAIN    ########################################
##################################################################################

my($seqfile,$treefile,$outfile,$runmode,$seqtype,$CodonFreq,$model,$NSsites,$icode);
my($fix_kappa,$kappa,$fix_omega,$omega,$fix_alpha,$alpha,$Malpha,$ncatG,$RateAncestor);
my($method,$getSE,$clock,$verbose,$noisy,$Nb_of_Seqs);
my $cur_params = ""; my $accepted = 0; my $rejected = 0; my $cycle = 0;
$outfile = ""; my $cur_lnL = 0.0; my $prop_lnL = 0.0; my $PI = 3.141592653;

chdir "$curdir";

random_set_seed(time(),int(time()/$PI));

open(MCMCPARAMSFILE,">$curdir/mcmc_params.out") || die "\ncould not create mcmc_params.out ... $!\n";
open(PROPPARAMFILE,">$curdir/mcmc_time.out") || die "\ncould not create mcmc_prop_params.out ... $!\n";

my $loc_time = localtime(); print "\nLocal time on your $OStype system is $loc_time\n"; # just for the fun of it!
print PROPPARAMFILE "MCMC started on $loc_time\npwd: $curdir\n";

if($start_from_mle){
	print "\nGetting MLEs to initialize MCMC...\n";
}else{
	print "\nInitializing MCMC from random starting points...\n";
}
initialize_mcmc();
print "\nInitialization stage completed.\n";

if($model || $NSsites == 1 || $NSsites == 2 || $NSsites == 4 || $NSsites == 5 || $NSsites == 6 || $NSsites > 8){
	print "\n\n\aThis is only tested for site-models M0, M3, M7 and M8!\n\n";
}elsif($model_M2a && ($NSsites == 3) && !($ncatG == 3)){
	print "\n\n\aRunning model M2a requires that NSsites & ncatG be set to 3 in codeml.ctl!\n\n";
}else{
	my $start_time = (times)[0];
	print "Starting the MCMC...";
	if($NSsites == 3 || $NSsites == 8){
		open(PSSFILE,">>$curdir/mcmc_PostPs.out") || die "\ncould not append to mcmc_PostPs.out ... $!\n";
		open(PostMeansFILE,">>$curdir/mcmc_PostMeans.out") || die "\ncould not append to mcmc_PostMeans.out ... $!\n";
	}
	while($accepted <= $N_SAMPLES){
		($cur_params,$cur_lnL) = do_the_mcmc($cur_params,$cur_lnL);
	}
	if($NSsites == 3 || $NSsites == 8){
		close PSSFILE; close PostMeansFILE;
	}
	print " ok.\n";
	print "Sampled $N_SAMPLES states every $THINNING steps of the MCMC.\n\a";
	print "                 ++++ CHECK CONVERGENCE ++++";
	my $end_time = (times)[0];
	my $running_time = $end_time - $start_time;
	print "\n$running_time seconds initialization time\n";
	print PROPPARAMFILE "This MCMC took $running_time seconds to run.\n";
	$loc_time = localtime(); print PROPPARAMFILE "MCMC stopped on $loc_time\n";
}

close MCMCPARAMSFILE;
close PROPPARAMFILE;

##################################################################################
##############################    SUBS    ########################################
##################################################################################
sub update_brlen{
	my $my_cur_params = shift @_;
	my $my_prop_params = ""; my $my_log_ratio = 0.0; my $interval_width = 0.5;# 0.5;
	my $min_brlen = 0.00001; my $max_brlen = 100.0;  # for unif[.001,100] on brlen
	my $ExpDistParam = 10.0;

	my @my_cur_params = split(/\s+/,$my_cur_params);

	my $cur_i = $my_cur_params[my $param_id = int(random_uniform(1,0,2*$Nb_of_Seqs-3))];
	my $prop_i = abs($cur_i + (random_uniform() * $interval_width) - ($interval_width/2.0));
	if($prop_i>$max_brlen){$prop_i = $max_brlen - ($prop_i - $max_brlen);}
	if($prop_i<$min_brlen){$prop_i = ($min_brlen-$prop_i) + $min_brlen;}
	$my_cur_params[$param_id] = $prop_i;
	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}
	$my_log_ratio += 0.0; # proposal ratio (unif distr)
	# $my_log_ratio += 0.0; # prior ratio (unif distr)
	$my_log_ratio += ($ExpDistParam * $cur_i - $ExpDistParam * $prop_i);; # prior ratio (exp(10) distr)

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub update_kappa{
	my $my_cur_params = shift @_;
	my $my_prop_params = "";     my $my_log_ratio = 0.0;    my $interval_width = 1.0;
	my $min_kappa = 0.01; my $max_kappa = 100.0;  # unif[.001,100] on kappa

	my @my_cur_params = split(/\s+/,$my_cur_params);

	my $cur_i = $my_cur_params[2*$Nb_of_Seqs-3];
	my $prop_i = abs($cur_i + (random_uniform() * $interval_width) - ($interval_width/2.0));
	if($prop_i>$max_kappa){$prop_i = $max_kappa - ($prop_i - $max_kappa);}
	if($prop_i<$min_kappa){$prop_i = ($min_kappa-$prop_i) + $min_kappa;}
	$my_cur_params[2*$Nb_of_Seqs-3] = $prop_i;
	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}
	$my_log_ratio += 0.0; # proposal ratio (unif distr)
	$my_log_ratio += 0.0; # prior ratio (unif distr)

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub update_omega{
	my $my_cur_params = shift @_;
	my $my_prop_params = "";     my $my_log_ratio = 0.0;    my $interval_width = 0.5; # 1.0;
	my $min_omega = 0.01; my $max_omega = 100.0;  # unif[.001,50] on omega
	my($cur_i,$prop_i,$whichOmega2update);

	my @my_cur_params = split(/\s+/,$my_cur_params);

	if($NSsites == 0){
		$cur_i = $my_cur_params[2*$Nb_of_Seqs-2];
		$prop_i = abs($cur_i + (random_uniform() * $interval_width) - ($interval_width/2.0));
		if($prop_i>$max_omega){$prop_i = $max_omega - ($prop_i - $max_omega);}
		if($prop_i<$min_omega){$prop_i = ($min_omega-$prop_i) + $min_omega;}
		$my_cur_params[2*$Nb_of_Seqs-2] = $prop_i;
		$my_log_ratio += 0.0; # proposal ratio (unif distr)
		$my_log_ratio += 0.0; # prior ratio (unif distr)
	}elsif($NSsites == 3){
		$whichOmega2update = int( random_uniform()*$ncatG );
		if($whichOmega2update == ($ncatG-1)){
			$interval_width *= 10.0;# 5.0;
		}
		elsif(!$whichOmega2update){
			$interval_width /= 1.;# 5.0;
		}
		if($model_M2a && $whichOmega2update == 1){
			$prop_i = 1.0;
		}else{
			$cur_i = $my_cur_params[(2*$Nb_of_Seqs-3)+1+($ncatG-1)+$whichOmega2update];
			$prop_i = abs($cur_i + (random_uniform() * $interval_width) - ($interval_width/2.0));
			# $prop_i = abs(random_normal(1,$cur_i,$interval_width));
			if(!$whichOmega2update){ # rate category 0, the smallest
				$max_omega = $my_cur_params[(2*$Nb_of_Seqs-3)+1+($ncatG-1)+($whichOmega2update+1)];
			}elsif($whichOmega2update == ($ncatG-1)){ # largest rate category
				$min_omega = $my_cur_params[(2*$Nb_of_Seqs-3)+1+($ncatG-1)+($whichOmega2update-1)];
			}else{ # all that inbetween
				$max_omega = $my_cur_params[(2*$Nb_of_Seqs-3)+1+($ncatG-1)+($whichOmega2update+1)];
				$min_omega = $my_cur_params[(2*$Nb_of_Seqs-3)+1+($ncatG-1)+($whichOmega2update-1)];
			}
			if($prop_i>$max_omega){$prop_i = $max_omega - ($prop_i - $max_omega);}
			if($prop_i<$min_omega){$prop_i = ($min_omega-$prop_i) + $min_omega;}
		}
		$my_cur_params[(2*$Nb_of_Seqs-3)+1+($ncatG-1)+$whichOmega2update] = $prop_i;

		# $my_log_ratio += lnNormalPDF($prop_i,$cur_i,$interval_width)-lnNormalPDF($cur_i,$prop_i,$interval_width); # proposal ratio (normal distr)
		$my_log_ratio += 0.0; # proposal ratio (unif distr)
		$my_log_ratio += 0.0; # prior ratio (unif distr)
	}elsif($NSsites == 8){
		$interval_width *= 5.0;
		$cur_i = $my_cur_params[2*$Nb_of_Seqs+1];
		$prop_i = abs($cur_i + (random_uniform() * $interval_width) - ($interval_width/2.0));
		if($prop_i>$max_omega){$prop_i = $max_omega - ($prop_i - $max_omega);}
		if($prop_i<$min_omega){$prop_i = ($min_omega-$prop_i) + $min_omega;}
		$my_cur_params[2*$Nb_of_Seqs+1] = $prop_i;
		$my_log_ratio += 0.0; # proposal ratio (unif distr)
		$my_log_ratio += 0.0; # prior ratio (unif distr)
	}else{
		last;
	}
	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub update_omegaFreqs{ # for M3 & M2a
	my $my_cur_params = shift @_;
	my $my_prop_params = "";  my $my_log_ratio = 0.0;  my $omegaFtune = 0.2;
	my $min_freq = 0.0001; my $max_freq = 0.9999;
	my $eps = 0.0000001; my $alpha_curFreq;
	my($sum_freqs,$class,$x,$y,@cur_freqs,@prop_freqs,@alpha,@dirichlet_params);
	my $priorPi = $ncatG;
	my($lhr,$cx,$cy);

	my @my_cur_params = split(/\s+/,$my_cur_params);

	# there are $ncatG-1 indep freqs, starting from 2*$Nb_of_Seqs-2
	$sum_freqs = 0.0;
	for($class=0;$class<$ncatG-1;$class++){
		$cur_freqs[$class] = $my_cur_params[2*$Nb_of_Seqs-2+$class];
		$sum_freqs += $cur_freqs[$class];
	}
	$cur_freqs[$ncatG-1] = 1.0 - $sum_freqs;

	# propose freqs from Dirichlet distribution
	my $g = Statistics::GammaDistribution->new();
	for(my $class=0;$class<$ncatG;$class++){
		$alpha[$class] = 10;# 3000; # 50.0;
		$alpha_curFreq = $alpha[$class] * $cur_freqs[$class];
		if($alpha_curFreq < $eps){$alpha_curFreq = $eps;}
		$dirichlet_params[$class] = $alpha_curFreq;
	}
	@prop_freqs = $g->dirichlet_dist(@dirichlet_params);

	for(my $class=0;$class<$ncatG;$class++){
		if($prop_freqs[$class]>$max_freq){$prop_freqs[$class] = $max_freq - ($prop_freqs[$class] - $max_freq);}
		if($prop_freqs[$class]<$min_freq){$prop_freqs[$class] = ($min_freq-$prop_freqs[$class]) + $min_freq;}
	}
	$sum_freqs = 0.0;
	for(my $class=0;$class<$ncatG;$class++){
		$sum_freqs += $prop_freqs[$class];
	}
	for(my $class=0;$class<$ncatG;$class++){
		$prop_freqs[$class] /= $sum_freqs;
	}

	for($class=0;$class<$ncatG-1;$class++){
		$my_cur_params[2*$Nb_of_Seqs-2+$class] = $prop_freqs[$class];
	}
	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}

	# proposal ratio: Dirichlet distribution
	$lhr = 0.0;
	for($class=0;$class<$ncatG;$class++){
		$cx = $alpha[$class] * $cur_freqs[$class];
		$cy = $alpha[$class] * $prop_freqs[$class];
		$lhr += ($cy-1.0)*log($cur_freqs[$class]) - ($cx-1.0)*log($prop_freqs[$class]) - lgamma($cy) + lgamma($cx);
	}
	$my_log_ratio += $lhr;

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub update_p0{
	my $my_cur_params = shift @_;
	my $my_prop_params = "";     my $my_log_ratio = 0.0;    my $interval_width = 0.2;
	my $min_freq = 0.0001; my $max_freq = 0.9999;
	my($sum_freqs,$class,$x,$y,@cur_freqs,@prop_freqs,@alpha,@dirichlet_params);
	my $priorPi = 2.0; my $eps = 0.0000001; my $alpha_curFreq;
	my($lhr,$cx,$cy);

	my @my_cur_params = split(/\s+/,$my_cur_params);

	# parameter p0: frequency of beta(p,q) in M8
	$cur_freqs[0] = $my_cur_params[2*$Nb_of_Seqs-2];
	$cur_freqs[1] = 1.0 - $cur_freqs[0];

	# propose freqs from Dirichlet distribution
	my $g = Statistics::GammaDistribution->new();
	for(my $class=0;$class<2;$class++){
		$alpha[$class] = 10.;# 2000.0; # 50.0; # 1.0;
		if($cur_freqs[$class] < $eps){$cur_freqs[$class] = $eps;}
		$alpha_curFreq = $alpha[$class] * $cur_freqs[$class];
		if($alpha_curFreq < $eps){$alpha_curFreq = $eps;}
		$dirichlet_params[$class] = $alpha_curFreq;
	}
	@prop_freqs = $g->dirichlet_dist(@dirichlet_params);
	for(my $class=0;$class<2;$class++){
		if($prop_freqs[$class]>$max_freq){$prop_freqs[$class] = $max_freq - ($prop_freqs[$class] - $max_freq);}
		if($prop_freqs[$class]<$min_freq){$prop_freqs[$class] = ($min_freq-$prop_freqs[$class]) + $min_freq;}
	}
	$sum_freqs = 0.0;
	for(my $class=0;$class<2;$class++){
		$sum_freqs += $prop_freqs[$class];
	}
	for(my $class=0;$class<2;$class++){
		$prop_freqs[$class] /= $sum_freqs;
	}

	$my_cur_params[2*$Nb_of_Seqs-2] = $prop_freqs[0];
	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}

	# proposal ratio: Dirichlet distribution
	$lhr = 0.0;
	for($class=0;$class<2;$class++){
		$cx = $alpha[$class] * $cur_freqs[$class];
		$cy = $alpha[$class] * $prop_freqs[$class];
		$lhr += ($cy-1.0)*log($cur_freqs[$class]) - ($cx-1.0)*log($prop_freqs[$class]) - lgamma($cy) + lgamma($cx);
	}
	$my_log_ratio += $lhr;

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub update_p{
	my $my_cur_params = shift @_;
	my $my_prop_params = "";     my $my_log_ratio = 0.0;    my $interval_width = 10.; # 0.5;
	my $min_pq = 0.01; my $max_pq = 15.0;  # unif[.001,15] on pq
	my($cur_i,$randu,$prop_i);

	my @my_cur_params = split(/\s+/,$my_cur_params);

	# parameter p of beta(p,q)
	if($NSsites == 7){
		$cur_i = $my_cur_params[2*$Nb_of_Seqs-2];
		$randu = (random_uniform() * $interval_width) - ($interval_width/2.0);
		$prop_i = $cur_i + $randu;
		if($prop_i>$max_pq){$prop_i = $max_pq - ($prop_i - $max_pq);}
		if($prop_i<$min_pq){$prop_i = ($min_pq-$prop_i) + $min_pq;}
		$my_cur_params[2*$Nb_of_Seqs-2] = $prop_i;
		$my_log_ratio += 0.; # proposal ratio (unif distr)
		$my_log_ratio += 0.; # prior ratio (unif distr)
	}elsif($NSsites == 8){
		$cur_i = $my_cur_params[2*$Nb_of_Seqs-1];
		$randu = (random_uniform() * $interval_width) - ($interval_width/2.0);
		$prop_i = $cur_i + $randu;
		if($prop_i>$max_pq){$prop_i = $max_pq - ($prop_i - $max_pq);}
		if($prop_i<$min_pq){$prop_i = ($min_pq-$prop_i) + $min_pq;}
		$my_cur_params[2*$Nb_of_Seqs-1] = $prop_i;
		$my_log_ratio += 0.0; # proposal ratio (unif distr)
		$my_log_ratio += 0.0; # prior ratio (unif distr)
	}

	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub update_q{
	my $my_cur_params = shift @_;
	my $my_prop_params = "";     my $my_log_ratio = 0.0;    my $interval_width = 10.; # 0.5;
	my $min_pq = 0.01; my $max_pq = 15.0;  # unif[.001,15] on pq
	my($cur_i,$randu,$prop_i);

	my @my_cur_params = split(/\s+/,$my_cur_params);

	# parameter q of beta(p,q)
	if($NSsites == 7){
		$cur_i = $my_cur_params[2*$Nb_of_Seqs-1];
		$randu = (random_uniform() * $interval_width) - ($interval_width/2.0);
		$prop_i = $cur_i + $randu;
		if($prop_i>$max_pq){$prop_i = $max_pq - ($prop_i - $max_pq);}
		if($prop_i<$min_pq){$prop_i = ($min_pq-$prop_i) + $min_pq;}
		$my_cur_params[2*$Nb_of_Seqs-1] = $prop_i;
		$my_log_ratio += 0.0; # proposal ratio (unif distr)
		# $my_log_ratio += $ExpPriorParam*$cur_i - $ExpPriorParam*$prop_i; # proposal ratio (exp distr)
		$my_log_ratio += 0.0; # prior ratio (unif distr)
	}elsif($NSsites == 8){
		$cur_i = $my_cur_params[2*$Nb_of_Seqs];
		$randu = (random_uniform() * $interval_width) - ($interval_width/2.0);
		$prop_i = $cur_i + $randu;
		if($prop_i>$max_pq){$prop_i = $max_pq - ($prop_i - $max_pq);}
		if($prop_i<$min_pq){$prop_i = ($min_pq-$prop_i) + $min_pq;}
		$my_cur_params[2*$Nb_of_Seqs] = $prop_i;
		$my_log_ratio += 0.0; # proposal ratio (unif distr)
		# $my_log_ratio += $ExpPriorParam*$cur_i - $ExpPriorParam*$prop_i; # proposal ratio (exp distr)
		$my_log_ratio += 0.0; # prior ratio (unif distr)
	}

	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub do_the_mcmc{
	my $kur_params = shift @_;
	my $kur_lnL = shift @_;
	my($prop_params,$prop_prior_log_ratio,$alphaMH,$rndu);

	# propose parameters
	if($NSsites == 0){
		if($cycle == 1){
			($prop_params,$prop_prior_log_ratio) = update_kappa($kur_params);
		}elsif($cycle == 2  && !$set_brlens_2_mles){
			($prop_params,$prop_prior_log_ratio) = update_brlen($kur_params);
		}else{
			($prop_params,$prop_prior_log_ratio) = update_omega($kur_params);
			$cycle = 0;
		}
	}elsif($NSsites == 3){
		if($cycle == 0){
			($prop_params,$prop_prior_log_ratio) = update_omega($kur_params);
		}elsif($cycle == 1){
			($prop_params,$prop_prior_log_ratio) = update_omegaFreqs($kur_params);
		}elsif($cycle == 2  && !$set_brlens_2_mles){
			($prop_params,$prop_prior_log_ratio) = update_brlen($kur_params);
		}else{
			($prop_params,$prop_prior_log_ratio) = update_kappa($kur_params);
			$cycle = 0;
		}
	}elsif($NSsites == 7){
		if($cycle  == 1){
			($prop_params,$prop_prior_log_ratio) = update_kappa($kur_params);
		}elsif($cycle == 2){
			($prop_params,$prop_prior_log_ratio) = update_p($kur_params);
		}elsif($cycle == 3  && !$set_brlens_2_mles){
			($prop_params,$prop_prior_log_ratio) = update_brlen($kur_params);
		}else{
			($prop_params,$prop_prior_log_ratio) = update_q($kur_params);
			$cycle = 0;
		}
	}elsif($NSsites == 8){
		if($cycle == 0){
			($prop_params,$prop_prior_log_ratio) = update_omega($kur_params);
		}elsif($cycle == 1){
			($prop_params,$prop_prior_log_ratio) = update_kappa($kur_params);
		}elsif($cycle == 2){
			($prop_params,$prop_prior_log_ratio) = update_p0($kur_params);
		}elsif($cycle == 3){
			($prop_params,$prop_prior_log_ratio) = update_q($kur_params);
		}elsif($cycle == 4  && !$set_brlens_2_mles){
			($prop_params,$prop_prior_log_ratio) = update_brlen($kur_params);
		}else{
			($prop_params,$prop_prior_log_ratio) = update_p($kur_params);
			$cycle = 0;
		}
	}else{
		last;
	}

	open(INFILE,">$curdir/in.codeml") || die "\ncould not create in.codeml ... $!\n";
	print INFILE "-1 $prop_params";
	close INFILE;
	# compute lnL (run codeml) & get results from outfile
	if($OStype eq "Darwin"){
		system "/Volumes/Home/stephane/bin/codeml_old $ctlfile > screen.out";
	}else{
		system "codeml_old $ctlfile > screen.out";
	}
	($prop_params,$prop_lnL) = get_params("$outfile");

	if(!$prop_lnL){next;}

	# Hastings ratio: accept / reject proposed state
	my $sum = ($prop_lnL - $kur_lnL) + $prop_prior_log_ratio;
	if($sum < -50.0){
		$alphaMH = 0.0;
	}elsif($sum > 0.0){
		$alphaMH = 1.0
	}else{
		$alphaMH = exp($sum);
	}

	if(random_uniform() < $alphaMH){
		$accepted ++;
		$kur_lnL = $prop_lnL; $kur_params = $prop_params; $kur_params =~ s/^\s+//;
		$cycle++;
		if($accepted % $THINNING == 0){
			get_post_proba_of_PSS($accepted) if ($NSsites == 3 || $NSsites == 8);
			print MCMCPARAMSFILE "$accepted\t$prop_lnL $prop_params\n";
			my $cur_accept_rate = int(100.0 * ($accepted/($accepted+$rejected))) ;
			print "\n#$accepted\t$kur_lnL\t$kur_params ($cur_accept_rate%)";
		} # thinning
	}else{$rejected++;}
	return($kur_params,$kur_lnL);
}
##################################################################################
sub initialize_mcmc{
	my($init_rand,$params,$lnL,$j); my $my_rand_params = "";
	# read control file
	open(CTLFILE,"$ctlfile") || die "\ncould not open $ctlfile ... $!\n";
	while(<CTLFILE>){
		chomp(my $curline = $_);
		if($curline =~ /^ +seqfile/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $seqfile = $tmp2[0];
			}else{$seqfile = $tmp1[1];}
		}
		elsif($curline =~ /^ +treefile/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $treefile = $tmp2[0];
			}else{$treefile = $tmp1[1];}
		}
		elsif($curline =~ /^ +outfile/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $outfile = $tmp2[0];
			}else{$outfile = $tmp1[1];}
		}
		elsif($curline =~ /^ +noisy/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $noisy = $tmp2[0];
			}else{$noisy = $tmp1[1];}
		}
		elsif($curline =~ /^ +verbose/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $verbose = $tmp2[0];
			}else{$verbose = $tmp1[1];}
		}
		elsif($curline =~ /^ +runmode/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $runmode = $tmp2[0];
			}else{$runmode = $tmp1[1];}
		}
		elsif($curline =~ /^ +seqtype/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $seqtype = $tmp2[0];
			}else{$seqtype = $tmp1[1];}
		}
		elsif($curline =~ /^ +CodonFreq/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $CodonFreq = $tmp2[0];
			}else{$CodonFreq = $tmp1[1];}
		}
		elsif($curline =~ /^ +model/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $model = $tmp2[0];
			}else{$model = $tmp1[1];}
		}
		elsif($curline =~ /^ +NSsites/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $NSsites = $tmp2[0];
			}else{$NSsites = $tmp1[1];}
		}
		elsif($curline =~ /^ +icode/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $icode = $tmp2[0];
			}else{$icode = $tmp1[1];}
		}
		elsif($curline =~ /^ +fix_kappa/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $fix_kappa = $tmp2[0];
			}else{$fix_kappa = $tmp1[1];}
		}
		elsif($curline =~ /^ +kappa/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $kappa = $tmp2[0];
			}else{$kappa = $tmp1[1];}
		}
		elsif($curline =~ /^ +fix_omega/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $fix_omega = $tmp2[0];
			}else{$fix_omega = $tmp1[1];}
		}
		elsif($curline =~ /^ +omega/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $omega = $tmp2[0];
			}else{$omega = $tmp1[1];}
		}
		elsif($curline =~ /^ +fix_alpha/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $fix_alpha = $tmp2[0];
			}else{$fix_alpha = $tmp1[1];}
		}
		elsif($curline =~ /^ +alpha/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $alpha = $tmp2[0];
			}else{$alpha = $tmp1[1];}
		}
		elsif($curline =~ /^ +Malpha/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $Malpha = $tmp2[0];
			}else{$Malpha = $tmp1[1];}
		}
		elsif($curline =~ /^ +ncatG/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $ncatG = $tmp2[0];
			}else{$ncatG = $tmp1[1];}
		}
		elsif($curline =~ /^ +clock/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $clock = $tmp2[0];
			}else{$clock = $tmp1[1];}
		}
		elsif($curline =~ /^ +getSE/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $getSE = $tmp2[0];
			}else{$getSE = $tmp1[1];}
		}
		elsif($curline =~ /^ +RateAncestor/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $RateAncestor = $tmp2[0];
			}else{$RateAncestor = $tmp1[1];}
		}
		elsif($curline =~ /^ +method/){
			$curline =~ s/\s//g;
			my @tmp1 = split(/=/,$curline);
			if($tmp1[1] =~ /\*/){
				my @tmp2 = split(/\*/,$tmp1[1]); $method = $tmp2[0];
			}else{$method = $tmp1[1];}
		}
	}
	close CTLFILE;
	# print "\n$seqfile\n$treefile\n$outfile\n$runmode\n$seqtype\n$CodonFreq\n$model\n$NSsites\n$icode\n$fix_kappa\n$kappa\n$fix_omega\n$omega\n$fix_alpha\n$alpha\n$Malpha\n$ncatG";
	
	$Nb_of_Seqs = Get_ns("$curdir/$seqfile");
	
	if($start_from_mle){
		# run codeml to get "good" starting values
		unlink "$curdir/in.codeml" if -e "$curdir/in.codeml";
		if($OStype eq "Darwin"){
			system "/Volumes/Home/stephane/bin/codeml_old $ctlfile";
		}else{
			system "codeml_old $ctlfile";
		}

		# write in.codeml file with starting values
		open(INFILE,">$curdir/in.codeml") || die "\ncould not create in.codeml ... $!\n";
		($params,$lnL) = get_params("$outfile");
		print INFILE "-1 $params";
		close INFILE;
	}else{ # start from random starting point
		# branch lengths
		for($j=0;$j<2*$Nb_of_Seqs-3;$j++){
			$init_rand = random_uniform(1,0.01,1);
			$my_rand_params .= "$init_rand\t";
		}
		# kappa
		$init_rand = random_uniform(1,0.01,5);
		$my_rand_params .= "$init_rand\t";
		# other parameters, model dependent...
		if($NSsites == 3){ # for M3
			# omega freqs
			for($j=0;$j<$ncatG-1;$j++){
				$init_rand = 1./$ncatG;
				$my_rand_params .= "$init_rand\t";
			}
			# omega vals
			for ($j=0;$j<$ncatG;$j++){
				$my_rand_params .= "1.00\t";
			}
		}elsif($NSsites == 8){ # for M8
			# p0, freq of beta(p,q)
			$init_rand = random_uniform(1,0.01,0.99);
			$my_rand_params .= "$init_rand\t";
			# p and q
			$init_rand = random_uniform(1,0.05,2);
			$my_rand_params .= "$init_rand\t";
			$init_rand = random_uniform(1,0.05,2);
			$my_rand_params .= "$init_rand\t";
			# omega
			$init_rand = random_uniform(1,0.25,2);
			$my_rand_params .= "$init_rand\t";
		}elsif($NSsites == 7){ # for M7
			# p and q
			$init_rand = random_uniform(1,0.05,2);
			$my_rand_params .= "$init_rand\t";
			$init_rand = random_uniform(1,0.05,2);
			$my_rand_params .= "$init_rand\t";
			# omega
			$init_rand = random_uniform(1,0.25,2);
			$my_rand_params .= "$init_rand\t";
		}else{ # for M0, it's only omega
			$init_rand = random_uniform(1,0.5,2);
			$my_rand_params .= "$init_rand\t";
		}

		open(INFILE,">$curdir/in.codeml") || die "\ncould not create in.codeml ... $!\n";
		print INFILE "-1 $my_rand_params";
		close INFILE;
		if($OStype eq "Darwin"){
			system "/Volumes/Home/stephane/bin/codeml_old $ctlfile";
		}else{
			system "codeml_old $ctlfile";
		}
		($params,$lnL) = get_params("$outfile");
	}

	# other initializations
	$cur_lnL = $lnL;
	$cur_params = $params;
	$cur_params =~ s/^\s+//;
	if($NSsites == 3 || $NSsites == 8){
		open(PSSFILE,">$curdir/mcmc_PostPs.out") || die "\ncould not create mcmc_PostPs.out ... $!\n";
		open(PostMeansFILE,">$curdir/mcmc_PostMeans.out") || die "\ncould not create to mcmc_PostMeans.out ... $!\n";
		get_post_proba_of_PSS(0);
		close PSSFILE; close PostMeansFILE;
	}
	print MCMCPARAMSFILE "0\t$lnL $params\n";
}
#################################################################################
sub Get_ns{
	my $infile = shift @_;

	open(MYDNATREE,"$infile") || warn "\ncould not open $infile ... $!\n";
	chomp (my $line = <MYDNATREE>);
	$line =~ s/^\s+//;
	my($ns,$ls)=split(/\s+/, $line);
	close MYDNATREE;
	return $ns;
}
##################################################################################
sub get_post_proba_of_PSS{
	my $mcmc_step = shift @_;
	my $proba_list = "";
	my $post_means_list = "";
	my $do_line_after_next_one = 0; my $do_next_line = 0; my @this_line;
	open(RSTFILE,"$curdir/rst") || die "\ncould not open rst file ... $!\n";
	# my $nb_of_params = scalar(@my_cur_params);
	while(<RSTFILE>){
		chomp(my $curline = $_);
		if($curline =~ /used as reference/){
			$do_line_after_next_one = 1;
		}
		elsif($do_line_after_next_one){
			$do_line_after_next_one = 0;
			$do_next_line = 1;
		}
		elsif($do_next_line){
			$curline =~ s/^\s+//;
			@this_line = split(/\s+/,$curline);
			if(scalar(@this_line) > 5){
				$proba_list .= "$this_line[scalar(@this_line)-1]\t";
				$post_means_list .= "$this_line[scalar(@this_line)-2]\t";
			}
		}
		elsif($curline =~ /Positively selected sites/){
			$do_next_line = 0;
		}
	}
	print PSSFILE "$mcmc_step\t$proba_list\n";
	print PostMeansFILE "$mcmc_step\t$post_means_list\n";
	close RSTFILE;
}
##################################################################################
sub get_params{
	my $res_file = shift @_;
	my $my_lnL = my $skip_next_line = my $get_next_line = 0; my $my_params;
	my $line_found = 0;
	open(RESFILE,"$curdir/$res_file") || die "\ncould not open $res_file ... $!\n";
	while(<RESFILE>){
		chomp(my $curline = $_);
		if($curline =~ /^lnL/){
			$line_found = 1;
			my @tmp1 = split(/\:/,$curline);
			my @tmp2 = split(/\+/,$tmp1[3]);
			$tmp2[0] =~ s/\s//g;
			$my_lnL = $tmp2[0];
			$skip_next_line = 1;
		}
		elsif($skip_next_line){
			$get_next_line = 1;
			$skip_next_line = 0;
		}
		elsif($get_next_line){
			$my_params = $curline;
			$get_next_line = 0;
		}
	}
	if(!$line_found){
		print "\n\ncodeml failed! MCMC step jumped.\n\n\a";
		last;
		$my_lnL = 0;
	}
	close RESFILE;
	$my_params =~ s/\s+/  /g;
	# print PROPPARAMFILE "$my_params\n";
	return ($my_params,$my_lnL);
}
##################################################################################
sub my_min{
	my $a = shift @_; my $b = shift @_; my $min;
	if($a<$b){$min = $a;}else{$min = $b;}
    return $min;
}
##################################################################################
sub my_max{
	my $a = shift @_; my $b = shift @_; my $max;
	if($a>$b){$max = $a;}else{$max = $b;}
    return $max;
}
##################################################################################
sub lnNormalPDF{
	my $x = shift @_; my $mu = shift @_; my $var = shift @_;
	return -($x-$mu)*($x-$mu)/(2.*$var)-.5*log(2.*$PI*$var);
}
##################################################################################
sub rescale_brlen{
	my $my_cur_params = shift @_;
	my $my_prop_params = "";     my $my_log_ratio = 0.0;
	my $min_brlen = 0.00001; my $max_brlen = 100.0;  # unif[.001,100] on brlen
	my $ExpPriorParam = 1.0/10.0;                     # for exp() prior
	my $tuning = 0.2;

	my @my_cur_params = split(/\s+/,$my_cur_params);
	my $lenFactor = exp($tuning * (random_uniform() - 0.5));

	my $cur_i = $my_cur_params[my $param_id = int(random_uniform(1,0,2*$Nb_of_Seqs-3))];
	my $prop_i = $cur_i * $lenFactor;
	if($prop_i>$max_brlen){$prop_i = $max_brlen - ($prop_i - $max_brlen);}
	if($prop_i<$min_brlen){$prop_i = ($min_brlen-$prop_i) + $min_brlen;}
	foreach my $i (@my_cur_params){
		$my_prop_params .= "$i  ";
	}
	$my_log_ratio += 0.0; # proposal ratio (unif distr)
	# $my_log_ratio += $ExpPriorParam*$cur_i - $ExpPriorParam*$prop_i; # proposal ratio (exp() distr)
	$my_log_ratio += 0.0; # prior ratio (unif distr)

	return ($my_prop_params,$my_log_ratio);
}
##################################################################################
sub lgamma {  # per code from numerical recipies
  my $xx = $_[0];
  my ($j, $ser, $stp, $tmp, $x, $y);
  my @cof = (0.0, 76.18009172947146, -86.50532032941677,
	     24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
	     -0.5395239384953e-5);
  $stp = 2.5066282746310005;
    
  $x = $xx; $y = $x;
  $tmp = $x + 5.5;
  $tmp = ($x+0.5)*log($tmp)-$tmp;
  $ser = 1.000000000190015;
  foreach $j ( 1 .. 6 ) {
    $y+=1.0;
    $ser+=$cof[$j]/$y;
  }
  return $tmp + log($stp*$ser/$x);
}
##################################################################################
sub gamma { 
  return exp(&lgamma ($_[0]));
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
##############################   DA END   ########################################
##################################################################################
