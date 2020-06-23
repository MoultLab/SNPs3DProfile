#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use Getopt::Long;

#### take in a list of mutation, two possibilities:
#### 1) When mutation is precomputed, compile feature, run SVM, and output score
#### 2) When mutation is not precomputed, do PSI-BLAST, compute and update table, run SVM, and output score



if ($#ARGV < 0) {
	print STDERR "1) $0 options\n";
	print STDERR "--dbuser (mysql username)\n--dbpass (mysql password)\n--adpass (mysql moult-db admin password, currently not in use)\n--log (logfile)\n--workpath (working path)\n--list (full path and name of the mutation list, no header line)\n--run_psi_blast (boolean 1 or 0, run psi-blast or skip, default 0)\n--model (SVM model, defaut HUWL2)\n--output (result filename, defaut pred_HWUL2.txt)\n--table (name of the entropy table, default SNPs3D_Profile_2017_precomputed_psimtx_500_1_0_JAN2018)\n--blastbin (path to the psiblast, default \"/moulthome/shared/apps/ncbi-blast-2.2.30+/bin/\")\n--blastdb (path to the blast database, default \"/moulthome/yini/share/BLASTDB/nr/nr\")\n--mysqldb (the mysql database to store precomputed features, default SNPs3D_2017)\n";
	print STDERR "\n2) $0 take in a tab-delimited list of mutation\n";
	print STDERR "format: type (default: 0, e.g. disease, benign ...), refseq, mutation (e.g. A235G)\n";
	print STDERR "\n3) $0 output a list of mutations with SNPs3D_Profile score\n";
	exit -1;
}

my $model = "HWUL2";
my $dbserver;
my $mysqldb;
my $table = "SNPs3D_Profile_2017_precomputed_psimtx_500_1_0_JAN2018";
my $dbuser;
my $dbpass;
my $blastbin;
my $blastdb;
my $rpb = 1;
my $path; 
my $list;
my $output = "pred_HWUL2.txt";
my $log = "log.profile_pipeline";

## Get options
GetOptions('dbserver=s' => \$dbserver, 'mysqldb=s' => \$mysqldb,'table=s' => \$table, 'dbuser=s' => \$dbuser, 'dbpass=s' => \$dbpass,'blastbin=s' => \$blastbin,'blastdb=s' => \$blastdb, 'run_psi_blast=i' => \$rpb, 'path=s' => \$path, 'list=s' => \$list, 'output=s' => \$output, 'log=s' => \$log); 

my $bin = "$path/operate/SNPs3D_Profile";
my $modelpath = "$path/model/"; 
my $working = "$path/bin";	
	
if (not defined $blastbin)
{
	$blastbin = "$path/third_party_apps/ncbi-blast-2.2.30+/bin";
}

if (not defined $blastdb)
{
	$blastdb = "$path/profile_data/BLASTDB/nr";
}
print "BLASTDB is $blastdb \n BLASTbin is $blastbin\n";
open (my $l, ">$working/$log") || die "Cannot open $working/$log: $!";

## Read list of proteins
print $l "READING LIST\n";
my %list;
my %missing;
open (my $lf, $list) || die "Cannot open $list: $!";		
while (<$lf>) {
	chomp;								
	my @line = split(/\s+/, $_);							
	my @mut = split(//, $line[2]);
	my $wt = shift @mut;
	my $mt = pop @mut;
	my $rnum = join("", @mut);
	my $key = "$line[1]\_$line[2]";  
	$list{$key}->{type} = $line[0];
	$list{$key}->{refseq} = $line[1];
	$list{$key}->{wt} = $wt;
	$list{$key}->{mt} = $mt;
	$list{$key}->{rnum} = $rnum;
	$missing{$line[1]} = 1;	
}
close $lf || die;
print "Created list of mutations\n";

## Create list of non-computed proteins 
print $l "\nCHECKING NON-PRECOMPUTED PROTEINS\n";
my $dbh = DBI->connect("DBI:mysql:$mysqldb:$dbserver", $dbuser, $dbpass, {RaiseError=>1, AutoCommit=>1});
my $sth = $dbh->prepare_cached("SELECT COUNT(*) FROM $table WHERE seq_ac = ?");
print "List of non-computed proteins:\n"; 
foreach my $p (keys %missing) {
	$sth->execute($p);
	my @q = $sth->fetchrow();	
	$sth->finish();				
	if (defined $q[0] && $q[0] > 0) {
		delete $missing{$p};
	}
	else {
		print $l "need to run psi-blast for:\t$p\n";
	}
}
foreach my $p (keys %missing) {
	print "$p\n";
}

## Runs psi-blast and computes entropy on missing mutations 
if ($rpb ne 0) {
	print $l "RUNNING PSI-BLAST AND COMPUTING ENTROPY FOR NEW PROTEINS\n";
	foreach my $p (keys %missing) {
		my $cmd1 = "perl $bin/run_psiblast.pl $dbuser $dbpass $p $working $table $mysqldb $blastbin $blastdb $dbserver";
		if (defined $blastbin && $blastbin ne "") {
			$cmd1 .= " $blastbin"; 
		}
		if (defined $blastdb && $blastdb ne "") {
			$cmd1 .= " $blastdb";
		} 
		print"Running psiblast on $p\n";
		system($cmd1);
		system("perl $bin/compute_entropy.pl $p log.psipsr_$p $dbpass $working $table $dbuser $mysqldb $path $dbserver | tee $working/log.compute_entropy");
	}
	$rpb = 0;
}


## Use PSSM and entropy information compile and analyze features with SVM    
if ($rpb eq 0) {
	print $l "COMPILING FEATURES\n";
	system("perl $bin/compile_features.pl $dbuser $dbpass $working $list $mysqldb $table $dbserver");	
	print $l "SVM_CLASSIFY\n";
	system("$path/third_party_apps/svm_light/svm_classify $working/total.data $modelpath/$model $working/pred_$model\_raw > log.svm_classsify_$model");
	print $l "WRITING RESULTS\n";
	system("perl $bin/glue.pl $working/total.index pred_$model\_raw $output");
}
print "Completed analysis using SVM\n";

## Clean up
$dbh->disconnect();
close $l || die "Cannot close $working/$log: $!";
