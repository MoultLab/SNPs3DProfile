#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;

#### Yizhou Yin, 25JAN2018 ####
#### take in a list of mutation and compile features from precomputed table
if ($#ARGV < 0) {
	print STDERR "$0 <mysql username> <mysql password> <working path> <full path and name of the list>\n";
	print STDERR "take in a tab-delimited list of mutation\n";
	print STDERR "format: type (default: 0), refseq, mutation (e.g. A235G)\n";
	exit -1;
}

my $dbpass = $ARGV[1];					
my $dbuser = $ARGV[0];
my $working = $ARGV[2];
my $rawf = $ARGV[3];
my $mysqldb = $ARGV[4];
my $table = $ARGV[5]; 
my $dbserver = $ARGV[6];
system("mkdir -p $working");
my $output1 = "$working/total.data";
my $output3 = "$working/total.index";
my $log = "$working/log.compile_error";

my $dbh = DBI->connect("DBI:mysql:$mysqldb:$dbserver", $dbuser, $dbpass, {RaiseError=>1, AutoCommit=>1}) || die "Cannot connect database ".DBI->errstr;
my $sth = $dbh->prepare_cached("SELECT * FROM $table WHERE seq_ac = ? AND rnum = ? AND res = ?");

open (my $r, $rawf) || die;
open (my $o1, ">$output1") || die;
open (my $o3, ">$output3") || die;
open (my $l, ">$log") || die;
#<$r>t;
print $o3 "dataset_type Refseq mutation PSSM_score entropy_mean entropy_std entropy_Z entropy\n";
while (<$r>) {
	chomp;
	my @line = split(/\t/, $_);
	my $type = $line[0];
	my $refseq = $line[1];
	my @mut = split(//, $line[2]);
	my $wt = shift @mut;
	my $mt = pop @mut;
	my $rnum = join("", @mut);
	$sth->execute($refseq, $rnum, $wt) || die "Cannot execute: ".$sth->errstr;
	my $q0 = $sth->fetchrow_hashref();
	if (! defined $q0->{'entropy'} || ! defined $q0->{$mt} || ! defined $q0->{'std_entropy'} || ! defined $q0->{'mean_entropy'} || ! defined $q0->{'Z_entropy'}) {
		print $l "missing_mutation:\t$_\n";
		print "Missing information for mutation $wt$rnum$mt in $refseq\n";
	}
	else {	
		
		if (! defined $type || $type eq "") {
			$type = 0;
		}
		my $pssm = $q0->{$mt};
		my $mean = $q0->{'mean_entropy'};
		my $std = $q0->{'std_entropy'};
		my $Z = $q0->{'Z_entropy'};
		my $entropy = $q0->{'entropy'};
		
		print $o1 "$type 1:$pssm 2:$mean 3:$std 4:$Z 5:$entropy\n";
		print $o3 "$type $refseq ".$wt.$rnum.$mt." 1:$pssm 2:$mean 3:$std 4:$Z 5:$entropy\n";
		print "Added the compiled features of mutation $wt$rnum$mt from $refseq into file \n";
	}
	$sth->finish();
}
close $l || die;
close $o3 || die;
close $o1 || die;
close $r || die;

$dbh->disconnect();
