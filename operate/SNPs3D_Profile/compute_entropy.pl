#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;

##################################################
## use updated psi-blast (blast+, 500 max hits) ##
## to compute entropy for all HGMD, PCSSM and   ##
## CCSM genes                                   ##
## -------------------------------------------- ##
## PSSM uses the new one:                       ##
## yini_new_snps3d_test.psimtx_250_3_extra@moult-db ##
## -------------------------------------------- 3#
## see relevant documents in:                   ##
## /moulthome/yini/work/project/scripts/self/sn ##
## ps3d_modification/psi_profile/               ##
## -------------------------------------------- ##
## Yizhou Yin, Moult lab, 17JAN2016             ##
## -------------------------------------------- ##
## Modified to integrate into the SNPs3D pipeline ##
## Yizhou, 25JAN2018                            ##
##################################################

my $refseq = $ARGV[0];
my $log = $ARGV[1];
my $dbpass = $ARGV[2];
my $working = $ARGV[3];
my $table = $ARGV[4];
my $dbuser = $ARGV[5];
my $mysqldb = $ARGV[6];
my $path;
BEGIN { $path = $ARGV[7] };
my $dbserver = $ARGV[8];
use lib "$path/operate/SNPs3D_Profile";
use PsiPsr;
use Column;

my $entropy_directory = "$working/prof_lib/entropy";
my $psi_path = "$working/prof_lib/psi/";
my $dbh = DBI->connect("DBI:mysql:$mysqldb:$dbserver", $dbuser, $dbpass, {RaiseError=>1, AutoCommit=>1}) || die "Cannot connect database: ".DBI->errstr;
my $sth0 = $dbh->prepare_cached("SELECT res FROM $table WHERE seq_ac = ? AND rnum = ?");
my $sth = $dbh->prepare_cached("UPDATE $table SET entropy = ?, relative_entropy = ?, alignment_depth = ?, relative_alignment_depth = ?, mean_entropy = ?, std_entropy = ?, Z_entropy = ?, L_outlier = ?, H_outlier = ?, accepted_hits = ?, align = ? WHERE seq_ac = ? AND rnum = ? AND res = ?");

## Computes entropy values of each sequence
system("mkdir -p $working/prof_lib/entropy");
my $parser = PsiPsr->new("$psi_path/$refseq.psi", "new4", 1, "$working/prof_lib/entropy/$log");
my $cols = Column->new($parser);
my $seq_ac = $refseq;
my $mean = $cols->{"mean"};
my $std = $cols->{"std"};
my $L_outlier = $parser->{"L_outlier"};
my $H_outlier = $parser->{"H_outlier"};
my $hits = scalar keys %{$parser->{"selected"}};
my $mean_depth = 0;

# check if query index is right #
my @sorted = sort {$a <=> $b} keys %{$parser->{"query"}};
if (($#sorted+1) ne $cols->{"col"}) {
	print "$refseq:\tparser->{query} ".($#sorted + 1)." not the same as column->{mtx} ".$cols->{"col"}."\n";
}
else {
	for (my $i = 0; $i < $cols->{"col"}; $i++) {
		$mean_depth += $cols->{"depth"}->{$i};
	}
	if ($cols->{"col"} eq 0) {
		print "$seq_ac:\treading error, 0 residue\n";
	}
	else {
		$mean_depth = $mean_depth / $cols->{"col"};
		
		## update to precomputed table
		$dbh->begin_work;	## temporarily set AutoCommit to 0
		my $status = eval {
			my $whole = "ok";	## any error residue?
			for (my $i = 0; $i <= $#sorted; $i++) {
				my $rnum = $sorted[$i];

				$sth0->execute($seq_ac, $rnum);
				my @q0 = $sth0->fetchrow();
				$sth0->finish();
				
				my $res = $q0[0];
				my $entropy = $cols->{"entropy"}->{$i};
				my $Z = $cols->{"Z"}->{$i};
				my $align = $cols->{"align"}->{$i};
				my $rl_entropy = "";
				my $depth = $cols->{"depth"}->{$i};
				my $rl_depth = "";
				if ($entropy ne "" && $mean ne "" && $mean ne 0) {
					$rl_entropy = $entropy / $mean;
				}
				if ($mean_depth ne 0) {
					$rl_depth = $depth/ $mean_depth;
				}
				if ($entropy eq "") {
					$entropy = undef;
				}
				if ($Z eq "") {
					$Z = undef;
				}
				if (defined $mean && $mean eq "") {
					$mean = undef;
				}
				if (defined $std && $std eq "") {
					$std = undef;
				}
				if ($rl_entropy eq "") {
					$rl_entropy = undef;
				}
				if ($depth eq "") {
					$depth = undef;
				}
				if ($rl_depth eq "") {
					$rl_depth = undef;
				}
					
				$sth->execute($entropy, $rl_entropy, $depth, $rl_depth, $mean, $std, $Z, $L_outlier, $H_outlier, $hits, $align, $seq_ac, $rnum, $res);
				$sth->finish();
			}
			return $whole;	## any error residue or insertion error?
		};
		
		if ($@) {
			print $@;
			$dbh->rollback();
		}
		elsif ($status eq "nah") {
			print "error residue(s) in $seq_ac\n";
		}
		elsif ($status eq "ok") {
			$dbh->commit();
		}
	}
}
$dbh->disconnect();
