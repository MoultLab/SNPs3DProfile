#!/usr/bin/perl -w
use strict;
use warnings;
use LWP::Simple;
use Data::Dumper;
use DBI;

#### Yizhou Yin, 25JAN2018 ####
sub read_mtx;
my $dbuser = $ARGV[0];
my $dbpass = $ARGV[1];
my $refseq = $ARGV[2];
my $working = $ARGV[3];
my $table = $ARGV[4];
my $mysqldb = $ARGV[5];
my $blastbin = $ARGV[6];
my $blastdb = $ARGV[7];
my $dbserver = $ARGV[8]; 

#### get fasta file using NCBI eutil ####
#########################################
my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
my $url = $base."esearch.fcgi?db=protein&term=$refseq&usehistory=y"; print "$url\n"; 
my $out1 = get($url);
my @out1 = split(/\n/, $out1);
my $WebEnv;
my $QueryKey;
if ($out1[2] =~ m/\<WebEnv\>(.*)\<\/WebEnv\>/) {
	$WebEnv = $1;
}
if ($out1[2] =~ m/\<QueryKey\>(.*)\<\/QueryKey\>/) {
	$QueryKey = $1;
}

$url = $base."efetch.fcgi?db=protein&WebEnv=$WebEnv&query_key=$QueryKey&rettype=fasta&retmode=text";
my $out2 = get($url);
my @out2 = split(/\n/, $out2);
chomp $out2;
chomp $out2; 
system("mkdir -p $working/prof_lib/fasta/");
system("echo \'$out2\' > $working/prof_lib/fasta/$refseq");

##### run psi-blast #####
system("mkdir -p $working/prof_lib/psi/");
system("mkdir -p $working/prof_lib/mtx/");
system("$blastbin/psiblast -query $working/prof_lib/fasta/$refseq -db $blastdb -out $working/prof_lib/psi/$refseq.psi -evalue 0.001 -outfmt 4 -show_gis -num_descriptions 500 -num_alignments 500 -num_iterations 3 -max_hsps 1 -num_threads 12 -out_ascii_pssm $working/prof_lib/mtx/$refseq.mtx > $working/prof_lib/psi/log.psiblast_$refseq");
## system("$blastbin/psiblast -query $working/prof_lib/fasta/$refseq -db $blastdb -out $working/prof_lib/psi/$refseq.psi -evalue 0.001 -outfmt 4 -show_gis -num_descriptions 500 -num_alignments 500 -num_iterations 3 -max_hsps 1 -num_threads 12 -out_ascii_pssm $working/prof_lib/mtx/$refseq.mtx > $working/prof_lib/psi/log.psiblast_$refseq");

##### update the PSSM numbers from the mtx file #####
my $dbh = DBI->connect("DBI:mysql:$mysqldb:$dbserver", $dbuser, $dbpass, {RaiseError=>1, AutoCommit=>1});
my $sth = $dbh->prepare_cached("INSERT IGNORE INTO $table (seq_ac, rnum, res, A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)") || die $dbh->errstr;

## in case you want to store new entropy in a new table
$dbh->do("CREATE TABLE IF NOT EXISTS $table (
		seq_ac VARCHAR(20) NOT NULL,
		rnum MEDIUMINT NOT NULL,
		res CHAR(1) NOT NULL,
		entropy FLOAT,
		relative_entropy FLOAT,
		alignment_depth SMALLINT,
		relative_alignment_depth FLOAT,
		mean_entropy FLOAT,
		std_entropy FLOAT,
		Z_entropy FLOAT,
		A TINYINT NOT NULL,
		R TINYINT NOT NULL,
		N TINYINT NOT NULL,
		D TINYINT NOT NULL,
		C TINYINT NOT NULL,
		Q TINYINT NOT NULL,
		E TINYINT NOT NULL,
		G TINYINT NOT NULL,
		H TINYINT NOT NULL,
		I TINYINT NOT NULL,
		L TINYINT NOT NULL,
		K TINYINT NOT NULL,
		M TINYINT NOT NULL,
		F TINYINT NOT NULL,
		P TINYINT NOT NULL,
		S TINYINT NOT NULL,
		T TINYINT NOT NULL,
		W TINYINT NOT NULL,
		Y TINYINT NOT NULL,
		V TINYINT NOT NULL,
		L_outlier SMALLINT,
		H_outlier SMALLINT,
		accepted_hits SMALLINT,
		align VARCHAR(510),
		DISOPRED2_pred TINYINT,
		DISOPRED2_conf TINYINT,
		DISOPRED2_idr_size SMALLINT,
		DISOPRED3_pred_2017 TINYINT,
		DISOPRED3_conf_2017 FLOAT,
		DISOPRED3_idr_size_2017 SMALLINT,
		DISOPRED3_pbdat_pred_2017 TINYINT,
		DISOPRED3_pbdat_conf_2017 FLOAT,
		DISOPRED3_pbdat_size_2017 SMALLINT,
		PRIMARY KEY (seq_ac, rnum, res))");

my %pssm = read_mtx("$working/prof_lib/mtx/$refseq.mtx");
my @sorted = sort {$a <=> $b} keys %pssm;
$dbh->begin_work;
eval {
	foreach my $rnum (@sorted) {
		$sth->execute($refseq, $rnum, $pssm{$rnum}->{WT}, @{$pssm{$rnum}->{PSSM}});
		$sth->finish();
	}
};
if ($@) {
	$dbh->rollback();
}
else {
	$dbh->commit();
}
$dbh->disconnect();


##### subroutein #####
sub read_mtx {
	my $file = shift @_;
	my %hash;
	open (my $f, $file) || die;
	<$f>; <$f>; <$f>;
	while (<$f>) {
		if ($_ !~ m/^[A-Z]/ && $_ !~ m/^$/ && $_ !~ m/Lambda/) {
			my @line = split(/\s+/, $_);
			if ($line[0] ne "") {
				unshift (@line, "");
			}
			$hash{$line[1]}->{WT} = $line[2];
			for (my $i = 3; $i < 23; $i++) {
				push(@{$hash{$line[1]}->{PSSM}}, $line[$i]);
			}
		}
	}
	close $f || die;
	return %hash;
}

