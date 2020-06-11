#!/usr/bin/perl -w
use strict;
use warnings;

### 26JAN2018, YY ###
### adding mutation names to SNPs3D Profile results ###

my $ttl = $ARGV[0];
my $pred = $ARGV[1];
my $output = $ARGV[2];

open (my $t, $ttl) || die "Cannot open $ttl: $!";
open (my $p, $pred) || die "Cannot open $pred: $!";
open (my $o, ">$output") || die "Cannot open $output: $!";
<$t>;
print $o "type\trefseq\tmutation\tSNPs3D_Profile_score\n";
while (<$t>) {
	my $tline = $_;
	chomp $tline;
	my $pline = <$p>;
	chomp $pline;
	
	my @tline = split(/\s+/, $tline);
	print $o "$tline[0]\t$tline[1]\t$tline[2]\t$pline\n";
}
close $o || die "Cannot close $output: $!";
close $p || die "Cannot close $pred: $!";
close $t || die "Cannot close $ttl: $!";
