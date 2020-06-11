package PsiPsr;

use strict;
use warnings;
use Exporter qw(import);
use Data::Dumper;
our @EXPORT_OK;

our $LOW_IDEN = 0.3;
our $HIGH_IDEN = 0.9;
our $COVERAGE = 0.5;
our $LINE_LENGTH = 60;

# ======================== #
# read psi-blast output    #
# filter and parse to      #
# multiple alignment       #
# with qualified sequences #
# ------------------------ #
# Yizhou Yin, Moult Lab    #
# OCT 2014                 #
# ======================== #

# constructor new #
# my $psiparser = Psipsr->new($psi_file, $outfm, $round, $logfile);
# accepted outfm: 1) old6
# variables: %{$self->{"full"}, %{$self->{"selected"}}, %{$self->{"text"}}, $self->{"log"}, $self->{"blk"}, $self->{"substr"}, $self->{"L_outlier"}, $self->{"H_outlier"}, %{$self->{"query"}}, @{$self->{"selected_key"}}, @{$self->{"list"}}, 
# methods: new; process; reorder; indexquery;
#          reader_o6, filter_o6, id_o6, length_o6, verbose_o6, findpos_o6, exception_handle,
#          reader_n4, findpos_n4
sub new {
	my $class = shift @_;
	my $self = {};
	bless($self, $class);
	
	$self->{"query"} = {};	## query sequence index
	$self->{"full"} = {};	## full alignment, including query
	$self->{"selected"} = {};	## hits to be selected, excluding query, could be null
	$self->{"selected_key"} = [];	## sorted key in selected list, excluding query, could be null
	$self->{"list"} = [];	## list of hit names, including query
	#$self->{"sum"} = {};	## statistics of selected hits
	#$self->{"aln"} = {};	## alignment of selected hits
	$self->{"text"} = {};	## alignment text row, including query
	$self->{"blk"} = 0;	## number of chunks in alignment
	$self->{"substr"} = 0;	## number of characters before alignment sequence
	#$self->{"aln_len"} = 0; ## length of one alignment line
	$self->{"L_outlier"} = 0;	## number of hits that fail the low_iden filter
	$self->{"H_outlier"} = 0;	## number of hits that fail the high_iden filter
	
	my $file = shift @_;
	my $outfm = shift @_;
	my $round = shift @_;
	my $logfile = shift @_;
		
	open (my $l, ">>$logfile") || die;
	$self->{"log"} = $l;
	
	my $signal = $self->exception_handle($file);			## handle exception when no hits were found
	if ($signal) {
		$self->process($file, $outfm, $round);
		close $l || die;
		undef $self->{"log"};
	}
	
	return $self;
}


# handle exception when no hits were found
# my $signal = $self->exception_handle($file);
sub exception_handle {
	my $self = shift @_;
	my $file = shift @_;
	my $l = $self->{"log"};

	my $flag1 = `grep \"Results from round\" $file`;
	chomp $flag1;
	my $flag2 = `grep \"No hits found\" $file`;
	chomp $flag2;
	if ($flag1 eq "" || $flag2 ne "") {
		print $l "No hits were found in $file\n";
		close $l || die;
		$self = {};
		return 0;
	}
	else {
		return 1;
	}
}


# process call and distribute reader #
# $self->process($file, $outfm, $round);
sub process {
	my $self = shift @_;
	if ($_[1] eq 'old6') {	## old psi-blast -m 6, -b = -v
		$self->reader_o6($_[0], $_[2]);
		$self->filter_o6();
		$self->indexequery();
		
		my @key = keys %{$self->{"selected"}};
		if ($#key >= 0) {
			@{$self->{"selected_key"}} = $self->reorder($self->{"list"}, \@key);
		}
	}
	
	if ($_[1] eq 'new4') {	## new psi-blast (blast+) -outfmt 4, -num_descriptions = -num_alignments, -max_hsps 1
		$self->reader_n4($_[0], $_[2]);
		$self->filter_o6();	## use the same method as old6
		$self->indexequery();
		
		my @key = keys %{$self->{"selected"}};
		if ($#key >= 0) {
			@{$self->{"selected_key"}} = $self->reorder($self->{"list"}, \@key);
		}
	}
}


# reorder a list based on a given list #
# @target must be a subset of @template #
# @new = $self->reorder(\@template, \@target);
sub reorder {
	my $self = shift @_;
	my $tpref = shift @_;
	my $tgref = shift @_;
	
	my %order;
	for (my $i = 0; $i <= $#{$tpref}; $i++) {
		$order{$tpref->[$i]} = $i;
	}
	
	my %tmp;
	for (my $i = 0; $i <= $#{$tgref}; $i++) {
		$tmp{$tgref->[$i]} = $order{$tgref->[$i]};
	}
	
	my @sorted = sort {$tmp{$a} <=> $tmp{$b}} keys %tmp;
	return @sorted;
}


# index query sequence
# $self->indexquery();
sub indexequery {
	my $self = shift @_;
	my $start; my $end;
	if ($self->{"text"}->{'QUERY'}->{0}->[0] =~ m/^QUERY\s+(\d+)/) {
		$start = $1;
	}
	if ($self->{"text"}->{'QUERY'}->{0}->[-1] =~ m/\s+(\d+)$/) {
		$end = $1;
	}
	my $n = $start-1;
	foreach my $pos (split(//, $self->{"full"}->{'QUERY'}->{0})) {
		if ($pos =~ m/[A-Z]/) {
			$n++;
			$self->{"query"}->{$n} = $pos;
		}
	}
	my $lg = $self->{"log"};
	if ($n ne $end) {
		print $lg "start:$start\tend:$end\tn:$n\n";
		$self->{"query"} = {};
	}
}


# read legacy psi-blast -m 6 output #
# $self->reader_o6($file, $round);
sub reader_o6 {
	my $self = shift @_;
	my $file = shift @_;
	my $rd = shift @_;
#	my $blk = 0;
	
	# mark the reading start and end points
	my $on = 0;
	my $go = -1;
	my $stop = -1;
	my $count = 0;
	
	open (my $f, $file) || die;
	while (<$f>) {
		$count++;
		if ($_ =~ m/^Results from round\s+(\d+)/) {
			if ($1 eq $rd) {
				$on = 1;
			}
			else {
				$on = 0;
			}
		}
		if ($on && $_ =~ m/^QUERY/ && $go < 0) {
			$go = $count - 1;
			chomp $_;
			$self->findpos_o6($_);
		}
		if ($on && $_ =~ m/^Searching|^  Database:/ && $stop < 0) {
			$stop = $count - 1;
		}
		if ($go > 0 && $stop > 0) {
			last;
		}
	}
	close $f || die;
	
	# truncate file
	my $tail = $stop - $go + 1;
	system("head -n $stop $file | tail -n $tail > psi_temp");
	
	# read temp file
	open (my $t, "psi_temp") || die;
	{
		local $/ = "";
		my %tmp;
		my $ct = 0;
		while (<$t>) {
			$ct++;
			my @array = split(/\n/, $_);
			$self->{"blk"}++;
			for (my $i = 0; $i <= $#array; $i++) {
				my @line = split(/\s+/, $array[$i]);
				my $content = substr $array[$i], $self->{"substr"};
				my @content = split(/\s+/, $content);
#print substr $array[$i], $self->{"substr"}, $self->{"aln_len"};
#print "\n";
#exit -1;
				$self->{"full"}->{$line[0]}->{$i} .= $content[0];
				push(@{$self->{"text"}->{$line[0]}->{$i}}, $array[$i]);

				if ($ct < 2) {
					if (! defined $tmp{$line[0]}) {
						push(@{$self->{"list"}}, $line[0]);
						$tmp{$line[0]} = 1;
					}
				}
			}
		}
	}
	close $t || die;
}


# filter old psi-blast outfmt 6 #
# Algorithm: 1) take the top one of the duplicated hits as representative
#            2) remove sequences that have less than 30% sequence identity with "QUERY"
#            3) Start from the "QUERY", check every possible pair and remove the second one if sequence identity between the pair is more than 90%.
# $self->filter_o6();
sub filter_o6 {
	my $self = shift @_;
	my $l = $self->{"log"};
	my @close = ('QUERY');
	my %closeid;
	$closeid{'QUERY'} = 1;
	foreach my $key (keys %{$self->{"full"}}) {
		if ($key ne 'QUERY') {

## !!! the algorithm below rank the duplicat hits of the same refseq by its sequence identity to the query sequence, this make more sense to me but is not the same as in the previous codes, so will not be used here --- Yizhou ##
=begin			
			
			my $rep = 0;
			my $max = 0;
			my $largest = 0;	
			my @largest = (0,0,0);
			foreach my $ky (keys %{$self->{"full"}->{$key}}) {
				my ($idp, $same, $length) = $self->id_o6('QUERY', 0, $key, $ky);
				if ($idp > $largest[0]) {
					@largest = ($idp, $same, $length);
					
				}
				if ($idp <= $HIGH_IDEN && $idp >= $LOW_IDEN && $idp > $max) {
					$max = $idp;
					$rep = $ky;
				}
			}
			if ($rep ne 'none') {
				@{$self->{"selected"}->{$key}} = ($rep, $max);
			}
			else {
				print $l "No representative for $key with largest identity at $largest[0]\t$largest[1]\t$largest[2]\n";
			}
=cut
#######################################################################################
			
			## this algorithm only take the top of the duplicate hits (presumably ranked by e-value as the only representative, and if it doesn't pass the filters, the algorithm won't try its other duplicates ##
			my @sorted = sort {$a <=> $b} keys %{$self->{"full"}->{$key}};
			my ($idp, $same, $length) = $self->id_o6('QUERY', 0, $key, $sorted[0]);

#if (! defined $same) {								## DEBUG
#print "key\t$key\nsorted[0]\t$sorted[0]\nidp\t$idp\nlength\t$length\n\n";	## DEBUG
#for (my $i = 0; $i < $self->{"blk"}; $i ++) {					## DEBUG
#print $self->{"text"}->{"QUERY"}->{0}->[$i]."\n";				## DEBUG
#print $self->{"text"}->{"$key"}->{$sorted[0]}->[$i]."\n";			## DEBUG
#print "\n\n";									## DEBUG
#}										## DEBUG
#exit -1;									## DEBUG
#}										## DEBUG

			if ($idp >= $LOW_IDEN) {
				push(@close, $key);
				$closeid{$key} = $idp;
			}
			else {
				$self->{"L_outlier"}++;
				print $l "Representative $sorted[0] for $key failed 30% filter with sequence identity $idp\t--\t$same\t$length\n";
			}
			
		}
	}
	
	## make a cluster excluding one of two sequences if sequence identity larger than 90%
	## really don't like the way this piece of code was organized, but have to stay with the original code ##
	my %ignore;
#print scalar @close; print "\n"; ## DEBUG
#print scalar @{$self->{"list"}}; print "\n"; ## DEBUG
#print Dumper(@close);			## DEBUG
#print "\n\n";				## DEBUG
#print Dumper(@{$self->{"list"}});	## DEBUG
	my @sorted_close = $self->reorder($self->{"list"}, \@close);
	foreach my $cl (@sorted_close) {
		my $start = 0;
		if (exists $ignore{$cl}) {
			next;
		}
		foreach my $check (@{$self->{"list"}}) {
			if ($cl eq $check) {
				$start = 1;
				next;
			}
			if (exists $ignore{$check}) {
				next;
			}
			elsif (not $start) {
				next;
			}
			my @sortcl = sort {$a <=> $b} keys %{$self->{"full"}->{$cl}};
			my @sortch = sort {$a <=> $b} keys %{$self->{"full"}->{$check}};
			my ($id, $same, $length) = $self->id_o6($cl, $sortcl[0], $check, $sortch[0]);
#print $id."\n"; ## DEBUG
			if ($id > $HIGH_IDEN) {
				$ignore{$check} = 1;
				print $l "Representative $sortch[0] for $check failed 90% filter with sequence identity $id\t--\t$same\t$length to $cl\_$sortcl[0]\n";
			}
		}
	}
#print scalar keys %ignore; print "\n"; ## DEBUG
	foreach my $cl (@sorted_close ) {
		my @sort = sort {$a <=> $b} keys %{$self->{"full"}->{$cl}};
		if (! exists $ignore{$cl} && $cl ne 'QUERY') {
			@{$self->{"selected"}->{$cl}} = ($sort[0], $closeid{$cl});
		}
		else {
			$self->{"H_outlier"}++;
		}
	}
}



# calculate sequence identity on old blast output -m 6
# equation: sequence identity = identical positions / (aligned positions(without external gap from either of the two sequences) + internal gaps)
# !! this equation is pretty unique. The more common ones are:
# !! 1) NCBI blast: sequence identity = identical positions / (aligned positions(blunt end) + internal gaps)
# !! 2) sequence identity = identical positions / aligned positions (without any gaps)
# my $iden = $self->id_o6($A_acc, $A_rep, $B_acc, $B_rep);
sub id_o6 {
	my $self = shift @_;
	my ($q_acc, $q_rep, $t_acc, $t_rep) = @_;
	my $l = $self->{"log"};
	
	my @q = split(//, $self->{"full"}->{$q_acc}->{$q_rep});
	my @t = split(//, $self->{"full"}->{$t_acc}->{$t_rep});
#	my $lenq = $self->length_o6($q_acc, $q_rep);
	if ($#q ne $#t) {
		print $l "error in blunt-end alignment $q_acc\_$q_rep and $t_acc\_$t_rep!\n";
		print $l "$#q - $#t\n";
		print $l $self->{"full"}->{$q_acc}->{$q_rep}."\n";
		print $l $self->{"full"}->{$t_acc}->{$t_rep}."\n";
		return 0;
	}
	else {
		############# chop off the "one-strand" ends ##############
		my $start = 0;
		my $stop = $#q;
		my $mark = 1;
		for (my $i = 0; $i <= $#q; $i++) {
			if ($mark && ($q[$i] eq '-' || $t[$i] eq '-')) {
				$start++;
			}
			elsif ($mark && ($q[$i] ne '-' && $t[$i] ne '-')) {
				$mark = 0;
				last;
			}
			else {
				last;
			}
		}
		$mark = 1;
		for (my $j = $#q; $j >= 0; $j--) {
			if ($mark && ($q[$j] eq '-' || $t[$j] eq '-')) {
				$stop--;
			}
			elsif ($mark && ($q[$j] ne '-' && $t[$j] ne '-')) {
				$mark = 0;
				last;
			}
			else {
				last;
			}
		}	
		###########################################################
		
		my $same = 0;
		my $lenq = $stop - $start + 1;
		for (my $i = $start; $i <= $stop; $i++) {
			if ($q[$i] ne '-' && $q[$i] eq $t[$i]) {
				$same++;
			}
			if ($q[$i] eq '-' && $t[$i] eq '-') {
				$lenq--;
			}
		}
		return (($same/$lenq), $same, $lenq);
		#return ((int(($same/$lenq)*100)/100), $same, $lenq);	## really don't like the this implementation, should use the original floating number ($name/$lenq), but must keep this way if need to be consistent with the previous code -- Yizhou
	}
}


# calculate sequence length for old blast output format -m 6
# my $length = $self->length_o6($acc, $rep);
sub length_o6 {
	my $self = shift @_;
	my ($acc, $rep) = @_;
	my @str = split(//, $self->{"full"}->{$acc}->{$rep});
	my $len = 0;
	for (my $i = 0; $i <= $#str; $i++) {
		if ($str[$i] =~ m/[A-Z]/) {
			$len++;
		}
	}
	return $len;
}


# output alignment and statistics of the stored selected hits from old blast -m 6
# $self->verbose_o6($filename);
sub verbose_o6 {
	my $self = shift @_;
	my $filename = shift @_;
	
	my @sorted_key = @{$self->{"selected_key"}};
	
	open (my $aln, ">$filename.aln") || die;
	open (my $sum, ">$filename.sum") || die;
	print $sum "QUERY\t0\t1\n";
	for (my $i = 0; $i <= $#{$self->{"text"}->{'QUERY'}->{0}}; $i++) {
		print $aln $self->{"text"}->{'QUERY'}->{0}->[$i]."\n";
		foreach my $key (@sorted_key) {
			if ($i eq 0) {
				print $sum $key."\t".$self->{"selected"}->{$key}->[0]."\t".$self->{"selected"}->{$key}->[1]."\n";
			}
			print $aln $self->{"text"}->{$key}->{$self->{"selected"}->{$key}->[0]}->[$i]."\n";
		}
		print $aln "\n";
	}
	close $sum || die;
	close $aln || die;
}


# find where does the alignment sequence start in old blast -m 6
# $self->findpos_o6($string);
sub findpos_o6 {
	my $self = shift @_;
	my $string = shift @_;
	if ($string =~ m/^QUERY(\s*\d*\s*)([A-Z]|-)/) {
		$self->{"substr"} = length($1) + 5;
	}
#	if ($string =~ m/([A-Z]|-)(\s*\d*)$/) {
#		if (defined $2) {
#			$self->{"aln_len"} = length($string) - $self->{"substr"} - length($2);
#		}
#		else {
#			$self->{"aln_len"} = length($string) - $self->{"substr"};
#		}
#	}
}


# read new psi-blast -outfmt 4 output and with -max_hsps 1 (one entry per hit) and -line_length 60 (default line length of alignment block)#
# $self->reader_n4($file, $round);
sub reader_n4 {
	my $self = shift @_;
	my $file = shift @_;
	my $rd = shift @_;
	my $l = $self->{"log"};
	
	# mark the reading start and end points, also the complete gi list
	my @gi = ("QUERY");
	my %gi = ("QUERY" => 0);
	my %rank = ("QUERY" => 0);
	my $on = 0;
	my $go = -1;
	my $stop = -1;
	my $count = 0;
	my $last_substr;
	my $qy;
	my $nqy = 0;
	
	open (my $f, $file) || die;
	while (<$f>) {
		$count++;
		if ($_ =~ m/^Results from round\s+(\d+)/) {
			if ($1 eq $rd) {
				$on = 1;
			}
			else {
				$on = 0;
			}
		}
		if ($on && $_ =~ m/^Query_1/) {
			$qy = $_;
			$nqy++;
		}
		if ($on && $_ =~ m/^gi\|(\d+)\|/) {
			if (! defined $gi{$1}) {
				$gi{$1} = 0;
				$rank{$1} = scalar @gi;
				push(@gi, $1);
			}
		}
		if ($on && $_ =~ m/^Query_1/ && $go < 0) {
			$go = $count - 1;
			chomp $_;
			$self->findpos_n4($_);
		}
		if ($on && $_ =~ m/^Lambda/ && $stop < 0) {
			$stop = $count - 3;
		}
		if ($go > 0 && $stop > 0) {
			last;
		}
	}
	close $f || die;
	@{$self->{"list"}} = @gi;	## copy whole hit-list
	chomp $qy;
	$last_substr = $self->findpos_n4($qy, 1);	## find "substr" for the last block
		
	# truncate file #
	my $tail = $stop - $go + 1;
	system("head -n $stop $file | tail -n $tail > psi_temp");
	
	# read temp file
	my $empty;
	for (my $i = 0; $i < $LINE_LENGTH; $i++) {
		$empty .= "-";
	}
	my $empty_last;
	for (my $i = 0; $i < ($last_substr-$self->{"substr"}); $i++) {
		$empty_last .= "-";
	}
	
	open (my $t, "psi_temp") || die;
	{
		local $/ = "";
		my %tmp;
		my $ct = 0;
		while (<$t>) {
			$ct++;
			my @array = split(/\n/, $_);
			$self->{"blk"}++;
			for (my $i = 0; $i <= $#array; $i++) {
				my @line = split(/\s+/, $array[$i]);
				if ($line[0] =~ m/Query_1/) {
					$line[0] = "QUERY";
				}
				my $head = substr $array[$i], 0, $self->{"substr"};
				$head =~ s/Query_1/QUERY  /g;
				my $content;
				my $tail;
				if ($ct eq $nqy) {	## adjustment for the last block
					$content = substr $array[$i], $self->{"substr"}, ($last_substr - $self->{"substr"});
					$tail = substr $array[$i], $last_substr;
				}
				else {
					$content = substr $array[$i], $self->{"substr"}, $LINE_LENGTH;
					$tail = substr $array[$i], ($self->{"substr"}+$LINE_LENGTH);
				}
					
				$content =~ s/\s/-/g;
				my $newline = $head.$content.$tail;
				
				if (! defined $gi{$line[0]}) {
					print $l "missing $line[0] in gi-list while reading\n";
				}
				else {	## -max_hsps 1 guarantee the only entry of each hit is the one with the highest e-value ##
					if ($gi{$line[0]} < $self->{"blk"}) { ## assume the order in gi-list is the same as in alignment ##
						$self->{"full"}->{$line[0]}->{$rank{$line[0]}} .= $content;
						push(@{$self->{"text"}->{$line[0]}->{$rank{$line[0]}}}, $newline);
						$gi{$line[0]}++;
					}
					
				}
#if ($i eq 0 && $line[0] eq "QUERY") {		## DEBUG
#print "array[0]\t$array[$i]\n";		## DEBUG
#print "line[0]\t$line[0]\n";			## DEBUG
#print "substr\t".$self->{"substr"}."\n";	## DEBUG
#print "content\t$content\n";			## DEBUG
#print "newline\t$newline\n";			## DEBUG
#print "text[0]\t".$self->{"text"}->{$line[0]}->{$rank{$line[0]}}."\n";	## DEBUG
#}						## DEBUG
			}
			foreach my $hit (@gi) {
				if ($gi{$hit} < $self->{"blk"}) {
					if ($ct eq $nqy) {
						$self->{"full"}->{$hit}->{$rank{$hit}} .= $empty_last;
					}
					else {
						$self->{"full"}->{$hit}->{$rank{$hit}} .= $empty;
					}
					my $newline = sprintf '%*s', (-1 * $self->{"substr"}), $hit;
					if ($ct eq $nqy) {
						$newline .= $empty_last."  ";
					}
					else {
						$newline .= $empty."  ";
					}
					push(@{$self->{"text"}->{$hit}->{$rank{$hit}}}, $newline);
					$gi{$hit}++;
				}
			}
		}
	}
	close $t || die;
}


# fine where does the alignment sequence start in new blast -outfmt 4
# $self->findpos_n4($string);
# my $sub = $self->findpos_n4($string, 1);
sub findpos_n4 {
	my $self = shift @_;
	my $string = shift @_;
	if ($string =~ m/^Query_1(\s*\d*\s*)([A-Z]|-)/) {
		$self->{"substr"} = length($1) + 7;
	}
	if ($#_ >= 0 && $string =~ m/^(.*)  \d+$/) {
		return length($1);
	}
}


1;
