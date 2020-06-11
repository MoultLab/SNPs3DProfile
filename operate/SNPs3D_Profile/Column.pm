package Column;

use strict;
use warnings;
use Exporter qw(import);
use Data::Dumper;
use Reference;
our @EXPORT_OK;

our $COVERAGE = 0.5;

# ====================== #
# Take in an alignment   #
# and an indexed query   #
# then calculate entropy #
# for each column and    #
# other statistics       #
#                        #
# Work with PsiPsr.pm    #
# ---------------------- #
# Yizhou Yin, Moult lab  #
# OCT 2014               #
# ====================== #

# constructor new #
# my $columns = Column->new($psiparser);
# variables: %{$self->{"mtx"}, %{$self-{"query"}}, $self->{"std"}, $self->{"mean"}, %{$self->{"entropy"}}, %{$self->{"Z"}}
# methods: initialize, transpose, calculate, composition, empty, calculate_entropy, text, Verbose

sub new {
	my $class = shift @_;
	my $self = {};
	bless ($self, $class);
	
	my $ref = shift @_;
	
	$self->{"mtx"} = {};	## column matrix organized by row
	$self->{"MTX"} = {};	## column matrix capticalized and organized by column
	#$self->{"query"} = {%{$ref->{"query"}}};## copy $psiparser->{"query"}
	$self->{"std"} = -1;	## entropy std of the query, could be null
	$self->{"mean"} = -1;	## entropy mean of the query, could be null
	$self->{"entropy"} = {};	## entropy of every position in query, could be null
	$self->{"Z"} = {};	## entropy Z-score of every position in query, could be null
	$self->{"row"} = 0;	## total number of rows in column matrix, is at least 1
	$self->{"col"} = 0;	## total number of columns in column matrix
	$self->{"align"} = {};	## text alignment of each columns
	$self->{"depth"} = {};	## alignment depth of each columns, is at least 1
	$self->{"gaps"} = {};	## gaps number of each columns, is at least 0
	
	$self->initialize($ref);	## fill in alignment matrix
	$self->transpose();	## organize column matrix by column and capitalize
	$self->calculate();	## calculate entropy and related numbers
	$self->text();		## store alignment for each columns
	
	return $self;
}


# fill in alignment matrix
# $self->initialize($ref);
sub initialize {
	my $self = shift @_;
	my $ref = shift @_;
	my @idx;
	
	my @array = split(//, $ref->{"full"}->{'QUERY'}->{0});
	for (my $j = 0; $j <= $#array; $j++) {
		if ($array[$j] ne "-") {
			push(@idx, $j);
			$self->{"mtx"}->{0}->{$#idx} = $array[$j];
		}
	}
	$self->{"col"} = scalar keys %{$self->{"mtx"}->{0}};
	$self->{"row"}++;
	
	for (my $i = 0; $i <= $#{$ref->{"selected_key"}}; $i++) {
		@array = split(//, $ref->{"full"}->{$ref->{"selected_key"}->[$i]}->{$ref->{"selected"}->{$ref->{"selected_key"}->[$i]}->[0]});
		$self->{"row"}++;
		for (my $j = 0; $j <= $#idx; $j++) {
			$self->{"mtx"}->{$i+1}->{$j} = $array[$idx[$j]];
		}
	}
}


## capitalize the column matrix organized by row and organize by column
## $self->transpose();
sub transpose {
	my $self = shift @_;
	for (my $i = 0; $i < $self->{"row"}; $i++) {
		for (my $j = 0; $j < $self->{"col"}; $j++) {
			$self->{"MTX"}->{$j}->{$i} = $self->{"mtx"}->{$i}->{$j};
			$self->{"MTX"}->{$j}->{$i} =~ tr/[a-z]/[A-Z]/;
		}
	}
}


## calculate entropy and other numbers
## $self->calculate();
sub calculate {
	my $self = shift @_;
	my $good = 0;
	my $tmpsum = 0;
	my $tmpes2 = 0;
	
	## calculate entropy for each column
	for (my $i = 0; $i < $self->{"col"}; $i++) {
		my %comp = $self->composition($i);
		$self->{"depth"}->{$i} = $comp{'depth'};
		$self->{"gaps"}->{$i} = $comp{'-'};
		if (($comp{'depth'}/$self->{"row"}) < $COVERAGE || $comp{'depth'} < 2) {
			$self->{"entropy"}->{$i} = "";
		}
		else {
			$self->{"entropy"}->{$i} = $self->calculate_entropy(\%comp);
			$good++;
			$tmpsum += $self->{"entropy"}->{$i};
			$tmpes2 += ($self->{"entropy"}->{$i}) ** 2;
		}
	}
	
	## calculate mean and std
	if ($good > 0) {
		$self->{"mean"} = $tmpsum / $good;
		$self->{"std"} = sqrt(($good * (($tmpes2 / $good) - (($self->{"mean"}) ** 2))) / ($good - 1));
	}
	else {
		$self->{"mean"} = "";
		$self->{"std"} = "";
	}
	
	## calculate Z-score
	for (my $i = 0; $i < $self->{"col"}; $i++) {
		if ($self->{"mean"} ne "" && $self->{"std"} ne "" && $self->{"entropy"}->{$i} ne "" && $self->{"std"} ne 0) {	## must guarantee all stats are meaningfull and std non-zero if assume a normal distribution
			$self->{"Z"}->{$i} = ($self->{"entropy"}->{$i} - $self->{"mean"}) / $self->{"std"};
		}
		else {
			$self->{"Z"}->{$i} = "";
		}
	}
}


## composition
## my %composition = $self->composition($col);
sub composition {
	my $self = shift @_;
	my $col = shift @_;
	
	my %aa = %Reference::standard_AA;
	foreach my $a (keys %aa) {
		$aa{$a} = 0;
	}
	$aa{'-'} = 0;
	$aa{'depth'} = 0;
	
	for (my $i = 0; $i < $self->{"row"}; $i++) {
		if (exists $aa{$self->{"MTX"}->{$col}->{$i}}) {
			$aa{$self->{"MTX"}->{$col}->{$i}}++;
			if ($self->{"MTX"}->{$col}->{$i} =~ m/[A-Z]/) {
				$aa{'depth'}++;
			}
		}
	}
	
	return %aa;
}


## calculate column entropy given the composition hash
## my $entropy = $self->calculate_entropy(\%composition);
sub calculate_entropy {
	my $self = shift @_;
	my $ref = shift @_;
	my $entropy = 0;
	foreach my $aa (keys %{$ref}) {
		if ($aa ne '-' && $aa ne 'depth') {
			my $p = $ref->{$aa} / $ref->{'depth'};
			if ($p ne 0) {
				$entropy += ($p * (log($p)/log(2)));
			}
		}
	}
	$entropy *= -1;
	return $entropy;
}


## store text of alignment of each column
## $self->text();
sub text {
	my $self = shift @_;
	for (my $i = 0; $i < $self->{"col"}; $i++) {
		my @tmp;
		for (my $j = 0; $j < $self->{"row"}; $j++) {
			push(@tmp, $self->{"MTX"}->{$i}->{$j});
		}
		$self->{"align"}->{$i} = join("", @tmp);
	}
}


## for print out
## $self->Verbose($outputfilehandler);
sub Verbose {
	my $self = shift @_;
	my $l = shift @_;
	print $l "mean entropy:\t$self->{\"mean\"}\n";
	print $l "std entropy:\t$self->{\"std\"}\n";
	my $e = 0;
	for (my $i = 0; $i < $self->{"col"}; $i++) {
		print $l ($i+1)."\t$self->{\"entropy\"}->{$i}\t$self->{\"Z\"}->{$i}\t$self->{\"align\"}->{$i}\n";
		if ($self->{"entropy"}->{$i} ne "") {
			$e++;
		}
	}
	print $l "number of residues with existing entropy:\t$e\n";
}


1
