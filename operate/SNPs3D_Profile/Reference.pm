package Reference;
use strict;
use warnings;
use base qw(Exporter);
our @EXPORT_OK;

### copied from /moulthome/yini/work/project/scripts/self/utility ###
### 26JAN2018 YY ###

#%standard_AA
#%standard_aa
#%standard_AAC
#%standard_aaC
#%standard_sasa_Chothia
#%standard_sasa_Sander
#%standard_sasa_Thornton
#% codon_PN;
#% codon_NP;
sub codon_NP;


our %standard_AA = (	'A' => 'Ala',
			'C' => 'Cys',
			'D' => 'Asp',
			'E' => 'Glu',
			'F' => 'Phe',
			'G' => 'Gly',
			'H' => 'His',
			'I' => 'Ile',
			'K' => 'Lys',
			'L' => 'Leu',
			'M' => 'Met',
			'N' => 'Asn',
			'O' => 'Pyl',
			'P' => 'Pro',
			'Q' => 'Gln',
			'R' => 'Arg',
			'S' => 'Ser',
			'T' => 'Thr',
			'U' => 'Sec',
			'V' => 'Val',
			'W' => 'Trp',
			'Y' => 'Tyr'
);

our %standard_aa = (	'Ala' => 'A',
			'Cys' => 'C',
			'Asp' => 'D',
			'Glu' => 'E',
			'Phe' => 'F',
			'Gly' => 'G',
			'His' => 'H',
			'Ile' => 'I',
			'Lys' => 'K',
			'Leu' => 'L',
			'Met' => 'M',
			'Asn' => 'N',
			'Pyl' => 'O',
			'Pro' => 'P',
			'Gln' => 'Q',
			'Arg' => 'R',
			'Ser' => 'S',
			'Thr' => 'T',
			'Sec' => 'U',
			'Val' => 'V',
			'Trp' => 'W',
			'Tyr' => 'Y'
);

our %standard_AAC = (	'A' => 'ALA',
                        'C' => 'CYS',
                        'D' => 'ASP',
                        'E' => 'GLU',
                        'F' => 'PHE',
                        'G' => 'GLY',
                        'H' => 'HIS',
                        'I' => 'ILE',
                        'K' => 'LYS',
                        'L' => 'LEU',
                        'M' => 'MET',
                        'N' => 'ASN',
			'O' => 'PYL',
                        'P' => 'PRO',
                        'Q' => 'GLN',
                        'R' => 'ARG',
                        'S' => 'SER',
                        'T' => 'THR',
			'U' => 'SEC',
                        'V' => 'VAL',
                        'W' => 'TRP',
                        'Y' => 'TYR'
);

our %standard_aaC = (	'ALA' => 'A',
                        'CYS' => 'C',
                        'ASP' => 'D',
                        'GLU' => 'E',
                        'PHE' => 'F',
                        'GLY' => 'G',
                        'HIS' => 'H',
                        'ILE' => 'I',
                        'LYS' => 'K',
                        'LEU' => 'L',
                        'MET' => 'M',
                        'ASN' => 'N',
			'PYL' => 'O',
                        'PRO' => 'P',
                        'GLN' => 'Q',
                        'ARG' => 'R',
                        'SER' => 'S',
                        'THR' => 'T',
			'SEC' => 'U',
                        'VAL' => 'V',
                        'TRP' => 'W',
                        'TYR' => 'Y'
);


our %standard_sasa_Chothia = (	'A' => 115,
				'C' => 135,
				'D' => 150,
				'E' => 190,
				'F' => 210,
				'G' => 75,
				'H' => 195,
				'I' => 175,
				'K' => 200,
				'L' => 170,
				'M' => 185,
				'N' => 160,
				'P' => 145,
				'Q' => 180,
				'R' => 225,
				'S' => 115,
				'T' => 140,
				'V' => 155,
				'W' => 255,
				'Y' => 230
);


our %standard_sasa_Sander = (	'A' => 106,
				'B' => 160,	## stands for D or N
				'C' => 135,
				'D' => 163,
				'E' => 194,
				'F' => 197,
				'G' => 84,
				'H' => 184,
				'I' => 169,
				'K' => 205,
				'L' => 164,
				'M' => 188,
				'N' => 157,
				'P' => 136,
				'Q' => 198,
				'R' => 248,
				'S' => 130,
				'T' => 142,
				'V' => 142,
				'W' => 227,
				'X' => 180,	## stands for undetermined
				'Y' => 222,
				'Z' => 196	## stands for Q or X
);


## all Thornton standard sasa were calculated from NACCESS output file examples, except W ##
our %standard_sasa_Thornton = (	'A' => 108,
				'C' => 160,
				'D' => 140,
				'E' => 172,
				'F' => 199,
				'G' => 80,
				'H' => 183,
				'I' => 175,
				'K' => 201,
				'L' => 179,
				'M' => 194,
				'N' => 144,
				'P' => 136,
				'Q' => 179,
				'R' => 239,
				'S' => 117,
				'T' => 139,
				'V' => 151,
				'W' => 241,	## get from "Protein Bioinformatics: From sequqnce to Function (By M. Michael Gromiha) page 81"
				'Y' => 213
);


our %codon_PN = (	'A' => ['GCT', 'GCC', 'GCA', 'GCG'],
			'R' => ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
			'N' => ['AAT', 'AAC'],
			'D' => ['GAT', 'GAC'],
			'C' => ['TGT', 'TGC'],
			'Q' => ['CAA', 'CAG'],
			'E' => ['GAA', 'GAG'],
			'G' => ['GGT', 'GGC', 'GGA', 'GGG'],
			'H' => ['CAT', 'CAC'],
			'I' => ['ATT', 'ATC', 'ATA'],
			#'START' => ['ATG'],
			'L' => ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
			'K' => ['AAA', 'AAG'],
			'M' => ['ATG'],
			'F' => ['TTT', 'TTC'],
			'P' => ['CCT', 'CCC', 'CCA', 'CCG'],
			'S' => ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
			'T' => ['ACT', 'ACC', 'ACA', 'ACG'],
			'W' => ['TGG'],
			'Y' => ['TAT', 'TAC'],
			'V' => ['GTT', 'GTC', 'GTA', 'GTG'],
			'X' => ['TAA', 'TGA', 'TAG']
);


sub codon_NP {
	my $ref = shift @_;
	my %data;
	foreach my $aa (keys %{$ref}) {
		foreach my $nn (@{$ref->{$aa}}) {
			my $udna = $nn;
			my $ldna = lc $udna;
			my $urna = $udna;
			$urna =~ s/T/U/g;
			my $lrna = lc $urna;
			$data{$udna} = $aa;
			$data{$ldna} = $aa;
			$data{$urna} = $aa;
			$data{$lrna} = $aa;
		}
	}
	return %data;
}

our %codon_NP = codon_NP(\%codon_PN);



1;	

