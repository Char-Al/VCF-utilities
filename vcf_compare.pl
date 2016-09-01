#!/usr/bin/perl -w
# Authors : Charles VAN GOETHEM & Pauline SARRAROLS

# Perl general libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use Cwd;
use Time::HiRes;

use Term::ANSIColor;



##########################################################################################
##########################################################################################


my $pwd = getcwd();


##########################################################################################
##########################################################################################


# Mandatory arguments
my @files;

# Optional arguments
my $output = $pwd."/";
#my $makeGraphs = '';

# General arguments
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;


## Parse options and print usage if there is a syntax error, 
## or if usage was explicitly requested.

GetOptions(
	'f|files=s'	 		=> \@files,
	'o|output=s'		=> \$output,
#	'g|makeGraphs'		=> \$makeGraphs,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'm|man' 			=> \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test arguments
pod2usage(colored(['Red'], "Error :
	Hello, I am so lost...
	I need at least 2 VCF on input.
	Thanks !", "\n")) unless (@files >= 2);

# Check the files
foreach my $file (@files) {
	pod2usage(colored(['Red'], "Error :
	Hi, I dont understand...
	The file '$file' does not exist !
	Please check that !", "\n")) unless (-e $file);
	
	pod2usage(colored(['Red'], "Error :
	Hi, I am sorry but there is somthing wrong...
	The file '$file' is not a regulary file !
	Can you verify this please.", "\n")) unless (-f $file);
}

##########################################################################################
##########################################################################################

my %HEADER;
my %VARIANTS;

foreach my $file (@files) {
	print STDERR "Treat file : '$file'\n";
	open(VCF,$file) or die ("Error !\n");
	#$INFOS{$file};
	$file = basename($file);
	while(<VCF>) {
		chomp($_);
		if(/^#/) {
			if(/INFO=<ID=([^,]*).*Description="([^"]*).*/) {
				$HEADER{$file}{"INFO"}{$1} = $2;
			}
			if(/FORMAT=<ID=([^,]*).*Description="([^"]*).*/) {
				$HEADER{$file}{"FORMAT"}{$1} = $2;
			}
		} else {
			my @line = split("\t");
#			print Dumper \@line;
			#$VARIANTS{$file}{"$line[0]_$line[1]_$line[3]_$line[4]"} = 1 ;
			push(@{$VARIANTS{"$line[0]_$line[1]_$line[3]_$line[4]"}{"files"}},$file) ;
			$VARIANTS{"$line[0]_$line[1]_$line[3]_$line[4]"}{"lines"} = $_ ;
		}
	}
}

##########################################################################################
##########################################################################################

my %COUNT_VARIANTS;

foreach my $key (keys(%VARIANTS)){
	if( @{$VARIANTS{$key}{"files"}} == 1) {
		$COUNT_VARIANTS{${$VARIANTS{$key}{"files"}}[0]} += 1;
	} else {
		my $name = "";
		foreach my $names (@{$VARIANTS{$key}{"files"}}) {
			$name .= $names."_AND_";
		}
		$name = substr($name, 0, -5);
		$COUNT_VARIANTS{$name} += 1;
		print $VARIANTS{$key}{"lines"}."\n";
	}
}




#print Dumper \%HEADER;
print Dumper \%VARIANTS;
print Dumper \%COUNT_VARIANTS;



