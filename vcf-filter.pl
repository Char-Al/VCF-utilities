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
use Scalar::Util qw(looks_like_number);

use feature "switch";

##########################################################################################
##########################################################################################


my $pwd = getcwd();


##########################################################################################
##########################################################################################


# Mandatory arguments
my @files;

# Optional arguments
my $output = "";
my $details;
my $pass;
my $config = "";
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
	'd|details'			=> \$details,
	'c|config=s'		=> \$config,
	'p|pass-only'		=> \$pass,
#	'g|makeGraphs'		=> \$makeGraphs,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'm|man' 			=> \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test arguments

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

# Check the output
if ($output ne "") {
	pod2usage(colored(['Red'], "Error :
	There is something wrong...
	The repertory '$output' does not exist !
	You may create it before !", "\n")) unless (-e $output);

	pod2usage(colored(['Red'], "Error :
	I dont think it is normal...
	The repertory '$output' is not a repertory !
	I will cry if you dont fix this :'(", "\n")) unless (-d $output);
} else {
	$output = $pwd."/";
	print STDERR colored(['Yellow'], "Warning : All the outputs will be create here :
	$output","\n");
}


##########################################################################################
##########################################################################################

my %hash_filter;
open(CONF,"$config") or die ("cannot open $config $!");
my $i = 1;
while(<CONF>){
	chomp();
	my @line = split("\t");

	$hash_filter{$line[4]}{"cond"} = $line[5];
	$hash_filter{$line[4]}{"step"} = $line[6];
	$hash_filter{$line[4]}{"rules"}{$i}{"id"} = $line[0];
	$hash_filter{$line[4]}{"rules"}{$i}{"val"} = $line[1];
	$hash_filter{$line[4]}{"rules"}{$i}{"is"} = $line[2];
	$hash_filter{$line[4]}{"rules"}{$i}{"for"} = $line[3];

	$i++;
}

#print Dumper \%hash_filter;


##########################################################################################
##########################################################################################

my %VARIANTS;

foreach my $file (@files) {
	print STDERR "Treat file : '$file'\n";
	open(VCF,$file) or die ("Error !\n");
	#$INFOS{$file};
	$file = basename($file);
	my @all_lines;
	while(<VCF>) {
		chomp($_);
		if($_ !~ /^#/) {
			my $line = $_;
			#print colored(['Blue'], "$_","\n");
			my %hash_line = parse_line($_);
			foreach my $key (sort keys(%hash_filter)) {
				my %boolean_rules = (
					1 => 0,
					0 => 0
				);
				#print Dumper \%boolean_rules;
				foreach my $rules (keys $hash_filter{$key}{"rules"}) {
					$boolean_rules{check_filter($hash_filter{$key}{"rules"}{$rules},\%hash_line)}++;
				}

				my $verify = check_condition($hash_filter{$key}{"cond"},$boolean_rules{1},$boolean_rules{0}) ;
				#print $verify."\n";
				#print Dumper \%boolean_rules;
				my $check ="";
				for ($hash_filter{$key}{"step"}) {
					when("stop")	{ if ($verify) {print $line."\n"; $check = "next_line"; } else { $check = "next_line"} }
					when("get")		{ if ($verify) {print $line."\n"; $check = "next_line"; } else { next ; }}
					when("next")	{ if($verify) { next; }else { $check = "next_line"} }
				}
				if ($check eq "next_line") {
					last;
				}


				#check_condition($hash_filter{$key}{"cond"},@boolean_rules);

				#if ($hash_filter{$key}{"step"} eq "get" ) {

					#print $_;
					#last;
				#} elsif ($hash_filter{$key}{"step"} eq "next" ) {
				#	if ($hash_filter{$key}{"cond"} eq "OR") {
				#	}
					#print $_;
					#next;
				#}
				#print $key."\n";
			}
			#print Dumper \%hash_line;
			#push(@all_lines,\%hash_line);
		}
	}
	#print Dumper \@all_lines;
	close(VCF);
}

##########################################################################################
##########################################################################################


sub parse_line {
	my $line = shift(@_);
	my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info, $format, @formats) = split("\t",$line);
	my @infos = split(";",$info);

	my %hash_line = (
		"chr"	=> $chr,
		"pos"	=> $pos,
		"ID"	=> $ID,
		"ref"	=> $ref,
		"alt"	=> $alt,
		"qual"	=> $qual,
		"filter"	=> $filter
	);
	#print Dumper \%hash_line;

	# Add each each formats
	my $i=1;
	my @all_formats = split(":",$format);
	#$hash_line{"format"}{"format"} = $format;
	foreach my $element (@formats) {
		my @elements = split(":",$element);
		for my $j (0..(@all_formats-1)){
			$hash_line{"format"}{"indiv".$i}{$all_formats[$j]} = $elements[$j];
		}
		$i++;
	}

	# Add each infos
	foreach my $value (@infos) {
		my @values = split("=",$value);
		$hash_line{"infos"}{$values[0]} = $values[1];
	}

	#print Dumper \%hash_line;
	return %hash_line;
}

sub check_filter {
	my $hash_1 = shift(@_);
	my $hash_2 = shift(@_);
	#my %line = shift(@_);
	my %hash_rule = %$hash_1;
	my %hash_line = %$hash_2;
#	print Dumper \%hash_rule;
#	print Dumper \%hash_line;

	for($hash_rule{"is"} )
	{
		when("num") {
			my $val = $hash_line{"infos"}{$hash_rule{"id"}};
			if ($val eq ".") {
				$val = 0;
			}
			for ($hash_rule{"for"}) {
				when("eq")		{ $val == $hash_rule{"val"} ? return 1 : return 0 ; }
				when("supeq")	{ $val >= $hash_rule{"val"} ? return 1 : return 0 ; }
				when("infeq")	{ $val <= $hash_rule{"val"} ? return 1 : return 0 ; }
				when("sup")		{ $val >  $hash_rule{"val"} ? return 1 : return 0 ; }
				when("inf")		{ $val <  $hash_rule{"val"} ? return 1 : return 0 ; }
				when("ne")		{ $val != $hash_rule{"val"} ? return 1 : return 0 ; }
				default			{return 0}
			}
		}
		when("str") {
			for ($hash_rule{"for"}) {
				#print $hash_line{"infos"}{$hash_rule{"id"}}."\n";
				when("eq")		{ $hash_line{"infos"}{$hash_rule{"id"}} eq $hash_rule{"val"} ? return 1 : return 0 ; }
				when("ne")		{ $hash_line{"infos"}{$hash_rule{"id"}} ne $hash_rule{"val"} ? return 1 : return 0 ; }
				default			{return 0}
			}
		}
	}
	#$hash_line{"infos"}{$hash_rule{"id"}};
}

sub check_condition {
	my $cond = shift(@_);
	my $true = shift(@_);
	my $false = shift(@_);

	#print "cond : $cond \n";
	#print "nb true : $true \n";
	#print "nb false : $false \n";

	for($cond ) {
		when("OR")	{ ($true >= 1) ? return 1 : 0 ; }
		when("AND")	{ ($false == 0) ? return 1 : 0 ; }
	}
	#$hash_line{"infos"}{$hash_rule{"id"}};
}

##########################################################################################
##########################################################################################

__END__

=pod

=encoding UTF-8

=head1 NAME

vcf-filter.pl - Filter a vcf file with a list of decision on a conf file

=head1 VERSION

version 0.02 - alpha test without bugs -- "last_but_not_least"

=head1 SYNOPSIS

compare_leon.pl  -f file.fastq (-f file2.fastq,file3.fastq) -d directory/with/some/fastq

=head1 DESCRIPTION

This script compress and uncompress automatically some fastQ file.

=head1 OPTIONS

=head2 General

	-h,--help		Print this help
	-m,--man		Open man page
	-v,--verbosity		Level of verbosity

=head2 Mandatory arguments

	-f,--file=file.fastq			Specify the fastq file you want use (possible multiple file)
	-d,--directory=path/to/directory	Specify a directory contains some fastq you want use (possible multiple directory)

=head2 Optional arguments

	-o,--output=repertory			You can specify the output repertory (default Current)
	-g,--makeGraphs					Generate automatically awesome graphs with R and ggplot2
	-p,--pass						Filter only on PASS variants

=head1 AUTHORS

=over 4

=item -
Charles VAN GOETHEM

=item -
Pauline SARRAROLS

=back

=cut
