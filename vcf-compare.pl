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
use File::stat;

use VCF;


##########################################################################################
##########################################################################################


my $pwd = getcwd();


##########################################################################################
##########################################################################################


# Mandatory arguments
my @files;

# Optional arguments
my $output = "";

# General arguments
my $man 		= 0;
my $help 		= 0;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.

GetOptions(
	'f|files=s'	 		=> \@files,
	'o|output=s'		=> \$output,
	'help|?' 			=> \$help,
	'm|man' 			=> \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test arguments
pod2usage(-verbose => 1, -message => colored(['Red'], "Error : two vcf at least are expected !", "\n")) unless (@files >= 2);

@files = split(/,/,join(',',@files));

foreach my $file (@files) {
	pod2usage(-verbose => 1, -message => colored(['Red'], "Error :	The file '$file' does not exist !", "\n")) unless (-e $file);
	pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The file '$file' is not a regulary file !", "\n")) unless (-f $file);
	my ($name, $dir, $ext) = fileparse($file,qw(.vcf));
	pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The file '$file' is not with 'vcf' extension.", "\n")) unless ($ext);
}

if ($output ne "") {
	pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The repertory '$output' does not exist !", "\n")) unless (-e $output);
	pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The repertory '$output' is not a repertory !", "\n")) unless (-d $output);
} else {
	$output = $pwd."/";
	print STDERR colored(['Yellow'], "Warning : All the outputs will be create here : '$output'","\n");
}


##########################################################################################
##########################################################################################


my @vcfs;
foreach my $file (@files) {
    my ($name, $dir, $ext) = fileparse($file,, qr/\.[^.]*/);
    #$file = basename($file);
    my $output_file = $output.$name.".tab";
    my $var = VCF->speed_comp($file);
	#print Dumper $var;
	push(@vcfs,$var);
}



print Dumper @vcfs;











#
# my %HEADER;
# my %VARIANTS;
#
# foreach my $file (@files) {
# 	print STDERR "Treat file : '$file'\n";
# 	open(VCF,$file) or die ("Error !\n");
# 	#$INFOS{$file};
# 	$file = basename($file);
# 	my @all_lines;
# 	while(<VCF>) {
# 		chomp($_);
# 		if(/^#/) {
# 			if(/INFO=<ID=([^,]*).*Description="([^"]*).*/) {
# 				$HEADER{$file}{"INFO"}{$1}{"desc"} = $2;
# 				$HEADER{$file}{"INFO"}{$1}{"line"} = $_;
# 			}
# 			if(/FORMAT=<ID=([^,]*).*Description="([^"]*).*/) {
# 				$HEADER{$file}{"FORMAT"}{$1}{"desc"} = $2;
# 				$HEADER{$file}{"FORMAT"}{$1}{"line"} = $_;
# 			}
# 		} else {
# 			my @line = split("\t");
# #			print Dumper \@line;
# 			#$VARIANTS{$file}{"$line[0]_$line[1]_$line[3]_$line[4]"} = 1 ;
# 			unless ($pass){
# 				push(@{$VARIANTS{"$line[0]_$line[1]_$line[3]_$line[4]"}{"files"}},$file) ;
# 				$VARIANTS{"$line[0]_$line[1]_$line[3]_$line[4]"}{"lines"} = $_ ;
# 			} else {
# 				if ($line[6] eq "PASS" && length($line[3]) == 1 && length($line[4]) == 1) {
# 					push(@{$VARIANTS{"$line[0]_$line[1]_$line[3]_$line[4]"}{"files"}},$file) ;
# 					$VARIANTS{"$line[0]_$line[1]_$line[3]_$line[4]"}{"lines"} = $_ ;
# 				}
# 			}
# 			my %hash_line = parse_line($_);
# 			#print Dumper \%hash_line;
# 			push(@all_lines,\%hash_line);
# 		}
# 	}
# 	#print Dumper \@all_lines;
# 	close(VCF);
# }
#
# ##########################################################################################
# ##########################################################################################
#
#
# my %TREAT_VARIANTS;
# my %VCF;
# foreach my $key (keys(%VARIANTS)){
# 	if( @{$VARIANTS{$key}{"files"}} == 1) {
# 		my $name = substr(${$VARIANTS{$key}{"files"}}[0],0,-4);
# 		$TREAT_VARIANTS{$name}{"count"} += 1;
# 		$TREAT_VARIANTS{$name}{"vcf"} .= $VARIANTS{$key}{"lines"}."\n";
# 	} else {
# 		my $name = "common_";
# 		foreach my $names (@{$VARIANTS{$key}{"files"}}) {
# 			$name .= substr($names,0,-4)."_AND_";
# 		}
# 		$name = substr($name, 0, -5);
# 		$TREAT_VARIANTS{$name}{"count"} += 1;
# 		$TREAT_VARIANTS{$name}{"vcf"} .= $VARIANTS{$key}{"lines"}."\n";
# 	}
# }
#
# # print Dumper \%HEADER;
# # print Dumper \%VARIANTS;
# # print Dumper \%TREAT_VARIANTS;
#
# foreach my $key (keys(%TREAT_VARIANTS)) {
# 	open(FILE, ">".$output.$key.".compare") or die("Cannot open $output$key");
# 	print FILE $TREAT_VARIANTS{$key}{"vcf"};
# 	print $key." : ".$TREAT_VARIANTS{$key}{"count"}."\n";
# }
#
# ##########################################################################################
# ##########################################################################################
#
#
# sub parse_line {
# 	my $line = shift(@_);
# 	my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info, $format, @formats) = split("\t",$line);
# 	my @infos = split(";",$info);
#
# 	my %hash_line = (
# 		"chr"	=> $chr,
# 		"pos"	=> $pos,
# 		"ID"	=> $ID,
# 		"ref"	=> $ref,
# 		"alt"	=> $alt,
# 		"qual"	=> $qual,
# 		"filter"	=> $filter
# 	);
# 	#print Dumper \%hash_line;
#
# 	# Add each each formats
# 	my $i=1;
# 	my @all_formats = split(":",$format);
# 	#$hash_line{"format"}{"format"} = $format;
# 	foreach my $element (@formats) {
# 		my @elements = split(":",$element);
# 		for my $j (0..(@all_formats-1)){
# 			$hash_line{"format"}{"indiv".$i}{$all_formats[$j]} = $elements[$j];
# 		}
# 		$i++;
# 	}
#
# 	# Add each infos
# 	foreach my $value (@infos) {
# 		my @values = split("=",$value);
# 		$hash_line{"infos"}{$values[0]} = $values[1];
# 	}
#
# 	#print Dumper \%hash_line;
# 	return %hash_line;
# }


##########################################################################################
##########################################################################################

__END__

=pod

=encoding UTF-8

=head1 NAME

compare_leon.pl - Compare Leon compression algorithm with gzip (generally used) for fastQ.

=head1 VERSION

version 0.01

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

=head1 AUTHORS

=over 4

=item -
Charles VAN GOETHEM

=item -
Pauline SARRAROLS

=back

=cut
