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
use File::stat;

use VCF;


##########################################################################################
##########################################################################################


my $pwd = getcwd();


##########################################################################################
##########################################################################################

# Mandatory arguments
my $file;
my $conf;

# Optional arguments
my $output = "";

# General arguments
my $man		= 0;
my $help	= 0;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.

GetOptions(
	'f|file=s'		=> \$file,
	'c|conf=s'		=> \$conf,
	'o|output=s'	=> \$output,
	'help|?'		=> \$help,
	'm|man'			=> \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test arguments
pod2usage(-verbose => 1, -message => colored(['Red'], "Error :	The file '$file' does not exist !", "\n")) unless (-e $file);
pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The file '$file' is not a regulary file !", "\n")) unless (-f $file);
my ($name, $dir, $ext) = fileparse($file,qw(.vcf));
pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The file '$file' is not with 'vcf' extension.", "\n")) unless ($ext);

pod2usage(-verbose => 1, -message => colored(['Red'], "Error :	The file '$conf' does not exist !", "\n")) unless (-e $conf);
pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The file '$conf' is not a regulary file !", "\n")) unless (-f $conf);
($name, $dir, $ext) = fileparse($conf,qw(.conf));
pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The file '$conf' is not with 'conf' extension.", "\n")) unless ($ext);



if ($output ne "") {
	pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The repertory '$output' does not exist !", "\n")) unless (-e $output);
	pod2usage(-verbose => 1, -message => colored(['Red'], "Error : The repertory '$output' is not a repertory !", "\n")) unless (-d $output);
} else {
	$output = $pwd."/";
	print STDERR colored(['Yellow'], "Warning : All the outputs will be create here : '$output'","\n");
}


##########################################################################################
##########################################################################################




($name, $dir, $ext) = fileparse($file, qr/\.[^.]*/);
my $output_tab = $output.$name.".filter.tab";
my $output_vcf = $output.$name.".filter.vcf";

my $header = VCF->read_header($file);
my $hash_cond = VCF->parseConfFile($conf);
my $filter = VCF->vcfCondFilter($header, $file, $hash_cond, $output_vcf, $output_tab);




#########################################################################################
##########################################################################################

__END__

=pod

=encoding UTF-8

=head1 NAME

vcf-cond-filter.pl - Filter a vcf based on a list of conditions and output on a tab separate file.

=head1 VERSION

version 0.01

=head1 USAGE

    perl vcf-cond-filter.pl  -f file1.vcf -c filter.conf [-o output_repertory]

=head1 DESCRIPTION

This script filter all variants of a vcf based on a list of a condition.
The variants filtered are print on a tab separate file.

=head1 OPTIONS

=head2 Mandatory arguments

	-f,--file file.vcf		vcf to convert
	-c,--conf filter.conf	rules of the filter

=head2 Optional arguments

	-o,--output /path/output	output repertory (Default : current)

=head2 General

	-h,--help	Print this help
	-m,--man	Open man page

=head1 AUTHORS

=over 1

=item -
Charles VAN GOETHEM (CHU - Montpellier)

L<c-vangoethem@chu-montpellier.fr>

=item -
Pauline SARRAROLS (CHU - Reims)

L<psararols@chu-reims.fr>

=item -
Charly MATHIEU (Université de Montpellier 1 - Montpellier)

=back

=head1 LICENSE

Copyright (c) 2016 c-vangoethem, p-sararols

This file is part of VCF-utilities.

VCF-utilities is free tools: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VCF-utilities is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VCF-utilities. If not, see <http://www.gnu.org/licenses/>.

=cut
