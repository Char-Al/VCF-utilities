#!/usr/local/bin/perl

=head1 INFO

This script merges VCF files from sequencing data. It can produce tsv files or/and vcf files.

  Usage: perl vcf-comparison.pl -f file1.vcf,file2.vcf -o t,v -g refgenome -d folder -s sample   
Returns: Tab-delimited or VCF file containing the common variants
		 Tab-delimited or VCF file containing the unique variants  
    Env: Unix/Mac
 Author: Pauline Sararols (CHU Reims) <psararols@chu-reims.fr>
 		 Charles VAN GOETHEM (CHU Montpellier) <c-vangoethem@chu-montpellier.fr>

Version: v1.0


Requirements :
Perl 

Test :
perl test.pl --vcf_file 16M00064M.vcf,16M00064M.vcf,full_variant_table_S1.vcf -o t,v -g hg18 -s patient1


VCF file
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOUR
1       5734214 4_1     G       G[4:131190349[  5       .       SVTYPE=BND;MATEID=4_2;IMPRECISE;CIPOS=0,484;CIEND=0,538;BKDIST=-1;OCC=1;TSRDS=HS26_07415:6:1204:11449:180778#145,HS26_07416:1:2203:6853:144722#145,HS26_07416:1:2206:18596:82992#145,HS8_07129:3:1107:18448:194691#145,HS8_07129:3:2206:17776:37193#145;SVCLASS=deletion        RC:PS   0:0     0:5

=head1 SYNOPSIS

perl vcf-comparison.pl -f file1.vcf,file2.vcf -o t,v -g refgenome -d folder -s sample

Options:

  --help                   h      Brief help message
  --man                    m      Manual
  --vcf_files              f      Type of output file, v for VCF or t for TSV (only)
  --output_file            o      Name of file to write output of analysis
  --ref_genome             g      Reference Genome
  --directory              d      New Directory to run the analysis
  --sample     			   s 	  Name of the sample

=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper; 



# Options:
my $opts = setup();

my $vcf_files       = $opts->{'vcf_file'},
my $output_file     = $opts->{'output_file'};
my $ref_genome      = $opts->{'ref_genome'};
my $directory       = $opts->{'directory'};
my $output_type     = $opts->{'type'};
my $sample          = $opts->{'sample'};


# Default Values : 
if (!defined($sample)) { $sample = ""; } 


### Variables 
my @files = split(",", $vcf_files);
my @format = split(",", $output_file);
my @array_file;

### Open the VCF files to catch information in a hash
foreach my $file (@files) {
	my $vcf ;
	my $vcf_file = $file ;
	my %dataVCF ;
	# Open VCF files and keep information in a hash 
	open($vcf, '<:encoding(UTF-8)', $vcf_file) or die "Could not open file '$vcf' $!";

	foreach my $line (<$vcf>)  {
		chomp $line;
	 	next if ( $line =~ /(^\s*$)|(^\")|(^-)|(^:)|(^#)/ );  							# Skip lines beginning with #, :, ", - and empty
	 	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split("\t",$line); 	# Split line on '\t'
	 	if ( ($chr =~ m/(chr)([0-9])/) || ($chr =~ m/(chrom)([0-9])/) ) {				# Homogeneisation of the chromosomes name
	 		$chr = $2 ;
		}
		my @list = ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) ; 				# Each line = an array (list)
		# Hash = Key(var) : Value(list of characteristics)
	 	my $variant = "$chr$pos$ref$alt" ;												# Name of the key in the hash for each line (variant)
	    $dataVCF{$variant} = \@list;													# Value = array with info of each variant
	}
	# Hash = Key(file) : Value(hash with var)
	#$id_file = "FILE_$i";																# Name of the key in the hash for each file
	#$hash_file{$id_file} = \%dataVCF;
	push @array_file, \%dataVCF ;														# Value = hash with variants containing in each file
	#$i++ ;
	close $vcf or die "Can't close $vcf $!";
}

#print Dumper (@array_file);


### Comparison between the VCF files
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
# dataVCF($chr, $pos, $id, $ref, $alt, $qual, $filter, $info)


### Create a hash containing variant and their occurrences 

my @list_of_var;
foreach my $file (@array_file) {														# Create an array with all the variant id in all files 
	for my $key (keys $file) {
		push @list_of_var, $key ;
	}
}

my %hash_count_var ;																	# Count the occurrence of each variant id in the array
foreach my $id (@list_of_var) {
	my $count = count_occurrences($id, @list_of_var) ;
	$hash_count_var{$id} = $count;
}


### Select variant with at least 2 occurrences 
my @common ;
my @unique ;

foreach my $var (keys %hash_count_var) {
	if ($hash_count_var{$var} > 1) {													# Common variants are those with more than one occurrence
		push (@common, $var) ;															# These common variants are put in the list @common
	}
	else {
		push (@unique, $var) ;															# The others variants are put in the list @unique, they are in only one VCF file
	}
}


### FOR COMMON VARIANTS

my $vcf_f ;
my $common_var_vcf = $sample.'common_variants.vcf';
my $tsv_f ;
my $common_var_tsv = $sample.'common_variants.tsv';


### Recover characteristics for common variants
my %table_common ;
foreach my $file (@array_file) {														# We are looking for characteristics of variants of the list @common
	for my $key (keys $file) {															# We recover these characteritics thanks to %hash_file
		if ( $key ~~ @common) {															# If the $key (var_id) is in the hash, we can keep it in another hash
			$table_common{$key} = $$file{$key};
		}
	}
} 


### Write final files
																					# List for already seen variant					
for my $letter (@format) {
	## TSV File : coords of variants + count of occurrences 
	#if (grep $_ eq "t", @format) {
	if ($letter eq "t") {																# Construction of a TSV file if "t" is put in @output_file
		open($tsv_f, '>', $common_var_tsv) or die "Could not open file '$tsv_f' $!";
		get_tsv_format($tsv_f) ;
		for my $var_common (sort keys %table_common) {
			print "$var_common\n" ;													 
			for my $var_hash (sort keys %hash_count_var) {									# To add the count of occurence in the TSV file, 
				if ($var_common eq $var_hash) {											# We look for identical variant id in both hash to print it in the file
					print $tsv_f "$table_common{$var_common}[0]\t$table_common{$var_common}[1]\t$table_common{$var_common}[2]\t$table_common{$var_common}[3]\t$table_common{$var_common}[4]\t$table_common{$var_common}[5]\t$table_common{$var_common}[6]\t$table_common{$var_common}[7]\t$hash_count_var{$var_hash}\n";
				}
			}														
		}
		close $tsv_f or die "Can't close $tsv_f $!";													
	}

	## VCF File : coords of variants  
	elsif ($letter eq "v") { 															# Construction of a VCF file if "v" is put in @output_file
		open($vcf_f, '>', $common_var_vcf) or die "Could not open file '$vcf_f' $!"; 
		get_vcf_format($vcf_f);	
		for my $var_com (sort keys %table_common) {											# Each variant characteristics are printed in the VCF file
			print $vcf_f "$table_common{$var_com}[0]\t$table_common{$var_com}[1]\t$table_common{$var_com}[2]\t$table_common{$var_com}[3]\t$table_common{$var_com}[4]\t$table_common{$var_com}[5]\t$table_common{$var_com}[6]\t$table_common{$var_com}[7]\n";												
		}
		close $vcf_f or die "Can't close $vcf_f $!";
	}
}



### FOR UNIQUE VARIANTS

## One VCF with unique variants
my $vcf_u ;
my $unique_var_vcf = $sample.'unique_variants.vcf';
my $tsv_n ;
my $unique_var_tsv = $sample.'unique_variants.tsv';
my %table_unique ;

foreach my $file (@array_file) {
	for my $key (keys $file) {
		if ( $key ~~ @unique) {
			# "$key\t$$file{$key}\n";
			$table_unique{$key} = $$file{$key};
		}
	}
} 

### Write final files

for my $letter (@format) {
	## TSV File : coords of variants + count of occurrences 
	#if (grep $_ eq "t", @format) {
	if ($letter eq "t") {
		open($tsv_n, '>', $unique_var_tsv) or die "Could not open file '$tsv_n' $!";
		get_tsv_format($tsv_n) ;
		for my $var_uniq (sort keys %table_unique) {													# 
			for my $var_hash (keys %hash_count_var) {
				if ($var_uniq eq $var_hash) {
					print $tsv_n "$table_unique{$var_uniq}[0]\t$table_unique{$var_uniq}[1]\t$table_unique{$var_uniq}[2]\t$table_unique{$var_uniq}[3]\t$table_unique{$var_uniq}[4]\t$table_unique{$var_uniq}[5]\t$table_unique{$var_uniq}[6]\t$table_unique{$var_uniq}[7]\t$hash_count_var{$var_hash}\n";
				}
			}														
		}
		close $tsv_n or die "Can't close $tsv_n $!";													
	}

	## VCF File : coords of variants  
	elsif ($letter eq "v") {
		open($vcf_u, '>', $unique_var_vcf) or die "Could not open file '$vcf_u' $!"; 
		get_vcf_format($vcf_u);	
		for my $var_u (sort keys %table_unique) {												# 
			print $vcf_u "$table_unique{$var_u}[0]\t$table_unique{$var_u}[1]\t$table_unique{$var_u}[2]\t$table_unique{$var_u}[3]\t$table_unique{$var_u}[4]\t$table_unique{$var_u}[5]\t$table_unique{$var_u}[6]\t$table_unique{$var_u}[7]\n";												
		}
		close $vcf_u or die "Can't close $vcf_u $!";
	}
}

## Sort files
#my $output_vcf = "sorted.".$common_var_vcf ;
#get_sorted_VCF($common_var_vcf, $output_vcf) ;


#my $output_com_tsv = "sorted.".$common_var_tsv ;
#get_sorted_TSV($common_var_tsv, $output_com_tsv) ;

#my $output_uni_tsv = "sorted.".$unique_var_tsv ;
#get_sorted_TSV($unique_var_tsv, $output_uni_tsv) ;


### VENN DIAGRAM




###  *** Functions ***


sub setup  {
    my %opts;
    my @random_args;
    GetOptions(
    'h|help'                   => \$opts{'h'},
    'm|man'                    => \$opts{'m'},
    'f|vcf_files=s'            => \$opts{'vcf_file'},
    'o|output_file=s'          => \$opts{'output_file'},
    'g|ref_genome=s'           => \$opts{'ref_genome'},  
    'd|directory=s'            => \$opts{'directory'},
    's|sample=s'               => \$opts{'sample'},
    '<>'                       => sub{push(@random_args,shift(@_));}
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  return \%opts;
}


sub count_occurrences {
	my ($word, @array) = @_ ;
	my $count_word = grep { $_ eq $word } @array;
	return $count_word ;
}


sub get_vcf_format {  
    my $vcf = $_[0] ;
    print $vcf "##fileformat=VCFv4.2
##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
##genome-build=$ref_genome
##technology-platform-name=Illumina
##technology-platform-read-type=pair
##FILTER=<ID=low_coverage,Description=\"Filtered due to low coverage\">
##FILTER=<ID=low_variant_fraction,Description=\"Filtered due to low variant fraction\">
##FILTER=<ID=homopolymer_region,Description=\"Filtered because in homopolymer region\">
##FILTER=<ID=problematic_region,Description=\"Filtered because of noisy region, e.g. CG-rich region or low complexity region\">
##FILTER=<ID=off_target,Description=\"Filtered because variant is outside of the target region\">
##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Indicates whether variant is SNP | INDEL.\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">
##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"#ref plus strand, #ref minus strand, #alt plus strand, #alt minus strand\">
##INFO=<ID=AD,Number=2,Type=Integer,Description=\"Allelic depth for the ref and alt alleles in the order listed\">
##INFO=<ID=DBXREF,Number=.,Type=String,Description=\"Colon-separated key-value pairs of overlaps in database e.g. dbSNP:rs838532,COSMIC:COSM28362\">
##INFO=<ID=OVLP,Number=.,Type=String,Description=\"Overlap with relevant database
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" ;
}


sub get_tsv_format {
	my $tsv = $_[0] ;
    print $tsv "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tOCCUR\n" ;
}






=end


sub compare_2_vcf {
	my (%dataVCF1, %dataVCF2) = @_ ;
	for my $var (keys %dataVCF1) {
		if (defined $dataVCF2{$var}) { 
			#print "$var1, $dataVCF1{$var}\n";
			if (grep { $_ eq "t" } @format) {
				print $tsv_f "$dataVCF1{$var}[0]\t$dataVCF1{$var}[1]\t$dataVCF1{$var}[2]\t$dataVCF1{$var}[3]\t$dataVCF1{$var}[4]\t$dataVCF1{$var}[5]\t$dataVCF1{$var}[6]\t$dataVCF1{$var}[7]\n";

			}		
			elsif (grep { $_ eq "v" } @format) {
				print $vcf_f "$dataVCF1{$var}[0]\t$dataVCF1{$var}[1]\t$dataVCF1{$var}[2]\t$dataVCF1{$var}[3]\t$dataVCF1{$var}[4]\t$dataVCF1{$var}[5]\t$dataVCF1{$var}[6]\t$dataVCF1{$var}[7]\n";
			}	
	    }  

	    else :
			if (grep { $_ eq "t" } @format) {
				print $tsv_n "$dataVCF1{$var}[0]\t$dataVCF1{$var}[1]\t$dataVCF1{$var}[2]\t$dataVCF1{$var}[3]\t$dataVCF1{$var}[4]\t$dataVCF1{$var}[5]\t$dataVCF1{$var}[6]\t$dataVCF1{$var}[7]\n";

			}		
			elsif (grep { $_ eq "v" } @format) {
				print $vcf_u "$dataVCF1{$var}[0]\t$dataVCF1{$var}[1]\t$dataVCF1{$var}[2]\t$dataVCF1{$var}[3]\t$dataVCF1{$var}[4]\t$dataVCF1{$var}[5]\t$dataVCF1{$var}[6]\t$dataVCF1{$var}[7]\n";
			}
		}
	}	
}


sub list_of_id {
	my %dataVCF = $_[0] ;
	my @list_of_id = (keys %dataVCF) ;
	return @list_of_id ;
}


sub uniq {
	my @list = $_[0];
	my %tab ;
	foreach my $item (@list) {
		$tab{$item} = 1;
	}
	my @uniq = (keys %tab) ;
	return @uniq ;
}



