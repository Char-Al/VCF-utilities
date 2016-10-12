#!/usr/local/bin/perl

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

use strict;
use warnings;
use Data::Dumper; # useful library to debug a programm ; print the structure and values of a variable (google it :D )


### Two VCF files
my $vcf1 ;
my $vcf_file1 = '/Users/Pau/Desktop/VCF_comparison/16M00064M.vcf';
my $vcf2 ;
#my $vcf_file2 = '/Users/Pau/Desktop/VCF_comparison/16M00064M.vcf';
my $vcf_file2 = '/Users/Pau/Desktop/VCF_comparison/full_variant_table_S1.vcf';
### One final VCF with common variants
my $vcf_f ;
my $common_var = '/Users/Pau/Desktop/VCF_comparison/common_variants.vcf';
### One VCF with non-common variants
my $vcf_n ;
my $noncommon_var = '/Users/Pau/Desktop/VCF_comparison/noncommon_variants.vcf';



### Variables needed for parsing vcf file
my %dataVCF1 ; 
my %dataVCF2 ; 


### Open VCF files and keep information in arrays

# Read the 1st file ans take information in a hash (Type {ChromPosRefAlt : (Chr, Pos, Ref, Alt)})   
open($vcf1, '<:encoding(UTF-8)', $vcf_file1) or die "Could not open file '$vcf_file1' $!";   

foreach my $line (<$vcf1>)  {
	chomp $line;
  	# Skip lines beginning with #, :, ", - and empty
 	next if ( $line =~ /(^\s*$)|(^\")|(^-)|(^:)|(^#)/ );  
 	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split("\t",$line); 			# splitting line on '\t'
 	my @list1 = ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) ;
 	if ( ($chr =~ m/(chr)([0-9])/) || ($chr =~ m/(chrom)([0-9])/) ) {
 		$chr = $2 ;
	}
 	my $variant = "$chr$pos$ref$alt" ;
    $dataVCF1{$variant} = \@list1;
    
}

close $vcf1 or die "Can't close $vcf1 $!";

#print Dumper(%dataVCF1);


# Read the 2nd file ans take information   
open($vcf2, '<:encoding(UTF-8)', $vcf_file2) or die "Could not open file '$vcf_file2' $!";

foreach my $line (<$vcf2>)  {
	chomp $line; 
  	# Skip lines beginning with #, :, ", - and empty
 	next if ( $line =~ /(^\s*$)|(^\")|(^-)|(^:)|(^#)/ );  
 	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split("\t",$line); 			# splitting line on '\t'
 	my @list2 = ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) ;
 	if ( ($chr =~ m/(chr)([0-9])/) || ($chr =~ m/(chrom)([0-9])/) ) {
 		$chr = $2 ;
	}
 	my $variant = "$chr$pos$ref$alt" ;
    $dataVCF2{$variant} = \@list2;
}

close $vcf2 or die "Can't close $vcf2 $!";

#print Dumper(%dataVCF2) ;

### Comparison between the two VCF files
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
# dataVCF($chr, $pos, $id, $ref, $alt, $qual, $filter, $info)

open($vcf_f, '>', $common_var) or die "Could not open file '$vcf_f' $!";
print $vcf_f "##fileformat=VCFv4.2
##reference=/sg/pipelineData/genome/Homo_sapiens/hg19_karyotypicOrdered/fasta/human_g1k_v37.fasta
##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
##genome-build=GRCh37
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

open($vcf_n, '>', $noncommon_var) or die "Could not open file '$vcf_n' $!"; 
print $vcf_n "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" ;


for my $var1 (keys %dataVCF1) {
	if (defined $dataVCF2{$var1}) { 
		print "$var1, $dataVCF1{$var1}\n";		
		print $vcf_f "$dataVCF1{$var1}[0]\t$dataVCF1{$var1}[1]\t$dataVCF1{$var1}[2]\t$dataVCF1{$var1}[3]\t$dataVCF1{$var1}[4]\t$dataVCF1{$var1}[5]\t$dataVCF1{$var1}[6]\t$dataVCF1{$var1}[7]\n";	
    }  

    else :
    	print $vcf_n "$dataVCF1{$var1}[0]\t$dataVCF1{$var1}[1]\t$dataVCF1{$var1}[2]\t$dataVCF1{$var1}[3]\t$dataVCF1{$var1}[4]\t$dataVCF1{$var1}[5]\t$dataVCF1{$var1}[6]\t$dataVCF1{$var1}[7]\n";	
}


close $vcf_f or die "Can't close $vcf_f $!";
close $vcf_n or die "Can't close $vcf_n $!";




=end


my $s = $word;
$s =~ s/\W//g;
my $k;
for (keys %hash){
    s/\W//g;
    if($_ eq $s){
        $k = $_;
        last;
    }
}
if(defined $k){
    # Do Something
}

$key = $word;
$key ~= s/\W//g;  # Any non-word characters are removed
if (defined $hash{$key}) { DoSomething; }






