package VCF;

use Data::Dumper;
use File::Basename;
use Hash::Merge qw(merge);


=head2 read_vcf
    About   : Reads a VCF and tranform it to an hash
    Usage   : my $x = read_vcf($file_name);
    Args    : name of the vcf file
=cut
sub read_vcf {
	my ($it,$file_name) = @_;

	print STDERR "Treat file : '$file_name'\n";
	open(VCF,$file_name) or die ("Error !\n");
	$file_name = basename($file_name);
	my @all_lines;
	my $hash_variants = {};
	my $hash_header = {};
	my $i = 0;
	while(<VCF>) {
		chomp($_);
		my ($hash_line,$type) = parse_line($_,$hash_header);
		if($type eq "HEADER") {
			$hash_header = merge($hash_header,$hash_line);
		} else {
			$i++;
			$hash_variants = merge($hash_variants,{ "variants" => {"variant_".$i => \$hash_line}});
		}
	}
	close(VCF);
	return merge($hash_header,$hash_variants);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

=head2 parse_line
    About   : Reads a VCF line and splits it into a hash.
    Usage   : my $x = parse_line($line);
    Args    : line to parse.
=cut
sub parse_line {
	my ($line,$hash_header) = @_;
	## Control line
	if ( !$line ) { return undef; }

	my %hash_line;
	if($line =~ /^#/) {
		%hash_line = parse_meta($line);
		return (\%hash_line,"HEADER");
	} else {
		my $ncols = scalar(split("\t",$line));
		if($ncols > 8) {
			#print Dumper \$hash_header;
			my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info, $format, @formats) = split("\t",$line);%hash_line = (
				"CHROM"		=> $chr,
				"POS"		=> $pos,
				"ID"		=> $ID,
				"REF"		=> $ref,
				"ALT"		=> $alt,
				"QUAL"		=> $qual,
				"FILTER"	=> $filter,
				"INFO"		=> parse_info($info, $hash_header),
				"FORMAT"	=> parse_format($format, $hash_header, @formats)
			);
		} else {
			my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info) = split("\t",$line);%hash_line = (
				"CHROM"		=> $chr,
				"POS"		=> $pos,
				"ID"		=> $ID,
				"REF"		=> $ref,
				"ALT"		=> $alt,
				"QUAL"		=> $qual,
				"FILTER"	=> $filter,
				"INFO"		=> parse_info($info, $hash_header)
			);
		}

		return (\%hash_line,"Var");

	}
	#print Dumper \%hash_line;
}

=head2 parse_info
    About   : Parse field info
    Usage   : my $x = parse_info($info);
    Args    : field to parse.
=cut
sub parse_info {
	my ($info,$hash_header) = @_;
	my @infos = split(";",$info);

	my %hash_info;

	foreach my $value (@infos) {
		my @values = split("=",$value);
		#print $values[0]."\n";
		if ($hash_header->{'meta-informations'}->{"INFO"}->{$values[0]}->{"Type"} eq "Flag") {
			$hash_info{$values[0]} = "TRUE";
		} else {
			$hash_info{$values[0]} = $values[1];
		}
	}

	return \%hash_info;
}
=head2 parse_format
    About   : Parse field format
    Usage   : my $x = parse_format($format);
    Args    : field to parse.
=cut
sub parse_format {
	my ($format,$hash_header,@formats) = @_;

	my @order = split(":",$format);

	my %hash_formats;
	my $j = 0;
	foreach my $value (@formats) {
		$j++;
		my @elements = split(":",$value);
		foreach my $i (0..scalar(@order)-1) {
			if ($elements[$i]) {
				$hash_formats{$hash_header->{'header'}[8+$j]}{$order[$i]} = $elements[$i] ;
			} else {
				$hash_formats{$hash_header->{'header'}[8+$j]}{$order[$i]} = "." ;
			}
		}
	}

	return \%hash_formats;
}
=head2 parse_meta
    About   : Parse meta information
    Usage   : my $x = parse_meta(line);
    Args    : line to parse
=cut
sub parse_meta {
	my $line=shift(@_);

	my %hash_meta;
	if ($line =~ /##([a-zA-Z0-9\-\.^,]*)=([a-zA-Z0-9\-\.^,]*)$/) {
		$hash_meta{"meta-informations"}{$1} = $2;
	} elsif($line =~ /##([^=]*)=<ID=([^,]*)(,Number=)*([^,]*)(,Type=)*([^,]*),Description="([^"]*).*/) {
		if($7) {$hash_meta{"meta-informations"}{$1}{$2}{"Description"} = $7};
		if($6) {$hash_meta{"meta-informations"}{$1}{$2}{"Type"} = $6};
		if($4) {$hash_meta{"meta-informations"}{$1}{$2}{"Number"} = $4};
		$hash_meta{"meta-informations"}{$1}{$2}{"line"} = $_;
	}  else {
		substr($_, 0, 1) = "";
		my @header = split("\t",$_);
		$hash_meta{"header"} = \@header;
	}
	return %hash_meta;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sub vcf2tab {
	my ($self,$vcf,$output) = @_;
	#print Dumper \$vcf;

	my $line = "";
	my @direct;
	my @infos;
	my @formats;

	foreach my $i (0..7-1) {
		$line .= $vcf->{"header"}[$i]."\t";
		push(@direct, $vcf->{"header"}[$i]);
	}
	foreach my $info (sort keys $vcf->{'meta-informations'}->{'INFO'}) {
		$line .= $info."\t";
		push(@infos, $info);
	}

	foreach my $format (sort keys $vcf->{'meta-informations'}->{'FORMAT'}) {
		push(@formats,$format);
	}

	foreach my $i (9..scalar(@{$vcf->{'header'}})-1) {
		foreach my $format (@formats) {
			$line .= $vcf->{'header'}[$i]."_".$format."\t";
		}
	}
	chop($line);
	print $line."\n";

	foreach my $variant (sort keys $vcf->{"variants"}) {
		$line = "";
		my $test = $vcf->{"variants"}->{$variant};
		foreach my $val (@direct){
			#print $val;
			$line .= ${$test}->{$val}."\t";
		}
		foreach my $info (@infos) {
			if (exists (${$test}->{"INFO"}->{$info})) {
				$line .= ${$test}->{"INFO"}->{$info}."\t";
			} else {
				$line .= ".\t";
			}
			#print Dumper %${$test}->{"INFO"};
		}

		foreach my $i (9..scalar(@{$vcf->{'header'}})-1) {
			my $sample = $vcf->{'header'}[$i];
			foreach my $format (@formats) {
				if (exists (${$test}->{"FORMAT"}->{$sample}->{$format})) {
					#print Dumper %${$test}->{"FORMAT"}->{$sample}->{$format};
					$line .= ${$test}->{"FORMAT"}->{$sample}->{$format}."\t";
				} else {
					$line .= ".\t";
				}
			}
		}
		print $line."\n";
	}

}

1;
