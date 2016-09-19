package VCF;

use Data::Dumper;
use File::Basename;
use Hash::Merge qw(merge);

use feature "switch";
no warnings 'experimental';


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
		$j++;
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
	if ( !$line ) { return undef; }

	my %hash_line;
	if($line =~ /^#/) {
		%hash_line = parse_meta($line,$cpt);
		return (\%hash_line,"HEADER");
	} else {
		my $ncols = scalar(split("\t",$line));
		if($ncols > 8) {
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
	my ($line) = @_;
	my %hash_meta;
	if ($line =~ /##([a-zA-Z0-9\-\.^,]*)=([a-zA-Z0-9\-\.^,]*)$/) {
		$hash_meta{"meta-informations"}{$1} = $2;
	} elsif($line =~ /##([^=]*)=<ID=([^,]*)(,Number=)*([^,]*)(,Type=)*([^,]*),Description="([^"]*).*/) {
			if ($1 eq "INFO") {
				$i++;
			}
			if($7) {$hash_meta{"meta-informations"}{$1}{$2}{"Description"} = $7};
			if($6) {$hash_meta{"meta-informations"}{$1}{$2}{"Type"} = $6};
			if($4) {$hash_meta{"meta-informations"}{$1}{$2}{"Number"} = $4};
			$hash_meta{"meta-informations"}{$1}{$2}{"line"} = $_;
			$hash_meta{"meta-informations"}{$1}{$2}{"order"} = $i;
	}  elsif($line =~ /^#CHROM/) {
		substr($_, 0, 1) = "";
		my @header = split("\t",$_);
		$hash_meta{"header"} = \@header;
	}
	return %hash_meta;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

=head2 vcf2tab
    About   : Convert a vcf (hash) into a tab separate file
    Usage   : my parse_meta($hash_vcf, $output);
    Args    :
		- hash of a vcf
		- output_file
=cut
sub vcf2tab {
	my ($self, $header, $meta, $variants, $output) = @_;

	open(FILE,">$output") or die "Cannot open file : '$output' $!";

	my $line = "";
	my @direct;
	my @infos;
	my @formats;

	foreach my $i (0..7-1) {
		$line .= $$header[$i]."\t";
		push(@direct, $$header[$i]);
	}
	foreach my $info (sort keys $meta->{'INFO'}) {
		$infos[$meta->{"INFO"}->{$info}->{'order'}-1] = $info;
	}
	$line .= join("\t",@infos)."\t";

	foreach my $format (sort keys $meta->{'FORMAT'}) {
		push(@formats,$format);
	}

	foreach my $i (9..scalar(@{$header})-1) {
		foreach my $format (@formats) {
			$line .= $$header[$i]."_".$format."\t";
		}
	}
	chop($line);
	print FILE $line."\n";

	foreach my $variant (sort keys $variants) {
		$line = "";
		my $hash_variant = $variants->{$variant};
		foreach my $val (@direct){
			$line .= ${$hash_variant}->{$val}."\t";
		}
		foreach my $info (@infos) {

			if (exists (${$hash_variant}->{"INFO"}->{$info})) {
				$line .= ${$hash_variant}->{"INFO"}->{$info}."\t";
			} else {
				$line .= ".\t";
			}
		}

		foreach my $i (9..scalar(@{$header})-1) {
			my $sample = $$header[$i];
			foreach my $format (@formats) {
				if (exists (${$hash_variant}->{"FORMAT"}->{$sample}->{$format})) {
					$line .= ${$hash_variant}->{"FORMAT"}->{$sample}->{$format}."\t";
				} else {
					$line .= ".\t";
				}
			}
		}
		print FILE $line."\n";
	}
	close(FILE);

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

=head2 vcfCondFilter
    About   :
    Usage   :
    Args    :
=cut
sub vcfCondFilter {
	my ($self,$vcf,$cond) = @_;
	my %filter;
	foreach my $variant (sort keys $vcf->{"variants"}) {
		my $hash_variant = $vcf->{"variants"}->{$variant};
		foreach my $key (sort keys($cond)) {
			my %boolean_rules = (
				1 => 0,
				0 => 0
			);

			foreach my $rules (sort keys $cond->{$key}->{"rules"}) {
				$boolean_rules{check_filter($cond->{$key}->{"rules"}->{$rules}, $vcf->{"variants"}->{$variant})}++;
			}

			my $verify = check_condition($cond->{$key}->{"cond"},$boolean_rules{1},$boolean_rules{0}) ;

			my $check ="";

			for ($cond->{$key}->{"step"}) {
				when("next")	{ if ($verify) { next; }									else { $check = "next_line"}	}
				when("get")		{ if ($verify) { $filter{$variant} = $vcf->{"variants"}->{$variant} ; $check = "next_line"; }	else { next ; }					}
				when("stop")	{ if ($verify) { $filter{$variant} = $vcf->{"variants"}->{$variant} ; $check = "next_line"; }	else { $check = "next_line"}	}
			}

			if ($check eq "next_line") {
				last;
			}
		}
	}
	return \%filter;
}

=head2 parseConfFile
    About   :
    Usage   :
    Args    :
=cut
sub parseConfFile {
	my ($self, $config) = @_;

	my %hash_filter;
	open(CONF,"$config") or die ("cannot open $config $!");
	my $i = 1;
	while(<CONF>){
		chomp();
		my @line = split("\t");

		if($line[1] eq "INFO" or $line[1] eq "FORMAT") {
			$hash_filter{$line[0]}{"rules"}{$i}{"target"}	= $line[1];
			$hash_filter{$line[0]}{"rules"}{$i}{"id"}		= $line[2];
			$hash_filter{$line[0]}{"rules"}{$i}{"is"}		= $line[3];
			$hash_filter{$line[0]}{"rules"}{$i}{"for"}		= $line[4];
			$hash_filter{$line[0]}{"rules"}{$i}{"val"}		= $line[5];

			$hash_filter{$line[0]}{"cond"} = $line[6];
			$hash_filter{$line[0]}{"step"} = $line[7];
		} else {
			$hash_filter{$line[0]}{"rules"}{$i}{"target"}	= "BASE";
			$hash_filter{$line[0]}{"rules"}{$i}{"id"}		= $line[1];
			$hash_filter{$line[0]}{"rules"}{$i}{"is"}		= $line[2];
			$hash_filter{$line[0]}{"rules"}{$i}{"for"}		= $line[3];
			$hash_filter{$line[0]}{"rules"}{$i}{"val"}		= $line[4];

			$hash_filter{$line[0]}{"cond"} = $line[5];
			$hash_filter{$line[0]}{"step"} = $line[6];
		}

		$i++;
	}
	return \%hash_filter;
}

=head2 check_filter
    About   :
    Usage   :
    Args    :
=cut
sub check_filter {
	my ($hash_rule, $hash_line) = @_;
	my $val = "";

	if($hash_rule->{"target"} eq "BASE"){
		$val = ${$hash_line}->{$hash_rule->{"id"}};
	} else {
		$val = ${$hash_line}->{$hash_rule->{"target"}}->{$hash_rule->{"id"}};
	}

	for($hash_rule->{"is"}) {
		when("Num") {
			for ($hash_rule->{"for"}) {
				when("==")	{ $val == $hash_rule->{"val"} ? return 1 : return 0 ; }
				when(">=")	{ $val >= $hash_rule->{"val"} ? return 1 : return 0 ; }
				when("<=")	{ $val <= $hash_rule->{"val"} ? return 1 : return 0 ; }
				when(">")	{ $val >  $hash_rule->{"val"} ? return 1 : return 0 ; }
				when("<")	{ $val <  $hash_rule->{"val"} ? return 1 : return 0 ; }
				when("!=")	{ $val != $hash_rule->{"val"} ? return 1 : return 0 ; }
				default		{return 0}
			}
		}
		when("Str") {
			for ($hash_rule->{"for"}) {
				when("==")	{ $val eq $hash_rule->{"val"} ? return 1 : return 0 ; }
				# when(">=")	{ $val >= $hash_rule->{"val"} ? return 1 : return 0 ; }
				# when("<=")	{ $val <= $hash_rule->{"val"} ? return 1 : return 0 ; }
				# when(">")	{ $val >  $hash_rule->{"val"} ? return 1 : return 0 ; }
				# when("<")	{ $val <  $hash_rule->{"val"} ? return 1 : return 0 ; }
				when("!=")	{ $val ne $hash_rule->{"val"} ? return 1 : return 0 ; }
				default		{return 0}
			}
		}
	}
}

=head2 check_condition
    About   :
    Usage   :
    Args    :
=cut
sub check_condition {
	my ($cond, $true, $false) = @_;

	for($cond ) {
		when("or")	{ ($true >= 1) ? return 1 : return 0 ; }
		when("and")	{ ($false == 0) ? return 1 : return 0 ; }
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

=head2 compareVcf
    About   :
    Usage   :
    Args    :
=cut
sub compareVcf {
	my ($self,@vcfs) = @_;

	compareSimple(@vcfs);
}

=head2 compareSimple
    About   :
    Usage   :
    Args    :
=cut
sub compareSimple {
	my ($self,@vcfs) = @_;

	compareSimple(@vcfs);
}

=head2 speed_comp
    About   : Reads a VCF line and splits it into a hash.
    Usage   : my $x = parse_line($line);
    Args    : line to parse.
=cut
sub speed_comp {
	my ($self,$file_name) = @_;

	print STDERR "Treat file : '$file_name'\n";
	open(VCF,$file_name) or die ("Error !\n");
	$file_name = basename($file_name);
	my @all_lines;
	my $hash_variants = {};
	my $hash_header = {};

	while(<VCF>) {
		chomp($_);
		my ($hash_line,$type) = parse_line($_,$hash_header);
		if($type eq "HEADER") {
			$hash_header = merge($hash_header,$hash_line);
		} else {
			my @alts = split(",",$hash_line->{"ALT"});
			foreach my $alt (@alts) {
				my $ID = $hash_line->{"CHROM"}."_".$hash_line->{"POS"}."_".$hash_line->{"REF"}."_".$alt;
				$hash_variants->{$ID} = $_;
			}
		}
	}
	close(VCF);
	return $hash_variants;
}
=head2 speed_read
    About   : Reads a VCF line and splits it into a hash.
    Usage   : my $x = parse_line($line);
    Args    : line to parse.
=cut
sub speed_read {
	my ($self,$file_name) = @_;

	print STDERR "Treat file : '$file_name'\n";
	open(VCF,$file_name) or die ("Error !\n");
	$file_name = basename($file_name);
	my @all_lines;
	my $hash_variants = {};
	my $hash_header = {};

	while(<VCF>) {
		chomp($_);
		my ($hash_line,$type) = parse_line($_,$hash_header);
		if($type eq "HEADER") {
			$hash_header = merge($hash_header,$hash_line);
		} else {
			my @alts = split(",",$hash_line->{"ALT"});
			foreach my $alt (@alts) {
				my $ID = $hash_line->{"CHROM"}."_".$hash_line->{"POS"}."_".$hash_line->{"REF"}."_".$alt;
				$hash_variants->{$ID} = $_;
			}
		}
	}
	close(VCF);
	return $hash_variants;
}





1;
