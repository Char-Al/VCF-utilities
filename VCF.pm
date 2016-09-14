package VCF;

use Data::Dumper;
=head2 parse_line
    About   : Reads a VCF line and splits it into a hash.
    Usage   : my $x = parse_line($line);
    Args    : line to parse.
=cut
sub parse_line {
	my ($it,$line) = @_;
	## Control line
	if ( !$line ) { return undef; }

	my %hash_line;
	if($line =~ /^#/) {
		%hash_line = parse_meta($line);
		return (\%hash_line,"HEADER");
	} else {
		my $ncols = scalar(split("\t",$line));
		if($ncols > 8) {
			my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info, $format, @formats) = split("\t",$line);%hash_line = (
				"chr"		=> $chr,
				"pos"		=> $pos,
				"ID"		=> $ID,
				"ref"		=> $ref,
				"alt"		=> $alt,
				"qual"		=> $qual,
				"filter"	=> $filter,
				"info"		=> parse_info($info),
				"format"	=> parse_format($format, @formats)
			);
		} else {
			my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info) = split("\t",$line);%hash_line = (
				"chr"		=> $chr,
				"pos"		=> $pos,
				"ID"		=> $ID,
				"ref"		=> $ref,
				"alt"		=> $alt,
				"qual"		=> $qual,
				"filter"	=> $filter,
				"info"		=> parse_info($info)
			);
		}

			return (\%hash_line,"Var");

	}
	print Dumper \%hash_line;
}

=head2 parse_info
    About   : Parse field info
    Usage   : my $x = parse_info($info);
    Args    : field to parse.
=cut
sub parse_info {
	my $info = shift(@_);
	my @infos = split(";",$info);

	my %hash_info;

	foreach my $value (@infos) {
		my @values = split("=",$value);
		$hash_info{$values[0]} = $values[1];
	}

	return \%hash_info;
}
=head2 parse_format
    About   : Parse field format
    Usage   : my $x = parse_format($format);
    Args    : field to parse.
=cut
sub parse_format {
	my ($format,@formats) = @_;

	my @order = split(":",$format);

	my %hash_formats;
	my $j = 0;
	foreach my $value (@formats) {
		$j++;
		my @elements = split(":",$value);
		foreach my $i (0..scalar(@order)-1) {
			$hash_formats{"sample_".$j}{$order[$i]} = $elements[$i]
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
		$hash_meta{"meta"}{$1} = $2;
	} elsif($line =~ /##([^=]*)=<ID=([^,]*)(,Number=)*([^,]*)(,Type=)*([^,]*),Description="([^"]*).*/) {
		if($7) {$hash_meta{"meta"}{$1}{$2}{"Description"} = $7};
		if($6) {$hash_meta{"meta"}{$1}{$2}{"Type"} = $6};
		if($4) {$hash_meta{"meta"}{$1}{$2}{"Number"} = $4};
		$hash_meta{"meta-informations"}{$1}{$2}{"line"} = $_;
	}  else {
		$hash_meta{"header"} = $_;
	}
	return %hash_meta;
}

1;
