#!/usr/local/bin/perl -w
# ContigLengthTable.pl
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set use
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

print STDERR "Correcting abundance for reference length.\n";

# Set Variables
## Options
my $opt_help;
my $input;
my $lengthfile;
my $output;

## Files
my $IN;
my $LENGTH;
my $OUT;

## Misc
my %contighash;

# Set options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	'l|lengthfile=s' => \$lengthfile,
	'o|output=s' => \$output
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open($IN, "<", "$input") || die "Unable to read $input: $!";
open($LENGTH, "<", "$lengthfile") || die "Unable to read $lengthfile: $!";
open($OUT, ">", "$output") || die "Unable to write to $output: $!";

# Save length information into a hash
while (my $line = <$LENGTH>) {
	chomp $line;
	my $contigname = (split /\t/, $line)[0];
	my $contiglength = (split /\t/, $line)[1];
	$contighash{$contigname} = $contiglength;
}

# Now run the filter
while (my $line = <$IN>) {
	chomp $line;
	my $contigname = (split /\t/, $line)[0];
	my $contighits = (split /\t/, $line)[1];
	my $sampleid = (split /\t/, $line)[2];
	my $contiglengthcall = $contighash{$contigname};
	my $correctedabund = 10^6 * $contighits / $contiglengthcall;
	print $OUT "$contigname\t$contighits\t$contiglengthcall\t$correctedabund\t$sampleid\n";
}

close($IN);
close($LENGTH);
close($OUT);

print STDERR "Finished length corrections.\n";
