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

print STDERR "Creating sequence length table.\n";

# Set Variables
## Options
my $opt_help;
my $input;
my $output;

## Files
my $IN;
my $OUT;

## Misc
my $flag = 0;
my $LineLength;
my $lineName;

# Set options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	'o|output=s' => \$output,
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open($IN, "<", "$input") || die "Unable to read $input: $!";
open($OUT, ">", "$output") || die "Unable to write to $output: $!";

# Now run the filter
while (my $line = <$IN>) {
	chomp $line;
	if ($line =~ /^\>/ && $flag==0) {
		($lineName = $line) =~ s/^\>//;
		$flag = 1;
		next;
	} elsif ($line !~ /^\>/ && $flag==1) {
		$LineLength = length $line;
		print $OUT "$lineName\t$LineLength\n";
		$flag = 0;
		$lineName = '';
	} else {
		die "Error in the fasta format: $!";
	}
}

close($IN);
close($OUT);

print STDERR "Finished creating sequence length table.\n";
