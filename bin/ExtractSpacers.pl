#!/usr/bin/perl
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan
# OVERVIEW: This is a simple script to extract the CRISPR
# spacers from piler-cr output into a fasta format.

# Set use
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
# And because I like timing myself
my $start_run = time();

# Set variables
my $opt_help;
my $IN;
my $OUT;
my $input;
my $output;
my $flag = 0;
my $name;
my $alterLine;
my $spacerSeq = 0;
my $SpaceCount = 0;
my $Header = 0;
my $location = 0;

# Set the options
GetOptions(
    'h|help' => \$opt_help,
    'i|datInput=s' => \$input,
    'o|fastaOutput=s' => \$output
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

# Open files
open($IN, "<", "$input") || die "Unable to read in $input: $!";
open($OUT, ">", "$output") || die "Unable to write to $output: $!";

# Parse the piler-cr output
foreach my $line (<$IN>) {
	chomp $line;
	if ($location =~ 0 & $flag =~ 0 & $line =~ /\>\S+/) {
		$name = $line;
		print STDERR "Getting spacers for $name.\n";
		$flag = 1;
		next;
	} elsif ($location =~ 0 & $flag =~ 1 & $line =~ /^\s+[0-9]+/) {
		$SpaceCount = $SpaceCount + 1;
		($alterLine = $line) =~ s/\s+/\t/g;
		$alterLine =~ s/^\t//;
		$spacerSeq = (split /\t/, $alterLine)[-1];
		print $OUT "$name"."_$SpaceCount\n$spacerSeq\n";
		next;
	} elsif ($location =~ 0 & $flag =~ 1 & $line =~ /^\===/) {
		if ($Header =~ 0) {
			$Header = 1;
			next;
		} elsif ($Header =~ 1) {
			# Reset when done with the array
			$Header = 0;
			$flag = 0;
			$name = 0;
			$SpaceCount = 0;
			$alterLine = 0;
		}
		next;
	} elsif ($line =~ /SUMMARY\sBY/) {
		$location = 1;
		next;
	} elsif ($line =~ /DETAIL\sREPORT/) {
		$location = 0;
	}
}

# Close files
close($IN);
close($OUT);

# See how long it took
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "Extracted spacers in $run_time seconds.\n";
