#!usr/bin/perl
# mergeClusters.pl
# Geoffrey Hannigan
# Patrick Schloss Lab
# University of Michigan

# Set use
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Set variables
my $IN;
my $OUT;
my $REF;
my $opt_help;
my $input;
my $reference;
my $output;

# Set the options
GetOptions(
    'h|help' => \$opt_help,
    'i|input=s' => \$input,
    'r|reference' => \$reference,
    'o|output=s' => \$output
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

# Open files
open($IN, "<", "$input") || die "Unable to read in $input: $!";
open($REF, "<", "$reference") || die "Unable to read in $reference: $!";
open($OUT, ">", "$output") || die "Unable to write to $output: $!";


