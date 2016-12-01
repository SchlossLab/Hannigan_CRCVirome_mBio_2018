#!usr/bin/perl
# ParseClusters.pl
# Geoffrey Hannigan
# Patrick Schloss Lab
# University of Michigan

# NOTE: Im going to just crank this out without a hash but
# for speed I will want to write it with one.

# Set use
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Set variables
my $IN;
my $OUT;
my $opt_help;
my $input;
my $output;
my $cluster;
my $uniqueid;
my $acc;
my $orgname;
my $protid;
my $percentid;
my $percentidname;

# Set the options
GetOptions(
    'h|help' => \$opt_help,
    'i|input=s' => \$input,
    'o|output=s' => \$output
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

# Open files
open($IN, "<", "$input") || die "Unable to read in $input: $!";
open($OUT, ">", "$output") || die "Unable to write to $output: $!";

while (my $line = <$IN>) {
	chomp $line;
	if ($line =~ /^>(.+)$/) {
		$cluster = $1;
		# Replace spaces with underscores
		$cluster =~ s/\s+/_/g;
	} elsif ($line !~ /^>.+/) {
		$line =~ s/\s+/\t/g;
		my $string = (split /\t/, $line)[2];
		print "$string\n";
		$uniqueid = $1 if ($string =~ />(\d+)\|sp\|.*\|(\S+)\;_(\S+)_\|\|_(\S+).../);
		print "$uniqueid\n";
		$acc = $2 if ($string =~ />(\d+)\|sp\|.*\|(\S+)\;_(\S+)_\|\|_(\S+).../);
		$orgname = $3 if ($string =~ />(\d+)\|sp\|.*\|(\S+)\;_(\S+)_\|\|_(\S+).../);
		$protid = $4 if ($string =~ />(\d+)\|sp\|.*\|(\S+)\;_(\S+)_\|\|_(\S+).../);
		$percentid = (split /\t/, $line)[3];
		$percentidname = "ref" if ($percentid =~ /\*/);
		$percentidname = (split /\t/, $line)[4] if ($percentid =~ /at/);
		print $OUT "$uniqueid\t$cluster\t$acc\t$orgname\t$protid\t$percentidname\n";
	}
}
