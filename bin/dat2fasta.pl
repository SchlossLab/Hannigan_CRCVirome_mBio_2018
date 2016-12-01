#! usr/bin/perl
# dat2fasta.pl
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
# And because I like timing myself
my $start_run = time();

# Set variables
my $IN;
my $OUT;
my $opt_help;
my $input;
my $output;
my $flag = 0;
my $formatVar;
my $sequence;
my $prot = '';
my $line;
my $SaveVariable;
my $genome;
my $secondVariable;
my $genomeSaver;
my $genename;
my $uniqcounter = 0;

# Set the options
GetOptions(
    'h|help' => \$opt_help,
    'd|datInput=s' => \$input,
    'f|fastaOutput=s' => \$output,
    'p|prot' => \$prot,
    'g|genome' => \$genome
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

# Open files
open($IN, "<", "$input") || die "Unable to read in $input: $!";
open($OUT, ">", "$output") || die "Unable to write to $output: $!";

# First option if the protein flag is used
if ($prot) {
    while ($line = <$IN>) {
	chomp $line;
        # Start the script by resetting the flag for each iteraction
        # within the file
        if ($flag == 0 && $line =~ /^ID\s+(\S+)\s/) {
            $SaveVariable = $1;
            $SaveVariable =~ s/\s/_/g;
            $flag = 0;
            $formatVar = 0;
    		$sequence = 0;
            next;
        } elsif ($flag == 0 && $line =~ /^AC\s+(\w.+)\;$/) {
            $secondVariable = "sp\|$1\|$SaveVariable";
            undef $SaveVariable;
            $flag = 0;
            next;
        } elsif ($flag == 0 && $line =~ /^OS\s+(\w.+$)/) {
            print $OUT ">$secondVariable $1\n" unless ($genome);
            $genomeSaver = "$secondVariable $1" if ($genome);
            $flag = 1;
            undef $secondVariable;
            next;
        # If this is not a genome protein file set
        } elsif ($flag == 1 && $line =~ /^\s+([A-Z\s]+[A-Z\s])$/ && !$genome) {
            $formatVar = $1;
            $formatVar =~ s/\s//g;
            $sequence = $formatVar;
            $flag = 2;
        } elsif ($flag == 2 && $line =~ /^\s+([A-Z\s]+[A-Z\s])$/ && !$genome) {
            $formatVar = $1;
            $formatVar =~ s/\s//g;
            $sequence = $sequence.$formatVar;
        # If these genes are from a genome file
        } elsif ($flag == 1 && $line =~ /^FT\s+\/product\="(.+)"$/ && $genome) {
            $genename = $1;
            $genename =~ s/\s/_/g;
            my $seqnumber = sprintf("%012d", $uniqcounter);
            print $OUT ">$seqnumber|$genomeSaver || $genename\n";
            ++$uniqcounter;
            $flag = 2;
        } elsif ($flag == 2 && $line =~ /^FT\s+\/translation\=\"([A-Z]+)$/ && $genome) {
            $formatVar = $1;
            $formatVar =~ s/\s//g;
            $sequence = $formatVar;
            $flag = 3;
        } elsif ($flag == 3 && $line =~ /^FT\s+([A-Z]+)/ && $genome) {
            $formatVar = $1;
            $formatVar =~ s/\s//g;
            $sequence = $sequence.$formatVar;
            print $OUT "$sequence\n" if ($line =~ /\"$/);
            $flag = 1 if ($line =~ /\"$/);
        } elsif ($flag == 1 && $line =~ /^SQ/ && $genome) {
            $flag = 2;
        } elsif ($flag == 2 && $line =~ /^\/\//) {
            print $OUT "$sequence\n";
            undef $genomeSaver;
            $flag = 0;
        }
    }
} else {
    while ($line = <$IN>) {
        chomp $line;
        # Start the script by resetting the flag for each iteraction
        # within the file
        if ($line =~ /^ID\s/) {
            print STDERR "Resetting counter...\n";
            $flag = 0;
            $formatVar = 0;
            $sequence = 0;
            next;
        } elsif ($flag =~ 0 & $line =~ /^OS\s+(\w.+$)/) {
            print STDERR "Name is $1\n";
            print $OUT "\>$1\n";
            $flag = 1;
            next;
        } elsif ($flag =~ 1 && $line =~ /^\s+([agct\s]+[agct])\s+[0-9]+$/) {
            $formatVar = $1;
            $formatVar =~ s/\s//g;
            $sequence = $formatVar;
            $flag = 2;
        } elsif ($flag =~ 2 && $line =~ /^\s+([agct\s]+[agct])\s+[0-9]+$/) {
            $formatVar = $1;
            $formatVar =~ s/\s//g;
            $sequence = $sequence.$formatVar;
        } elsif ($flag =~ 2 && $line =~ /^\/\//) {
            print $OUT "$sequence\n";
        } else {
            next;
        }
    }
}

close($IN);
close($OUT);

# See how long it took
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "Conversion completed in $run_time seconds.\n";
