## calculate_abundance_from_sam.pl
## Hannigan, Meisel, et. al
## The Human Skin Virome and its Interactions with the Host Microbiome

#!/usr/local/bin/perl -w

use strict;
use warnings;

my $usage = "Usage: perl $0 <INFILE> <OUTFILE>";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

if($infile =~ /\.sam$/) {
        open(IN, "<$infile") || die "Unable to open $infile: $!";
}
else {
        open(IN, "samtools view $infile |") || die "Unable to open $infile: $!";
}

open(OUT, ">$outfile") || die "Unable write to $outfile: $!";

# count each contig
my %contig_count;
while(my $line = <IN>) {
        if($line =~ /^@[A-Z]{2}/) {
                next;
        }
        my ($contig) = (split(/\t/, $line))[2];
        $contig_count{$contig}++;
}

# output
print OUT "contig\tcount\n";
#foreach my $contig (keys %contig_count) {}
while(my ($contig, $count) = each %contig_count) {
        print OUT "$contig\t$count\n";
}

close(IN);
close(OUT);