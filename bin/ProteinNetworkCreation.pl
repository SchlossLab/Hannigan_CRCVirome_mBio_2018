#! /usr/bin/perl
# ProteinNetworkCreation.pl
# Geoffrey Hannigan
# Patrick Schloss Lab
# University of Michigan

# Set use
use strict;
use warnings;
# Use the neo4j module to facilitate interaction
use REST::Neo4p;
# For documentation and whatnot
use Getopt::Long;
use Pod::Usage;

# Startup the neo4j connection using default location
# Be sure to set username and password as neo4j
# User = 2nd value, PW = 3rd value
eval {
	REST::Neo4p->connect('http://localhost:7474/', "neo4j", "neo4j");
};
ref $@ ? $@->rethrow : die $@ if $@;

# Set variables
my $opt_help;
my $phage;
my $bacteria;
my $dat;
my $n1;
my $n2;
my $flag;
my $formatVar;
my $formname;

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'b|bacteria=s' => \$bacteria,
	'p|phage=s' => \$phage,
	'd|dat=s' => \$dat
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(my $BACTERIA, "<", "$bacteria") || die "Unable to read in $bacteria: $!";
my $bacterialinecount = 0;
$bacterialinecount++ while <$BACTERIA>;
print STDERR "Bacteria line count is $bacterialinecount\n";
open(my $PHAGE, "<", "$phage") || die "Unable to read in $phage: $!";
my $phagelinecount = 0;
$phagelinecount++ while <$PHAGE>;
print STDERR "Bacteria line count is $phagelinecount\n";
open(my $DAT, "<", "$dat") || die "Unable to read in $dat: $!";

# Reset the files
seek $BACTERIA, 0, 0;
seek $PHAGE, 0, 0;

sub AddNodes {
	my ($fileInput, $label, $linecount) = @_;
	my $progcounter = 0;
	my $progress = 0;
	while (my $line = <$fileInput>) {
		if ($progcounter % 1000 == 0) {
			$progress = 100 * $progcounter / $linecount;
			print STDERR "\rNodes Processed: $progress\%";
		}
		++$progcounter;
		chomp $line;
		$line =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
		my $uniqueid = (split /\t/, $line)[0];

		# Skip if it has already been added
		my @n11 = REST::Neo4p->get_nodes_by_label( $uniqueid );
		die "Oh no! You have duplicate unique nodes: $!" if (scalar(@n11) gt 1);
		next if (@n11);

		my $clusterid = (split /\t/, $line)[1];
		my $acc = (split /\t/, $line)[2];
		my $name = (split /\t/, $line)[3];
		$n1 = REST::Neo4p::Node->new( {UniqueID => $uniqueid} );
		$n1->set_property( {ClusterID => $clusterid} );
		$n1->set_property( {Acccession => $acc} );
		$n1->set_property( {DataType => "ReferenceGenes"} );
		$n1->set_labels($label, $name);
	}
}

# Run the subroutines
print STDERR "\nPROGRESS: Creating Phage Nodes.\n";
AddNodes(\*$PHAGE, "Phage", $phagelinecount);
print STDERR "\nPROGRESS: Creating Bacteria Nodes.\n";
AddNodes(\*$BACTERIA, "Bacteria", $bacterialinecount);
print STDERR "\nPROGRESS: Establishing Relationships.\n";

my @phagenodes;
my @bacterianodes;

while (my $line = <$DAT>) {
	chomp $line;
	if ($line =~ /^ID\s/) {
		$flag = 0;
		$n1 = 0;
		$n2 = 0;
		$formatVar = 0;
		$formname = 0;
		next;
	} elsif ($flag =~ 0 & $line =~ /^OS\s+(\w.+$)/) {
		# File really should already be without spaces though
		($formname = $1) =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
		$formname =~ s/_$//;
		@phagenodes = REST::Neo4p->get_nodes_by_label( $formname );
		my $PhageLength = length scalar(@phagenodes);
		$flag = 1;
		next;
	} elsif ($flag =~ 1 & $line =~ /host=\"(.+)\"/) {
		(my $FullName = $1) =~ s/\s/_/g;
		$FullName =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
		$FullName =~ s/_$//;
		my @bacterianodes = REST::Neo4p->get_nodes_by_label( $FullName );
		my $BacteriaLength = length scalar(@bacterianodes);
		foreach my $phagenode (@phagenodes) {
			foreach my $bacterianode (@bacterianodes) {
				$phagenode->relate_to($bacterianode, 'LinkedGenes')->set_property({Literature => "TRUE"});
				print STDOUT "Relationship created from $bacterianode to $phagenode\n";
			}
		}
	} else {
		next;
	}
}
