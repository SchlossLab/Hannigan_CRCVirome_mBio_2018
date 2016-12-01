#!usr/bin/perl
# Metadata2graph.pl
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
	REST::Neo4p->connect('http://127.0.0.1:7474', "neo4j", "neo4j");
};
ref $@ ? $@->rethrow : die $@ if $@;

my $opt_help;
my $samples;
my $metadata;

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	's|samples=s' => \$samples,
	'm|metadata=s' => \$metadata
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(my $SAMPLES, "<", "$samples") || die "Unable to read in $samples: $!";
open(my $META, "<", "$metadata") || die "Unable to read in $metadata: $!";

my $phageid = 0;
my $abund = 0;
my $sampleid = 0;
my $n1;

foreach my $line (<$SAMPLES>) {
	chomp $line;
	$phageid = (split /\t/, $line)[0];
	print "Organism ID is $phageid\n";
	$abund = (split /\t/, $line)[2];
	$sampleid = (split /\t/, $line)[1];
	print "Sample ID is $sampleid\n";

	my @n11 = REST::Neo4p->get_nodes_by_label( $sampleid );
	my @n12 = REST::Neo4p->get_nodes_by_label( $phageid );
	# print scalar(@n12)."\n";

	# Ensure there are no duplicated nodes
    die "You have duplicate sample node IDs: $!" if (scalar(@n11) gt 1);
    die "You have duplicate phage node IDs: $!" if (scalar(@n12) gt 1);
    next if (scalar(@n12) eq 0);

	unless (@n11) {
		$n1 = REST::Neo4p::Node->new( {Name => $sampleid} );
		$n1->set_property( {Organism => 'SampleID'} );
		$n1->set_labels('SampleID',$sampleid);
	}

	@n11 = REST::Neo4p->get_nodes_by_label( $sampleid );
	@n12 = REST::Neo4p->get_nodes_by_label( $phageid );
	print scalar(@n11)."\n";
	print scalar(@n12)."\n";

	my $array1 = pop @n11;
    my $array2 = pop @n12;

	# Ensure there are no duplicated nodes
    die "You have duplicate sample node IDs: $!" if (scalar(@n11) gt 1);
    die "You have duplicate phage node IDs: $!" if (scalar(@n12) gt 1);

	$array1->relate_to($array2, 'Sampled')->set_property({Abundance => $abund});
}

my $disease;
my $studyid;
my $seconddisease;
my $mdatype;
my $bodylocation;
my $purificationtype;
my $location;
my $hosttype;
my $platform;

foreach my $line (<$META>) {
	chomp $line;
	my @linearray = split /\t/, $line;
	$studyid = $linearray[0];
	print STDERR "Study ID is $studyid\n";
	$sampleid = $linearray[2];
	$platform = $linearray[4];
	$disease = $linearray[5];
	print STDERR "Disease name is $disease\n";
	$mdatype = $linearray[7];
	$bodylocation = $linearray[8];
	$purificationtype = $linearray[9];
	$location = $linearray[10];
	$hosttype = $linearray[11];
	# Skip the header
	next if ($studyid eq "SRA_Study_s");
	print "Sample ID is $sampleid\n";

	# Get existing sample nodes
	my @n11 = REST::Neo4p->get_nodes_by_label( $sampleid );
	# Get existing disease nodes
	my @n12 = REST::Neo4p->get_nodes_by_label( $disease );
	# Get existing study nodes
	my @n13 = REST::Neo4p->get_nodes_by_label( $studyid );

	my $existingnodes = scalar(@n13);
	my $samplenodes = scalar(@n11);

	print STDERR "There are $existingnodes study ID nodes with this name.\n";
	print STDERR "There are $samplenodes sample ID nodes with this name.\n";

	# Ensure there are no duplicated nodes
    die "You have duplicate sample node IDs: $!" if (scalar(@n11) gt 1);
    die "You have duplicate disease node IDs: $!" if (scalar(@n12) gt 1);
    die "You have duplicate study node IDs: $!" if (scalar(@n13) gt 1);
    next if (scalar(@n11) eq 0);

    # Build nodes if the do not yet exist
    unless (@n12) {
		$n1 = REST::Neo4p::Node->new( {Name => $disease} );
		$n1->set_property( {Organism => 'Disease'} );
		$n1->set_labels('Disease',$disease);
	}
	unless (@n13) {
		print STDERR "Creating this StudyID node.\n";
		$n1 = REST::Neo4p::Node->new( {Name => $studyid} );
		$n1->set_property( {Organism => 'StudyID'} );
		$n1->set_labels('StudyID',$studyid);
	}

	@n11 = REST::Neo4p->get_nodes_by_label( $sampleid );
	@n12 = REST::Neo4p->get_nodes_by_label( $disease );
	@n13 = REST::Neo4p->get_nodes_by_label( $studyid );

	# Made array 1 the sample ID
	my $array1 = pop @n11;
	# Make array 2 the disease
    my $array2 = pop @n12;
    # Make array 3 the study
    my $array3 = pop @n13;

	# Ensure there are no duplicated nodes
    die "You have duplicate sample node IDs: $!" if (scalar(@n11) gt 1);
    die "You have duplicate phage node IDs: $!" if (scalar(@n12) gt 1);
    die "You have duplicate study node IDs: $!" if (scalar(@n13) gt 1);

    ###########################
    # Set the data into nodes #
    ###########################
    # Relate disease to sample
	$array2->relate_to($array1, 'Diseased')->set_property({Disease => "TRUE"});
	# Relate study to sample
	$array3->relate_to($array1, 'IncludedInStudy')->set_property({Study => "TRUE"});
	# Relate study to disease as well
	$array2->relate_to($array3, 'IncludedInStudy')->set_property({StudyDisease => "TRUE"});
	# Set sample properties
	# I think these could also be nodes but for now I am going to use them as edges
	$array1->set_property( {Platform => $platform} );
	$array1->set_property( {MDAtype => $mdatype} );
	$array1->set_property( {BodyLocation => $bodylocation} );
	$array1->set_property( {PurificationType => $purificationtype} );
	$array1->set_property( {Location => $location} );
	$array1->set_property( {Host => $hosttype} );
}
