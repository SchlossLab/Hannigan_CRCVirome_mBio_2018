#!usr/bin/perl
# BenchmarkDatabaseCreation.pl
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
# And because I like timing myself
my $start_run = time();

# Set variables
my $opt_help;
my $output;
my $input;
my $flag = 0;
my $n1;
my $n2;
my $sequence;
my $formatVar;
my $FullName;
my $Genus;
my $Species;
my $formname = 0;
my $Spacer;
my $crispr;
my $PhageTarget;
my $array1;
my $array2;
my $uniprot;
my $relate;
my $phage;
my $bacteria;
my $phageForm;
my $bacteriaForm;
my $blast;
my $pfam;
my $SpacerForm;
my $PhageTargetForm;
my $relnHit;
my $score;
my $scorenum;
my $blastx;
my $validation;

# Startup the neo4j connection using default location
# Be sure to set username and password as neo4j
# User = 2nd value, PW = 3rd value
eval {
    REST::Neo4p->connect('http://localhost:7474/',"neo4j","neo4j");
};
ref $@ ? $@->rethrow : die $@ if $@;

# Set the options
GetOptions(
    'h|help' => \$opt_help,
    'i|input=s' => \$input,
    'c|crispr=s' => \$crispr,
    'b|blast=s' => \$blast,
    'p|pfam=s' => \$pfam,
    'x|blastx=s' => \$blastx,
    'v|validation' => \$validation
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

# Open files
open(my $IN, "<", "$input") || die "Unable to read in $input: $!";
open(my $CRISPR, "<", "$crispr") || die "Unable to read in $crispr: $!";
open(my $BLAST, "<", "$blast") || die "Unable to read in $blast: $!";
open(my $PFAM, "<", "$pfam") || die "Unable to read in $pfam: $!";
open(my $BLASTX, "<", "$blastx") || die "Unable to read in $blastx: $!";

sub AddGenericFile {
    # Import file handle
    my ($fileInput, $label, $score) = @_;
    # print "Score is $score\n";
    foreach my $line (<$fileInput>) {
        # print $line;
        chomp $line;
        $Spacer = (split /\t/, $line)[0];
        $PhageTarget = (split /\t/, $line)[1];
        $scorenum = (split /\t/, $line)[2] unless ($score eq "FALSE");
        # Remove illegal characters
        ($SpacerForm = $Spacer) =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
        # print "Bacteria Form: $SpacerForm\n";
        ($PhageTargetForm = $PhageTarget) =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
        # print "Phage Form: $PhageTargetForm\n";
        my @n11 = REST::Neo4p->get_nodes_by_label( $PhageTargetForm );
        # print scalar(@n11)."\n";
        my @n12 = REST::Neo4p->get_nodes_by_label( $SpacerForm );
        # print scalar(@n12)."\n";
        # Create new phage target node if it does not exist
        unless (@n11) {
            $formname = $PhageTargetForm;
            # print STDERR "New CRISPR phage target is $formname\n";
            $n1 = REST::Neo4p::Node->new( {Name => $formname} );
            # print STDERR "$n1\n";
            $n1->set_property( {Organism => 'Phage'} );
            $n1->set_labels('Phage',$formname);
        }
        unless (@n12) {
            ($FullName = $Spacer) =~ s/\s/_/g;
            # print STDERR "New spacer host is $FullName\n";
            $Genus = (split /_/, $FullName)[0];
            $Species = $Genus."_".(split /_/, $FullName)[1];
            # Remove all non-standard characters from the variable names
            # now that underscore delimiter is not being used
            $FullName =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
            $Genus =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
            $Species =~ s/[^A-Z^a-z^0-9^\t]+/_/g;
            # print STDERR "CRISPR host genus is $Genus\n";
            # print STDERR "CRISPR host species is $Species\n";
            $n2 = REST::Neo4p::Node->new( {Name => $FullName} );
            $n2->set_property( {Genus => $Genus} );
            $n2->set_property( {Species => $Species} );
            $n2->set_property( {Organism => 'Bacterial_Host'} );
            $n2->set_labels('Bacterial_Host',$FullName);
        }
        
        # Then get the newly created nodes as arrays
        @n11 = REST::Neo4p->get_nodes_by_label( $PhageTargetForm );
        @n12 = REST::Neo4p->get_nodes_by_label( $SpacerForm );
        # print scalar(@n11)."\n";
        # print scalar(@n12)."\n";
        # Ensure there are no duplicated nodes
        die "You dont have only one phage node ID: $!" unless (scalar(@n11) eq 1);
        die "You dont have only one duplicate bacteria node ID: $!" unless (scalar(@n12) eq 1);
        # Get nodes into scalar variables
        $array1 = pop @n11;
        $array2 = pop @n12;

        # Determine whether relationship exists
        my @phage_reln = $array1->get_outgoing_relationships();
        $flag = 0;
        foreach my $RelnItr (@phage_reln) {
            my $PhageNode = $RelnItr->start_node;
            my $BacteriaNode = $RelnItr->end_node;
            if ($PhageNode eq $array1 && $BacteriaNode eq $array2) {
                $flag = 1;
                $relnHit = $RelnItr;
                last;
            } else {
                next;
            }
        }
        if ($score eq "FALSE") {
            if ($flag eq 0) {
                # This means I need to create a new relationship
                $array1->relate_to($array2, 'Infects')->set_property({$label => "TRUE"});
            } else {
                $relnHit->set_property({ $label => "TRUE" });
            }
        } else {
            if ($flag eq 0) {
                # This means I need to create a new relationship
                $array1->relate_to($array2, 'Infects')->set_property({$label => $scorenum});
            } else {
                $relnHit->set_property({ $label => $scorenum });
            }
        }
        
    }
}

print STDERR "\nRunning Data As Validation Dataset\n" if defined $validation;
print STDERR "\n\n\nProgress: Adding Experimentally Validated Interactions\n" if defined $validation;
AddGenericFile(\*$IN, "Interaction", "TRUE") if defined $validation;

print STDERR "\n\n\nProgress: Adding Predicted CRISPR Interactions\n";
AddGenericFile(\*$CRISPR, "CRISPR", "TRUE");

print STDERR "\n\n\nProgress: Adding Predicted BLAST Interactions\n";
AddGenericFile(\*$BLAST, "BLAST", "TRUE");

print STDERR "\n\n\nProgress: Adding Predicted Pfam Interactions\n";
AddGenericFile(\*$PFAM, "PFAM", "TRUE");

print STDERR "\n\n\nProgress: Adding Blastx comparisons\n";
AddGenericFile(\*$BLASTX, "BLASTX", "TRUE");


# See how long it took
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "Processed the file in $run_time seconds.\n";
