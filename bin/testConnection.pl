#!usr/bin/perl

# Set use
use strict;
use warnings;
use REST::Neo4p;

eval {
	REST::Neo4p->connect('http://127.0.0.1:7474', "neo4j", "neo4j");
};
ref $@ ? $@->rethrow : die $@ if $@;

my $n1;
my $n2;
my $r1;

$n1 = REST::Neo4p::Node->new( {name => 'Harry'} );
$n2 = REST::Neo4p::Node->new( {name => 'Sally'} );
$r1 = $n1->relate_to($n2, 'met');
$r1->set_property({ when => 'July' });

