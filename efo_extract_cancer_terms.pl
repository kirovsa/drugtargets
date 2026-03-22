#!/usr/bin/env perl
use strict;
use warnings;

use Bio::OntologyIO;
use Getopt::Long qw(GetOptions);

# Extract cancer-related terms from the EFO (Experimental Factor Ontology) OBO file.
#
# Reads an OBO file, finds the top-level "cancer" term (EFO:0000311), and traverses
# all descendant branches via is_a and other relationship edges using breadth-first search.
#
# Output: a tab-delimited file with two columns:
#   efo_id      - EFO identifier in EFO:XXXXXXX format
#   description - human-readable term label

my $input   = 'efo.obo';
my $output  = 'efo_cancer_terms.tsv';
my $root_id = 'EFO:0000311';

GetOptions(
  'input=s'  => \$input,
  'output=s' => \$output,
  'root=s'   => \$root_id,
) or die "Usage: $0 --input efo.obo --output efo_cancer_terms.tsv --root EFO:0000311\n";

my $io = Bio::OntologyIO->new(
  -format => 'obo',
  -file   => $input,
);

my $ont = $io->next_ontology();
die "Failed to read ontology from $input\n" unless $ont;

my $root = $ont->find_term(-identifier => $root_id);
die "Could not find root term $root_id in $input\n" unless $root;

my %visited;
my @queue = ($root);

while (@queue) {
  my $term = shift @queue;
  my $id = $term->identifier;
  next unless defined $id && length $id;
  next if $visited{$id}++;

  my @children;

  # Prefer the ontology helper if available (covers is_a; may cover more depending on implementation)
  eval {
    push @children, $ont->get_child_terms($term);
    1;
  };

  # Also try relationship objects if supported to include non-is_a edges.
  eval {
    my @rels = $ont->get_relationships(-term => $term);
    for my $rel (@rels) {
      my $subj = eval { $rel->subject_term } || undef;
      my $obj  = eval { $rel->object_term  } || undef;

      # Many OBO relationships are stored as child(subject) -> parent(object)
      # If current term is the parent (object), then the subject is a child.
      if ($obj && $obj->identifier && $obj->identifier eq $id && $subj) {
        push @children, $subj;
      }
    }
    1;
  };

  my %seen_child;
  for my $ch (@children) {
    next unless $ch;
    my $cid = $ch->identifier;
    next unless defined $cid && length $cid;
    next if $visited{$cid};
    next if $seen_child{$cid}++;
    push @queue, $ch;
  }
}

open my $fh, '>', $output or die "Cannot write $output: $!\n";
print $fh "efo_id\tdescription\n";

my @ids = sort { _efo_num($a) <=> _efo_num($b) || $a cmp $b } keys %visited;
for my $id (@ids) {
  my $t = $ont->find_term(-identifier => $id);
  my $name = $t ? ($t->name // '') : '';
  next if $name eq '';
  print $fh "$id\t$name\n";
}
close $fh;

print "Found " . scalar(keys %visited) . " cancer-related terms (including root).\n";
print "Wrote output to $output\n";

sub _efo_num {
  my ($id) = @_;
  return 9_999_999_999 unless defined $id;
  return $1 if $id =~ /^EFO:(\d{7})$/;
  return 9_999_999_999;
}