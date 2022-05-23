#!/usr/bin/env perl

# junction_mapping.pl
# Recover all original read headers after dereplicating
# and clustering barcodes in junction mapping.
# Also annotates cluster IDs

use v5.26;
use warnings;
use autodie;
use Text::CSV qw( csv );
#use DDP;

die "Usage: $0 <clustered.uc> <derep.uc>"
  unless @ARGV == 2;

my $cluster_file = shift;
my $derep_file   = shift;
my %centroids;

open my $cf, "<:encoding(utf8)", $cluster_file;
my $csv = Text::CSV->new(
    {
        binary    => 1,
        auto_diag => 1,
        sep       => "\t",
    }
);
$csv->column_names(
    qw(type cid slen_csize pident
      strand u1 u2 cigar qseqlab cseqlab)
);

while ( my $row = $csv->getline_hr($cf) ) {
    my ( $header, $derep_size ) = split( ";", $row->{'qseqlab'} );
    $derep_size = ( split /=/, $derep_size )[1];

    # If cluster record
    if ( $row->{'type'} eq 'C' ) {
        $centroids{ $row->{'cid'} }->{'csize'} = $row->{'slen_csize'};

        die "Error: centroid sequence header does not match"
          unless $header eq $centroids{ $row->{'cid'} }->{'cseqlab'};
    }

    # If centroid record, just save header and cid
    elsif ( $row->{'type'} eq 'S' ) {

        #p $row->{'qseqlab'}, as => "Cluster $row->{'cid'} centroid header";
        $centroids{ $row->{'cid'} }->{'calcsize'} += $derep_size;

        # Annotate cluster with centroid seq
        $centroids{ $row->{'cid'} }->{'cseqlab'} = $header;

        push $centroids{ $row->{'cid'} }->{'hlist'}->@*, $header;
    }

    # If hit record, add to list of headers
    elsif ( $row->{'type'} eq 'H' ) {

#p $derep_size, as => "Cluster $row->{'cid'} derep size of hit $row->{'qseqlab'}";

        $centroids{ $row->{'cid'} }->{'calcsize'} += $derep_size;
        push $centroids{ $row->{'cid'} }->{'hlist'}->@*, $header;
    }

}

close $cf;

#say np %centroids;

# Now we read the derep cluster file
my %derep_hits;
open my $df, "<:encoding(utf8)", $derep_file;
while ( my $row = $csv->getline_hr($df) ) {
    if ( $row->{'type'} eq 'H' ) {
        push $derep_hits{ $row->{'cseqlab'} }->@*, $row->{'qseqlab'};
    }
}

#p %derep_hits;

# For clusters with at least size 10, retrieve read headers for all raw reads
# which might have been collapsed into cluster members.
for my $cid (
    sort { $a <=> $b }
    grep { $centroids{$_}->{'csize'} >= 10 } keys %centroids
  )
{
    # Go through each member of the cluster
    for my $derep ( $centroids{$cid}->{'hlist'}->@* ) {

        # Print cluster member header and its cluster ID
        say "$derep\t$cid";

        # If this member is one of a derep set, do the same for each derep read
        if ( exists $derep_hits{$derep} ) {
            say "$_\t$cid" for $derep_hits{$derep}->@*;
        }
    }
}
