#!/usr/bin/perl

use strict;
use warnings;

my $usage = "
tavsta.pl <genomes.ta>

Replaces the step in the doTAcounts subroutine that relies on pass since
newer versions of pass no longer seem to work and this does the job much
quicker, as long as only unique sites are requested. The output is sorted
lexographically since that's the way Essentials does it.

";

die $usage unless @ARGV;

open (my $in, "<$ARGV[0]") or die "ERROR: Can't open $ARGV[0]: $!\n";
my %hash;
my $id;
#my $incount = 0;
while (my $line = <$in>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    if ($line =~ m/^>/){
        $id = substr($line, 1);
        next;
    }
    $line =~ s/\s//g;
    $line = uc($line);
    push @{$hash{$line}}, $id;
    #$incount++;
}
close ($in);
#print STDERR "incount: $incount\n";

my @keep;
foreach my $key (keys %hash){
    my @array = @{$hash{$key}};
    next if scalar @array > 1;
    push @keep, $array[0];
}
#print STDERR "outcount: ", scalar @keep, "\n";

@keep = sort{$a cmp $b}@keep;
foreach (@keep){
    print "1\t$_\n";
}
