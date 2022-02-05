#!/usr/bin/env perl

use strict;
use warnings;

# Declare and initialize variables
my @file_data = (  );
my $query = $ARGV[1];
my $tasitelength = $ARGV[2];
my $dna = '';
my $regexp = '';
my @locations = (  );

# Read in the file "sample.dna"

@file_data = get_file_data($ARGV[0]);

# Extract the DNA sequence data from the contents of the file
$dna = extract_sequence_from_fasta_data(@file_data);

# check for empty query
do {
    chomp $query;
    # Exit if empty query
    if ($query =~ /^\s*$/ ) {
        exit;
    }

# Perform the search in the DNA sequence
@locations = match_positions($query, $dna);

# Report the sequences and the positions to the user as fasta
if (@locations) {
    foreach (@locations) {
	print "$_\n";
    }
} else {
    print "A $query site is not in the DNA:\n";
}

exit;

#########################################################################
#
# Subroutine
#
# get_file_data
#
# A subroutine to get data from a file given its filename

sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    @filedata = <GET_FILE_DATA>;
    close GET_FILE_DATA;
    return @filedata;
}

###########################################################################
# Subroutine
#
# extract_sequence_from_fasta_data
#
# A subroutine to extract sequence data from a fasta file in an array

sub extract_sequence_from_fasta_data {

    my(@fasta_file_data) = @_;

    use strict;
    use warnings;

    # Declare and initialize variables
    my $sequence = '';

    foreach my $line (@fasta_file_data) {
        # discard blank line
        if ($line =~ /^\s*$/) {
            next;
        # discard comment line
        #} elsif($line =~ /^\s*#/) {
        #    next;
        # discard fasta header line
        } elsif($line =~ /^>/) {
            next;
        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }
    # remove non-sequence data (in this case, whitespace) from $sequence string
    $sequence =~ s/\s//g;
    return $sequence;
}

################################################################################
#
# Subroutine
#
# Find locations of a match of a regular expression in a string on the forward strand
#
# 
# return an array of sequences that contain the regular expression
#

sub match_positions {
    my($regexp, $sequence) = @_;
    use strict;

    # Declare variables
    #
    my @positions = (  );
    #
    # Determine positions and sequences surrounding regular expression matches
    #
    while ( $sequence =~ /$regexp/ig ) {
	#pos or substr are slow. according to strace, the string gets reread from start everytime.
	#push ( @positions, '>'.(pos($sequence) - length($&) + 1));
	#push ( @positions, substr($sequence, pos($sequence) -length($&), $tasitelength));
	#push ( @positions, '>'.(0-(pos($sequence) - length($&) + 1)));
	#push ( @positions, reverse_complement(substr($sequence, pos($sequence) -$tasitelength, $tasitelength)));
	
	#$1 does not work (??)
	push (@positions, '>'.($-[0]+1)) ;
	push (@positions, substr($sequence, $-[0], $tasitelength));
	push (@positions, '>-'.($-[0]+1)) ;
	push (@positions, reverse_complement(substr($sequence, $-[0]-$tasitelength+2, $tasitelength)));
    }
    return @positions;
}


################################################################################
# Subroutine
#
# subroutine for working out the reverse complement of a sequence 
# 
# Returns reverse complement

sub reverse_complement { 
    my $seq = shift @_; 
    my $rev = reverse $seq; # reverses string 
    $rev =~ tr/gatcGATC/ctagCTAG/;  # complement nts 
    return $rev; 
} 


}
