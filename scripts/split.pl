#!/usr/bin/perl
# Read a fasta file and extract the sequence data

use strict;
use warnings;

# Declare and initialize variables
my @file_data = (  );
my $dna = '';
# Read in the contents of the file "sample.dna"
@file_data = get_file_data($ARGV[0]);
# Extract the sequence data from the contents of the file "sample.dna"
$dna = extract_sequence_from_fasta_data(@file_data);
# Print the sequence in lines 20 characters long
print_sequence($dna, $ARGV[1]);
exit;

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

# extract_sequence_from_fasta_data
#
# A subroutine to extract FASTA sequence data from an array
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
        } elsif($line =~ /^\s*#/) {
            next;

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

# print_sequence
#
# A subroutine to format and print sequence data 
sub print_sequence {
    my($sequence, $length) = @_;
    my $posgenome;
    use strict;
    use warnings;
    # Print sequence in lines of $length
    for ( my $pos = 0 ; $pos < length($sequence) ; $pos += 1 ) {
	$posgenome = $pos+1;
	print ">$posgenome\n";
        print substr($sequence, $pos, $length), "\n";
    }
}
