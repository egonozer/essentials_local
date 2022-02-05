#!/usr/bin/env perl -w
use warnings;
use strict;

#use fgweb::general::general ;
#use fgweb::general::numerics ;
#use fgweb::general::cmd_input ;

use Bio::SeqIO;

# This script takes a GenBank file as input, and produces a
# Fasta file and a NCBI PTT file (protein table) as output.
# A PTT file is a line based, tab separated format with fixed column types.
# 
# Written by Torsten Seemann
# 18 September 2006
# Modified by Victor de Jager 17-November 2009 to acomodate exon information instead of full CDS length
# All exons are written as single entries, but with the same GI and description.
# These exons should be joined to one protein, either in forward or using the reverse complement

########################

my %ini ;

$ini{'-gbk'}{desc} = 'Subject GenBank file; the file containing the reference DNA sequence(s) and annotation(s).' ;
$ini{'-gbk'}{val}  = '' ;
$ini{'-gbk'}{checkisexisting} = 1 ;
$ini{'-gbk'}{checkisstring} = 1 ;
$ini{'-gbk'}{mandatory} = 1 ;

$ini{'-s'}{desc} = 'Subject fasta output file.' ;
$ini{'-s'}{val}  = '' ;
#$ini{'-s'}{checkisexisting} = 1 ;
$ini{'-s'}{checkisstring} = 1 ;
$ini{'-s'}{mandatory} = 1 ;

$ini{'-sa'}{desc} = 'Subject ptt file; the output file containing the annotation for the reference DNA sequence(s).' ;
$ini{'-sa'}{val}  = '' ;
#$ini{'-sa'}{checkisexisting} = 1 ;
$ini{'-sa'}{checkisstring} = 1 ;
$ini{'-sa'}{mandatory} = 1 ;
 
#global variables

my @PTTFILE; # will hold all lines for the PTT file.
my $gbk;
my @cds;
my @gene;
my @trna;
my @rrna;
my $seq;

#############################

#sub parseparams {       # shows the usage information if applicable
#    my $title = "\nBy Victor de Jager [vcitor.de.jager\@nbic.nl]\n\n" ;
#    
#    $title   .= "This script takes a GenBank file as input, and produces a Fasta file and PTT annotation file .\n" ;
#    $title   .= "The following steps are performed:\n" ;
#    $title   .= "1 Parses the genbank file.\n" ;
#    $title   .= "2. Creates a Fasta file\n" ;
#    $title   .= "3. Creates a PTT file\n" ;
#    cmd_input::check_cmd_input(\%ini, \@ARGV,$title) ;
#    
#    general::error('GenBank file does not exist') if ($ini{'-gbk'}{val} ne '' and !-e $ini{'-gbk'}{val}) ;
#    general::error('An output ptt file is mandatory') if ($ini{'-sa'}{val} eq '');
#    general::error('An output fasta file is mandatory') if ($ini{'-s'}{val} eq '');
#
#} # parseparams


sub tag {
   my($f, $tag) = @_;
   return '-' unless $f->has_tag($tag);
   return join(' ', $f->get_tag_values($tag));
}

sub create_header {
	my $seq = shift;
	my $cds = shift;
	my $pttfile = shift;
	
#	push (@PTTFILE, $seq->description." - 0..".$seq->length."\n") ;
#	push (@PTTFILE, scalar(@{$cds})." proteins\n");
	push (@PTTFILE, join("\t", qw(!Locus Start Stop Strand Length PID Gene Product))."\n") ;
	
}
sub parse_cds{
	my $cds=shift;

	for my $f (@{$cds}) {
    	my $gi = '-';
    	$gi = $1 if tag($f, 'db_xref') =~ m/\bGI:(\d+)\b/;
    	my $cog = '-';
    	$cog = $1 if tag($f, 'product') =~ m/^(COG\S+)/;
    	my $start;
    	my $end;
    	my $strand;
    	my $pttline;
   
    	#Test whether the cds is a split location
   
    	if ( $f->location->isa('Bio::Location::SplitLocationI')){
			#ignoring join() features.	
			for my $loc ( $f->location->sub_Location ) {
    			$start=$loc->start;
    			$end= $loc->end;
  		        #&create_columns($f,$start,$end,$gi,$cog);
  			}
		
   	 	}
   	 	else {
   			#We are dealing with a simple feature
   			$start=$f->start;
    			$end= $f->end;
   			&create_columns($f,$start,$end,$gi,$cog);
   		}		
   	 
	}
}

sub create_columns {
	
	my $f=shift;
	my $start=shift;
	my $end = shift;
	my $gi = shift;
	my $cog = shift;
	;
	
	my @col = (
     
     tag($f, 'locus_tag'),
     $start,
     $end,
     $f->strand >= 0 ? '+' : '-',
     (int(($f->length))),
     $gi,
     tag($f, 'gene'),
     tag($f, 'product'),
   );
   my $col = join("\t", @col);
   push (@PTTFILE,$col."\n"); 
}

sub write_ptt{
	
	my $file = shift;
	#open PTT,">$file" or general::error('Cannot create PTT file');
    open PTT,">$file" or die "Cannot create PTT file";
	
	for my $ptt (@PTTFILE){
		print PTT $ptt;
	}
	close(PTT);
}

sub write_fna{
	
	my $file = shift;
	my $seq = shift;
	my $fna = Bio::SeqIO->new( -file=> ">$file",
							   -format=>'fasta');
	$fna->write_seq($seq);
	
}

#################################
######
######    MAIN
######
#################################


print "Converting genbank to ptt file and fna file\n";
#&parseparams();
die "
ERROR: Enter the following three arguments:
1) path to genbank file
2) output genome fasta file name
3) output genome ptt file name
" unless @ARGV >= 3;
$ini{'-gbk'}{val} = $ARGV[0];
$ini{'-s'}{val} = $ARGV[1];
$ini{'-sa'}{val} = $ARGV[2];

$gbk = Bio::SeqIO->new(	-file=> $ini{'-gbk'}{val},
							-format=>'genbank');
$seq = $gbk->next_seq;
@cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;
@trna = grep { $_->primary_tag eq 'tRNA' } $seq->get_SeqFeatures;
@rrna = grep { $_->primary_tag eq 'rRNA' } $seq->get_SeqFeatures;
@gene = grep { $_->primary_tag eq 'gene' } $seq->get_SeqFeatures;

push (@cds, @trna);
push (@cds, @rrna);
push (@cds, @gene);

&create_header($seq,\@cds);
&parse_cds(\@cds);
&write_ptt($ini{'-sa'}{val});
&write_fna($ini{'-s'}{val}, $seq);


