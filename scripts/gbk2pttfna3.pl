#!/usr/bin/env perl -w
use warnings;
use strict;
use File::Basename;

# This script takes a GenBank file as input, and produces a
# Fasta file and a NCBI PTT file (protein table) as output.
# A PTT file is a line based, tab separated format with fixed column types.
# 
# Written by Torsten Seemann
# 18 September 2006
# Modified by Victor de Jager 17-November 2009 to acomodate exon information instead of full CDS length
# All exons are written as single entries, but with the same GI and description.
# These exons should be joined to one protein, either in forward or using the reverse complement
# Re-written by Egon Ozer 25 January 2023 to remove BioPerl dependency and filter out redundant records.

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
my @c_seqs;
my %crecs;

#############################

sub gbk_convert{
    my $file = shift;
    my $return_status = 0;
    my $shortfile = basename($file);
    return(1) unless -e $file;
    my $gbkin;
    if ($file =~ m/\.gz$/){
        open ($gbkin, "gzip -cd $file | ") or return(1);
    } else {
        return(5) if -B $file;
        open ($gbkin, "<", $file) or return(1);
    }
    my $loccount = 0;
    my $seqcount = 0;
    my ($c_id, $c_seq);
    my $is_prod;
    my @tags = ("-") x 6;
    my @ptags;
    my @ctg_order;
    my $reading = 1; # 1 = front material, 2 = annotations, 3 = sequence

    while (my $fline = <$gbkin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            next if $line =~ m/^\s*$/;
            if ($line =~ m/^LOCUS\s+\S*\s+\d+\sbp/){
                if ($reading == 2){ #no ORIGIN sequence record was found between LOCUS records
                    return (2);
                }
                if ($reading == 3){ #this shouldn't happen. Should see a `//` between records 
                    if ($c_seq and $c_id){
                        push @c_seqs, ([$c_id, $c_seq]);
                        $c_seq = "";
                        $reading = 1;
                        print "More than one sequence record found in gbk file. Only using first sequence.\n";
                        return(0); ### Added this to only use the first sequence record in the file, similar to prior essentials/gbk2pttfna behavior.
                    } else {
                        return (2);
                    }
                }
            }
            if ($line =~ m/^\/\//){ #reached the end of the file (or record)
                if ($c_seq and $c_id){
                    push @c_seqs, ([$c_id, $c_seq]);
                    $c_seq = "";
                    $reading = 1;
                    return(0); ### Added this to only use the first sequence record in the file, similar to prior essentials/gbk2pttfna behavior.
                } else {
                    return (2);
                }
            }
            if ($reading == 1){
                if ($line =~ m/^LOCUS\s+([^\s]+)/){
                    $seqcount++;
                    if ($line =~ m/^LOCUS\s+(\S+)\s+\d+ bp/){
                        $c_id = $1;
                    } else {
                        $c_id = "rec$seqcount";
                    }
                    push @ctg_order, $c_id;
                    next;
                }
                if ($line =~ m/^FEATURES\s+Location\/Qualifiers/){
                    $reading = 2;
                    next;
                }
            } elsif ($reading == 2){
                if ($line =~ m/^\s+(\S+)\s+(complement\()*[<>]*(\d+)<*\.\.[<>]*(\d+)>*\)*\s*$/){
                    $is_prod = "";
                    my ($type, $start, $stop) = ($1, $3, $4);
                    my $dir = "+";
                    $dir = "-" if $2;
                    unless ($type eq "source" or $type eq "gene" or $crecs{$c_id}{$start}{$stop}{$dir}){
                        @{$crecs{$c_id}{$start}{$stop}{$dir}} = ($type);
                    }
                    if ($type eq "CDS"){
                        $loccount++;
                    }
                    if ($ptags[0]){
                        while (@ptags){
                            my ($o_start, $o_stop, $o_dir) = (shift @ptags, shift @ptags, shift @ptags);
                            #if ($tags[5] eq "p"){ ## skip pseudogenes (or not, I guess)
                            #    delete($crecs{$c_id}{$o_start}{$o_stop}{$o_dir});
                            #} else {
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[0];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[1];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[3] = $tags[2];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[4] = $tags[3];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[5] = $tags[4];
                                $loccount++;
                            #}
                        }
                        @tags = ("-") x 6
                    }
                    @ptags = ($start, $stop, $dir);
                    if ($type eq "source" or $type eq "gene"){
                        undef @ptags;
                        @tags = ("-") x 6
                    }
                    next;
                }

                ## Handle split genes. Could either ignore, provide with suffixes, or ignore "introns" and just take
                ## the first and last coordinates in the set as start and stop (this is the old behavior of gbk2pttfna3)  
                ## Including (with suffixes) will be afffected by truncation, i.e. each "exon" will be 3' truncated 
                ## rather than just the end of the gene sequence.
                ## Will use first and last coordinates as start & stop.
                if ($line =~ m/^\s+(\S+)\s+(complement)*\S*join\(([<>]*\d+<*\.\.[<>]*\d+>*,[^\)]+)\)\S*\s*$/){
                    $is_prod = "";
                    my ($type, $coords) = ($1, $3);
                    my $dir = "+";
                    $dir = "-" if $2;
                    $coords =~ m/^(\d+).*\.(\d+)$/;
                    my ($start, $stop) = ($1,$2);
                    unless ($type eq "source" or $type eq "gene" or $crecs{$c_id}{$start}{$stop}{$dir}){
                        @{$crecs{$c_id}{$start}{$stop}{$dir}} = ($type);
                    }
                    if ($type eq "CDS"){
                        $loccount++;
                    }
                    if ($ptags[0]){
                        while (@ptags){
                            my ($o_start, $o_stop, $o_dir) = (shift @ptags, shift @ptags, shift @ptags);
                            #if ($tags[5] eq "p"){ ## skip pseudogenes (or not, I guess)
                            #    delete($crecs{$c_id}{$o_start}{$o_stop}{$o_dir});
                            #} else {
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[0];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[1];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[3] = $tags[2];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[4] = $tags[3];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[5] = $tags[4];
                                $loccount++;
                            #}
                        }
                        @tags = ("-") x 6
                    }
                    @ptags = ($start, $stop, $dir);
                    if ($type eq "source" or $type eq "gene"){
                        undef @ptags;
                        @tags = ("-") x 6
                    }
                    next;
                }


                if ($line =~ m/^ORIGIN\s*$/){
                    $is_prod = "";
                    if ($ptags[0]){
                        while (@ptags){
                            my ($o_start, $o_stop, $o_dir) = (shift @ptags, shift @ptags, shift @ptags);
                            #if ($tags[5] eq "p"){ ## skip pseudogenes (or not, I guess)
                            #    delete($crecs{$c_id}{$o_start}{$o_stop}{$o_dir});
                            #} else {
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[0];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[1];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[3] = $tags[2];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[4] = $tags[3];
                                ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[5] = $tags[4];
                                $loccount++;
                            #}
                        }
                    }
                    undef @ptags;
                    @tags = ("-") x 6;
                    $reading = 3;
                    next
                }
                ## tags: [0]=locus_tag, [1]=product, [2]=gene, [3]=gi, [4]=cog, [5]=pseudo
                if ($line =~ m/^\s+\/(\S+)=\"*([^"]*)\"*/){
                    $is_prod = "";
                    my ($key, $val) = ($1, $2);
                    if ($key eq "locus_tag"){
                        $tags[0] = $val;
                    }
                    if ($key eq "product"){
                        $tags[1] = $val;
                        $is_prod = 1;
                        $tags[4] = $1 if $val =~ m/^(COG\S+)/;
                    }
                    if ($key eq "gene"){
                        $tags[2] = $val;
                    }
                    if ($key eq "db_xref"){
                        my $gi = $1 if $val =~ m/\bGI:(\d+)\b/;
                        if ($gi){
                            if ($tags[3] eq "-"){
                                $tags[3] = $gi;
                            } else {
                                $tags[3] .= ",$gi";
                            }
                        }
                    }
                    next;
                }
                if ($line =~ m/^\s+\/pseudo\s*$/){
                    $tags[5] = "p";
                }
                if ($is_prod){
                    $line =~ s/^\s*//;
                    $line =~ s/"*\s*$//;
                    $tags[1] .= " $line";
                }
            } elsif ($reading == 3){
                $line =~ s/\d//g;
                $line =~ s/\s//g;
                $c_seq .= $line;
                next;
            }
        }
    }
    if ($c_seq and $c_id){
        push @c_seqs, ([$c_id, $c_seq]);
        $c_seq = "";
        $reading = 1;
    }
    close ($gbkin);
    return(0);
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

## read in genbank file
my $status = gbk_convert($ini{'-gbk'}{val});
die "ERROR: $0 died with status $status\n" if $status;

## output PTT
open (my $ptt,">$ini{'-sa'}{val}") or die "Cannot create PTT file";
print $ptt join("\t", qw(!Locus Start Stop Strand Length PID Gene Product))."\n";
my %seen_lid;
foreach my $contig (sort keys %crecs){
    foreach my $start (sort {$a <=> $b} keys %{$crecs{$contig}}){
        foreach my $stop (sort {$a <=> $b} keys %{$crecs{$contig}{$start}}){
            foreach my $dir (sort keys %{$crecs{$contig}{$start}{$stop}}){
                my ($type, $lid, $prod, $gene, $gi, $cog) = @{$crecs{$contig}{$start}{$stop}{$dir}};
                my $length = ($stop - $start) + 1;
                next if ($length < 0); ## the only time this should happen is if a gene spans the end / beginning of the sequence. Safest to just ignore.
                unless($lid){
                    print STDERR "$start, $stop, $dir\n";
                }
                unless ($lid eq "-" or $seen_lid{$lid}){
                    print $ptt "$lid\t$start\t$stop\t$dir\t$length\t$gi\t$gene\t$prod\n";
                }
                $seen_lid{$lid} = 1;
            }
        }
    }
}
close ($ptt);

## output fasta
open (my $fa, ">$ini{'-s'}{val}") or die "Cannot create fasta file";
print $fa ">$c_seqs[0][0]\n";
my $seqleng = length($c_seqs[0][1]);
for (my $j = 0; $j < $seqleng; $j += 60){
    my $segment = uc(substr($c_seqs[0][1], $j, 60));
    print $fa "$segment\n";
}



