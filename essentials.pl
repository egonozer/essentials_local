#!/usr/bin/env perl

my $version = "2.0";

## Simplified local version of Essentials that uses external barcode trimmer and sorter in TIS tools

## Calls scripts for parsing the genbank (module 1) and generating the insertion site flanking sequences (module 2), 
# parses configfile, downloads reads, does alignments  (module 3)  and runs statistical analysis module 4 (R)
# prerequisites: R, GNU tools, curl, pass. fq_all2std_noqual.pl, fastx_barcode_splitter_trimmer, gbk2pttfna.pl, tasite.pl
# fastx_barcode_splitter_trimmer requires Text::LevenshteinXS if you want to split your barcodes quickly

use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw ( catfile path );
use Cwd qw(cwd abs_path);


my $scriptpath = dirname(__FILE__);

## Set defaults
my $insertion           = "TA";
my $normalization       = "TMM";
my $librarysize         = "15000"; ## Remember this needs to be multipled by 2 later
my $statmethod          = "qCML";
my $loess               = "Yes";
my $adjustment          = "BH";
my $dispersion          = "tagwise";
my $smoothing           = 5;
my $ptt                 = "truncated.ptt";
my $zip                 = "No";
my $tasitelength        = 24; #mode read length. 17 for Goodman, 24 for Boll 
my $repeat              = "Yes"; ## This is not used anymore. Was used to determine the "-uniq" setting for pass aligner. I think unsetting it was a way to bypass removing non-unique TA site flanks, though why?

my $usage = "
essentials.pl

Required:
  -c or --config    Path to configuration file. Should contain path to read count
                    files in wiggle format (as output by INSeq_read_preprocess.pl)
                    group, and ID, separated by tabs. Group can be anything you
                    want to separate two groups to be compared, usually 'control'
                    and 'target'. Use 'ignore' as a group to skip lines.
                    Example:
                    /path/to/Pool1.wiggle   control Control_1
                    /path/to/Pool2.wiggle   target  Treatment_2
                    etc.
                    Will also autodetect old 8-column Essentials config files.
  -g or --gbk       Path to genbank-formatted sequence and annotation file. If
                    you are using a settings file (--settings) that contains the
                    genbank file path this can be omitted and will be ignored

Optional:
  -l or --libsize   Expected library size (default: $librarysize)
  --rdleng          Expected length of reads used for alignment. Options:
                    24 = Boll protocol (default)
                    17 = Goodman protocol
  --insert          Transposon insertion site type. Options are:
                    'TA' (default) or 'random' (NONFUNCTIONAL)
  --full            Use full gene lengths (default: use 5' truncated genes)
  --norepeat        Do not count repeats (default: repeats are counted)
  --noloess         Skip loess normalization of intra-pool read counts
                    (default: loess normalization is performed)
  --norm            Normalization method across pools. Options are:
                    'TMM' = trimmed mean of M-values (default)
                    'TMMwsp' = TMM with singleton pairing
                    'RLE' = relative log expression
                    'upperquartile' = scale factor from 75% quartile of counts
                    'none' = no normalization
  --stat            Statistical method used. Options are:
                    'qCML' = quantile-adjusted conditional maximum likelihood
                             (default)
                    'CR' = Cox-Reid profile-adjusted likelihood
  --disp            Dispersion estimates. Options are:
                    'tagwise' (default) or 'common'
  --smooth          Smoothing factor. (default: $smoothing)
  --adjust          p-value adjustment. Options are:
                    'BH' (default), 'holm', 'hochberg', 'hommel', 'bonferroni',
                    'BY', or 'none'
  --zip             Zip results files. Requires Archive::Zip perl module to be
                    installed. (default: results are not zipped)
  
  --settings        A settings file containing some or all of the above
                    settings as was used in previous versions of this program.
                    File settings will override any command line settings.

";

## input files
my $settingsfile;
my $genbank;
my $configfile;


use Getopt::Long;
Getopt::Long::Configure qw(gnu_getopt);

my ($noloess, $usezip, $usefull);
GetOptions(
    'settings=s'    => \$settingsfile,
    'gbk|g=s'       => \$genbank,
    'config|c=s'    => \$configfile,
    'rdleng=i'      => \$tasitelength,
    'insert=s'      => \$insertion,
    'norm=s'        => \$normalization,
    'libsize|l=i'   => \$librarysize,
    'stat=s'        => \$statmethod,
    'noloess'       => \$noloess, ## default is loess normalization
    'adjust=s'      => \$adjustment,
    'disp=s'        => \$dispersion,
    'smooth=i'      => \$smoothing,
    'nozip'         => \$usezip, ## default is no zip
    'full'          => \$usefull ## default is truncated genes
) or die $usage;

die $usage unless $configfile;
die $usage unless ($settingsfile or $genbank); ## Die later if there's no genbank file in the settings file.

## No longer used settings  
my ($barcodemismatch, $barcode_end, $transposon_end, $minhit, $strand, $minomicsreads);

if ($settingsfile){
    open (my $sin, "<settingsfile") or die "ERROR: Can't open $settingsfile: $!\n";
    while (my $line = <$sin>){
        chomp $line;
        if ($line =~ m/^([^=]+)=(.*)/){
            my ($type, $val) = ($1, $2);
            $genbank = $val if $type eq 'TEMPLATE';
            $insertion = $val if $type eq 'INSERTIONSITE';
            $normalization = $val if $type eq 'NORMALIZATION';
            $barcodemismatch = $val if $type eq 'BARCODEMISMATCH';
            $barcode_end = $val if $type eq 'BARCODE_END';
            $transposon_end = $val if $type eq 'TRANSPOSON_END';
            $minhit = $val if $type eq 'MINHIT';
            $tasitelength = $val if $type eq 'TASITELENGTH';
            $strand = $val if $type eq 'STRAND';
            $librarysize = 2 * $val if $type eq 'LIBRARYSIZE';
            $statmethod = $val if $type eq 'STATMETHOD';
            $adjustment = $val if $type eq 'ADJUSTMENT';
            $dispersion = $val if $type eq 'DISPERSION';
            $smoothing = $val if $type eq 'SMOOTHING';
            $ptt = $val if $type eq 'GENETRUNCATION';
            $loess = $val if $type eq 'LOESS';
            $repeat = $val if $type eq 'REPEAT';
            $minomicsreads = $val if $type eq 'MINOMICSREADS';
            $zip = $val if $type eq 'DOZIP';
        }
    }
}
## Sanity checking
die $usage unless ($genbank);
$tasitelength =~ s/\s//g;
die $usage unless $tasitelength =~ m/^\d+$/; ## only allow integers
$librarysize =~ s/\s//g;
die $usage unless $librarysize =~ m/^\d+$/; ## only allow integers
$insertion =~ s/\s//g;
die $usage unless ($insertion eq "TA" or $insertion = "random");
my %normopts = (
    'TMM' => 1,
    'TMMwsp' => 1,
    'RLE' => 1,
    'upperquartile' => 1,
    'none' => 1
);
$normalization =~ s/\s//g;
die $usage unless $normopts{$normalization};
$statmethod =~ s/\s//g;
die $usage unless ($statmethod eq "qCML" or $statmethod eq "CR");
$dispersion =~ s/\s//g;
die $usage unless ($dispersion eq "tagwise" or $dispersion eq "common");
$smoothing =~ s/\s//g;
die $usage unless $smoothing =~ m/^\d+$/; ## only allow integers
my %adjopts = (
    'BH' => 1,
    'holm' => 1,
    'hochberg' => 1,
    'hommel' => 1,
    'bonferroni' => 1,
    'BY' => 1,
    'none'  => 1
);
$adjustment =~ s/\s//g;
die $usage unless $adjopts{$adjustment};
$loess = "No" if $noloess;
$zip = "Yes" if $usezip;
$ptt = "genome.ptt" if $usefull;

######## Generate a settings output file if none was given. For record keeping #########

$librarysize = $librarysize * 2;

my $workdir = cwd;
## Define required executables
my $Rnew;
my $path = is_path("R");
if ($path){
    $Rnew = $path;
} else {
    die "ERROR: R is required but could not be found in your PATH.\n";
}
my $pass;
#my $pass  = "/Users/hauserlab/Applications/pass_v2.12/bin/pass";
my $curl  = "curl";
my $cat   = "cat";
my $cut   = "cut";
my $mv    = "mv";
my $grep  = "grep";
my $sed   = "sed";
my $tr    = "tr";
my $sort  = "sort";
my $uniq  = "uniq";
my $paste = "paste";
my $perl  = "$^X";

#fix if using random (probably a quick perl script):
my $revseq;
#my $revseq = '/Users/hauserlab/Applications/EMBOSS-6.4.0/emboss/revseq';

## perl scripts
my $gbk2fna3        = "$scriptpath/scripts/gbk2pttfna3.pl";
my $tasitescript 	= "$scriptpath/scripts/tasite.pl"; 
my $split		    = "$scriptpath/scripts/split.pl";
my $tavsta          = "$scriptpath/scripts/tavsta.pl";

# R scripts for stats
my $loessscript		= "$scriptpath/scripts/loess.R";
my $loessgenescript	= "$scriptpath/scripts/loess_genes.R";
my $edger           = "$scriptpath/scripts/edgeR.R";

my $archivefile = 'archive.zip';

######################################################################################
# Run programs down here. Comment out sections to speed up a rerun or to debug.

#cleanup previous run
unlink("allinsertions.txt");
unlink("samplenames.txt");
unlink (<*.error>);

# Parsing genbank
print "Parsing selected genbank file. \n";
qx("$perl" "$gbk2fna3" "$genbank" genome.fasta genome_all.ptt 2>> gbk2fna.error);

# sorting ptt file with genes and proteins on genename and then be annotation length
my @pttdata;
my @ordered_pttdata;
open(ALLPTT, "genome_all.ptt") or general::error("Can't open file: genome_all.ptt");
open (ALLPTT2, "> genome_all_sorted.ptt");

while (<ALLPTT>){
    chomp;
    push(@pttdata, $_);
} 

@ordered_pttdata = sort { (split '\t', $a)[0]  cmp (split '\t', $b)[0] || (length $b  <=> length $a) } @pttdata;
print ALLPTT2 join( "\n", @ordered_pttdata )."\n";
close ALLPTT;
close ALLPTT2;


# TODO: get rid of crap below and do the duplicate search on the array.
my $alllocus;
my $allstart;
my $allstop ;
my $allstrand ;
my $alllength ;
my $allpid ;
my $allgene ;
my $allproduct ;
my $alllocusold ;
my $allpttline;
my @splitallptt;
$alllocusold = 0;

open(ALLPTT3, "genome_all_sorted.ptt") or die ("Can't open file: genome_all_sorted.ptt");
open (PTT, "> genome.ptt");

while (<ALLPTT3>){
    $allpttline = $_ ;
    chomp $allpttline ;
    @splitallptt = split '\t', $allpttline;
    $alllocus = $splitallptt[0] ;
    $allstart = $splitallptt[1] ;
    $allstop = $splitallptt[2] ;
    $allstrand = $splitallptt[3] ;
    $alllength =  $splitallptt[4] ;
    $allpid = $splitallptt[5];
    $allgene = $splitallptt[6];
    $allproduct = $splitallptt[7];
    if ($alllocus ne $alllocusold) {
        print PTT "$alllocus\t$allstart\t$allstop\t$allstrand\t$alllength\t$allpid\t$allgene\t$allproduct\n";
    }
    $alllocusold = $alllocus;
}
close PTT;
close ALLPTT3;

# Truncate end of genes to filter out function-retaining C-terminal transposon insertions
truncategenes();

# Finding TA insertion site flanking sequences. 
if ($insertion eq "TA") {
    print "Searching TA sites and indexing insertion site flanking sequences ($tasitelength bp) of selected genome\n";
    qx("$perl" "$tasitescript" genome.fasta "$insertion" "$tasitelength" >genome.ta 2>> tasitescript.error);
}

# Fragmenting genomic DNA for detection of non-unique insertion site flanking sequences.
if ($insertion eq "random") {
    print "Indexing insertion site flanking sequences ($tasitelength bp) of selected genome\n";
    qx("$perl" "$split" genome.fasta "$tasitelength" >genome.ta_fw 2>>splitscript.error);
    qx("$revseq" -sequence genome.ta_fw -outseq genome.ta_rev);
    qx("$cat" genome.ta_fw genome.ta_rev > genome.ta);
}

# Run the analysis
get_and_convert();

# Run the statistics
doStats();

# zip everything 
if ($zip eq "Yes") {
    
    use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
    
    my $archive = Archive::Zip->new();
    my $dirfile ;
    
    unlink ($archivefile);
    opendir DIR, $workdir or die "cannot open dir $workdir: $!";
    my @files= readdir DIR;
    closedir DIR;
    foreach $dirfile (@files) {if ($dirfile ne "."){if ($dirfile ne ".."){$archive->addFile($dirfile);}}}
    unless ( $archive->writeToFileNamed($archivefile) == AZ_OK ) {     
        die 'write error';
    }

}


######################################################################################
# get_and_convert
# Main subroutine. Everything else is started from here. 
# This subroutine will parse the configfile, download the read files in fastq or
# fasta format and will convert all files to fasta. Also it will generate the required 
# files for other subroutines

sub get_and_convert {
    my @samples;
    my $barcode;
    my $inputformat;
    my $inputlink;
    my $filecount;
    my $lineconfig;
    my $transposon;
    my $sampletype;
    my $samplecount;
    my $samples;
    my $inputcompression;
    my $samplename;
    my $decompress;
    my $library ;
    my %config ;
    my @splot ;
    my @splut ;
    #unused    my $transposonmismatch;

    #print "Preparing headers of output files\n";
    writeheaders();
    
    # parse the configfile line by line. each line contains the link, the barcode, the transposon sequence, the sampletype (control or target)  and the input format
    # If the barcode or transposon is not present, use N as sequence. 
    
    print "\nParsing configfile\n";
    open(CONFIGFILE, $configfile) or die "Can't open file:" .$configfile;
    $filecount = 0;
    $samplecount = 0;
    while (<CONFIGFILE>){
        $lineconfig = $_ ;
        chomp $lineconfig ;
        @splot = split '\t', $lineconfig;
        ($inputlink, $barcode, $transposon, $sampletype, $library, $inputformat, $inputcompression, $samplename) = ("") x 8;
        if (scalar @splot == 3){
            ($inputlink, $sampletype, $samplename) = @splot;
        } else {
            $inputlink = $splot[0] ;
            $barcode = $splot[1] ;
            $transposon = $splot[2] ;
            $sampletype = $splot[3] ;
            $library =  $splot[4] ;
            $inputformat = $splot[5];
            $inputcompression = $splot[6];
            $samplename =  $splot[7];
        }
    
        next if ($inputlink eq "Link");
        next if ($inputlink eq "link");
        next if ($inputlink eq "#");
        next if ($inputlink eq "#link");
        next if ($inputlink eq "#Link");
        next if ($inputlink eq "");
        next if ($inputlink eq " ");
        $samplecount = $samplecount + 1;

        $filecount = $filecount +1;
        
        #define readfile and sample
        my $readfile = "reads".$filecount.".fasta";
        my $sample = "split_".$sampletype."_sample".$samplecount."";
        
        ## Process wiggle file
        print "Input file $inputlink is a wiggle file.\n";
        my $alltafile = "allta_".$sample.".counts.txt";
        open (my $in, "<$inputlink") or die "ERROR: Can't open $inputlink: $!\n";
        open (my $out, ">$alltafile") or die "ERROR: Can't open $alltafile for writing: $!\n";
        while (my $line = <$in>){
            chomp $line;
            my ($one, $two) = split("\t", $line);
            print $out "$two\t$one\n";
        }
        close ($in);
        close ($out);
        
        #create samplename file for later use
        if ($sampletype ne "ignore") {
            open (SAMPLENAME, ">> samplenames.txt");
            if ($samplename ne "") {print SAMPLENAME "".$samplename."_".$samplecount."_".$sampletype."\n";}
            if ($samplename eq "") {print SAMPLENAME "sample".$samplecount."_".$sampletype."\n";}
            close (SAMPLENAME);
        }
        
        #add sample to file with control and target samples for gene edgeR
        if ($sampletype ne "ignore") {
            my $countsgene = "gene_".$sample.".counts.txt";
            open (TARGETSGENES, ">> gene_targets.txt");
            print TARGETSGENES "".$countsgene."\t".$sampletype."\t".$sampletype."\t".$library."\n";
            close (TARGETSGENES);
            
            #add sample to file with all samples as control samples for essential genes edgeR
            open (TARGETSESSENTIALGENES, ">> essentialgenes_targets.txt");
            print TARGETSESSENTIALGENES "".$countsgene."\ttarget\ttarget\t".$library."\n";
            close (TARGETSESSENTIALGENES);
            
            #add sample to file with control and target samples for insertion site edgeR
            my $tacounts = "ta_".$sample.".counts.txt";
            open (TARGETSTA, ">> ta_targets.txt");
            print TARGETSTA "".$tacounts."\t".$sampletype."\t".$sampletype."\t".$library."\n";
            close (TARGETSTA);
            push(@samples, $sample);
        }
    }
    
    # align the fasta files of the samples to the insertion site flanking sequences, count the number of reads per insertion site and count the number of reads per gene
    pass_alignments(@samples);

    #calculate number of unique insertion site flanking sequences per gene
    doTAcounts();
    
    #annotate the insertion site flanking sequences
    doTAannotation();
    print "\nCounting total number of insertion site flanking sequences hit\n";
    countlines2("allinsertionsuniq.txt");
        
}#END get_and_convert

######################################################################################
# writeheaders
# this subroutine prepares some initial files.

sub writeheaders {
    #create initial gene_targets.txt readDGE object input file
    open (TARGETSGENE, '>gene_targets.txt');
    print TARGETSGENE "files \t group \t description \t library \n";
    close (TARGETSGENE);
    #create initial essentialgenes_targets.txt readDGE object input file
    open (TARGETSESSENTIALGENES, '>essentialgenes_targets.txt');
    print TARGETSESSENTIALGENES "files \t group \t description \t library \n";
    close (TARGETSESSENTIALGENES);
    #create initial ta_targets.txt readDGE object input file
    open (TARGETSTA, '>ta_targets.txt');
    print TARGETSTA "files \t group \t description \t library \n";
    close (TARGETSTA);
}# END writeheaders

#######################################################################################
# pass_alignments
# This subroutine executes the pass commands to align the reads to the insertion site flanking sequences or to the genome.  
# input is $sample.fasta, output of the aligner is named $tacounts

sub pass_alignments {
    my @samples=@_;
    my $sample;
    foreach $sample (@samples) {
        
        my $alltacounts = "allta_".$sample.".counts.txt";
        my $tacounts = "ta_".$sample.".counts.txt";
        my $wiggle = "ta_".$sample.".wiggle";
        
        #create header for ta counts file
        open (TACOUNTS, "> $tacounts");
        print TACOUNTS "counts\tnames\n";
        close (TACOUNTS);
        
        ## Removed alignments of fasta files using pass
        
        stat ($alltacounts);
        unless (-e _){die "sample ".$sample.".fasta was not created and wiggle file ".$alltacounts." was not found either.\n"}

        countwiggle($alltacounts);
        print "Retaining only the $librarysize insertion site flanking sequences with the highest numbers of reads\n";
        
        ## Sort read counts at TA sites from highest to lowest and keep first '$librarysize' sites
        qx("$cat" "$alltacounts" |"$sort" -n -r |head -n "$librarysize" > "$tacounts");
        
        ## Loess normalize the per-site read counts. Overwrites the input file of top raw read counts (ta_XXX.counts.txt)
    	if ($loess eq "Yes") {
            print "Loess normalizing ".$tacounts." for read count bias\n";
            qx("$Rnew" --vanilla --args "$tacounts" < "$loessscript" 2>>R.error);
        }    	         
    	                                                         
    	## Add the normalized read counts to a file containing normalized read counts from all the read pools in the analysis.  
        qx("$cat" "$tacounts" >>allinsertions.txt);
        ## This removes headers, removes - signs on positions, and sorts by site. Leave this here in case loess normalization is not performed. Otherwise unnecessary,  
        qx("$cat" "$tacounts" |"$grep" -v names|"$tr" -d '-'|"$sort" -n -k2 >$wiggle);
        
        qx("$cat" $ptt |"$sort" -k2 -n > genome.ptt.tmp);
        qx("mv" genome.ptt.tmp $ptt);
        print "Counting read frequencies of genes and insertion site flanking sequences of ".$sample."\n";
        # rapidly add all insertion site counts (the wiggle track of the genome) per gene, which positions are annotated in genome.ptt generated by gbk2fnaptt.
        my $featureend;
        my $featurestart;
        my $featureline;
        my $wiggleline;
        my $wigglepos;
        my $featurecount = 0;
        my $featurename;
        my $wiggledepth;
        my @splet;
        my @splyt;
        ## This file contains the total number of (normalized) reads per gene
        my $countsgene = "gene_".$sample.".counts.txt";
        my $lengthline ;
        my @buffer = (0) x 10000;
        my $filepos = 0 ;

        open (COUNTSGENE, "> $countsgene" ) or die("Can't open file: ".$countsgene."");
        print COUNTSGENE "counts\tnames\n";
        open(LOCATIONS, $ptt) or die("Can't open file: ".$ptt."");
        open(WIGGLE, "< $wiggle" ) or die("Can't open file: ".$wiggle."");
        while (<LOCATIONS>){
    	    $featureline = $_ ;
    	    chomp $featureline ;
    	    @splyt = split '\t', $featureline;
    	    $featurename =  $splyt[0];
    	    $featurestart = $splyt[1];
    	    $featureend = $splyt[2];
    	    next if ($featurename eq "!Locus");
    	    $featurecount = 0;
    	    $wiggledepth = 0;
    	    while (<WIGGLE>){
                $filepos = tell(WIGGLE);
                shift @buffer;
                push  @buffer, $filepos;
        	
                $wiggleline = $_ ;
                #$lengthline = $_;
                chomp $wiggleline ;
                @splet = split '\t', $wiggleline;
                $wigglepos = $splet[1];
                $wiggledepth = $splet[0];
                if ($featurestart <= $wigglepos && $featureend >= $wigglepos) {$featurecount = $featurecount + $wiggledepth;}
                if ($featureend < $wigglepos) {last}
    	    }
    	    print COUNTSGENE "$featurecount\t$featurename\n";
    	    seek(WIGGLE, $buffer[0], 0);
    	    #seek(WIGGLE, 0, 0);
    	    #alternative to speed up. Insertions sites within overlapping genes will not be counted correctly when using this
    	    #seek(WIGGLE, -length($lengthline), 1);
        }
        close(WIGGLE);
        close(LOCATIONS);
        close(COUNTSGENE);
        
        ## Loess normalize the per-gene read counts. Overwrites the input file (gene_XXX.counts.txt)
        if ($loess eq "Yes") {
            print "Loess normalizing ".$countsgene." for insertion frequency bias\n";
            qx("$Rnew" --vanilla --args $countsgene < "$loessgenescript" 2>>R.error);
        }
    }
}# End pass_alignments


#######################################################################################
# doTAcounts
# This rubroutine prepares unique insertion site counts per gene for essential gene statistics

sub doTAcounts {

    my $featureend;
    my $featurestart;
    my $featureline;
    my $wiggleline;
    my $wigglepos;
    my $featurecount;
    my $featurename;
    my $wiggledepth;
    my @splet;
    my @splyt;
    my $countsgene = "tavsgenes.txt";
    my $wiggle = "tavsta.wiggle";    
    my $wigglepos_remember=0;
    my $wiggledepth_remember=0;
                   
    print "\nCounting number of possible unique insertion site flanking sequences per gene\n";
    #perform pass alignment command, sort command and uniq command, do some magic with sed and tr to make it tab delimited.  

    if ($insertion eq "TA") { 
        print "\nCounting number of unique insertion site flanking sequences in genome.ta before removal of non-uniques \n";
        countlines("genome.ta");
        qx("$perl" "$tavsta" genome.ta > tavsta.txt); ## replaces pass alignment. Much faster and simpler anyway. Confirmed output was the same as previous pass commands
	}
    
    ## This is untouched. Will need to fix if we start using random transposon insertion libraries.
    if ($insertion eq "random") {
    	print "\nCounting number of unique insertion site flanking sequences in genome.ta before removal of non-uniques \n";
        countlines("genome.ta");
        print "\nDetecting non-unique insertion sites\n";
        my $uniqoption; ## Dummy since this option was deleted.
        print "".$pass." -i genome.ta -d genome.fasta -fle ".$minhit." -b ".$uniqoption."  -cpu 4 -s -gff -S 2 \n";
        qx("$pass" -i genome.ta -d genome.fasta -fle "$minhit" -b "$uniqoption"  -cpu 4 -j -g 3 -s -gff -S 2  2>/dev/null > passout_tavsta.txt);
        if ($transposon_end eq "bol") {
            qx("$cat" passout_tavsta.txt |"$grep" "	-	" |"$cut" -f 5 | "$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
            qx("$cat" passout_tavsta.txt |"$grep" "	+	" |"$cut" -f 4 | "$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' >> tavsta.txt);
	    }
        if ($transposon_end eq "eol") {
            qx("$cat" passout_tavsta.txt |"$grep" "	+	" |"$cut" -f 5 | "$sed" 's?^?-?' |"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
            qx("$cat" passout_tavsta.txt |"$grep" "	-	" |"$cut" -f 4 | "$sed" 's?^?-?' |"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' >> tavsta.txt);
	    }
        unlink("passout_tavsta.txt");
        
    	#qx("$pass" -i genome.ta -d genome.fasta -b "$uniqoption" -flc 0 -cpu 4 -gff -S 2 2>/dev/null |"$cut" -f 4 |"$sort" |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
	}

    #awful hack. plan to incorporate the removal of the - sign in the perl script below. I need to sort the file anyway.
    print "Counting number of unique insertion site flanking sequences in genome.ta after removal of non-uniques \n";
    countlines2("tavsta.txt");
    ## Remove header, remove negative signs, and sort by position 
    qx("$cat" tavsta.txt |"$grep" -v names|"$tr" -d '-'|"$sort" -n -k2 >tavsta.wiggle); ## I could do this in tavsta.pl, but will leave it here. Doesn't take long.

    # rapidly add all ta site counts (the wiggle track of the genome) per gene, which positions are annotated in genome.ptt.    
    open (COUNTSGENE, "> $countsgene" ) or die("Can't open file: ".$countsgene."");
    print COUNTSGENE "counts\tnames\n";
    open(LOCATIONS, $ptt) or die("Can't open file: ".$ptt."");
    open(WIGGLE, "< $wiggle" )  or die("Can't open file: ".$wiggle."");
    
    my @buffer = (0) x 10000;
    my $filepos = 0 ;
    
    while (<LOCATIONS>){
        $featureline = $_ ;
        chomp $featureline ;
        @splyt = split '\t', $featureline;
        $featurename =  $splyt[0];
        $featurestart = $splyt[1];
        $featureend = $splyt[2];
        $featurecount=0;
        next if ($featurename eq "!Locus");
        while (<WIGGLE>){
    	    $filepos = tell(WIGGLE);
    	    shift @buffer;
    	    push  @buffer, $filepos;
    	    
    	    $wiggleline = $_ ;
    	    #$lengthline = $_ ;
    	    chomp $wiggleline ;
    	    @splet = split '\t', $wiggleline;
    	    $wigglepos = $splet[1];
    	    $wiggledepth = $splet[0];
            if ($featurestart <= $wigglepos && $featureend >= $wigglepos) {$featurecount = $featurecount + $wiggledepth;}
            if ($featureend < $wigglepos) {last}
            
        }
        print COUNTSGENE "$featurecount\t$featurename\n";
        seek(WIGGLE, $buffer[0], 0);
        #seek(WIGGLE, 0, 0);
    }
    close(WIGGLE);
    close(LOCATIONS);
    close(COUNTSGENE);

    open (TARGETSESSENTIALGENES, '>>essentialgenes_targets.txt');
    print TARGETSESSENTIALGENES "tavsgenes.txt\tcontrol\tcontrol\tinsertionsites\n";
    close (TARGETSESSENTIALGENES);
}# END doTAcounts

#######################################################################################
# doTAannotation
# This rubroutine annotates the insertion site flanking sequences

sub doTAannotation {

    my $featureend;
    my $featurestart;
    my $featureline;
    my $wiggleline;
    my $wigglepos;
    my $featurecount;
    my $featurename;
    my $wiggledepth;
    my $feature3;
    my $feature4;
    my $feature5;
    my $feature6;
    my $feature7;
    my $wigglename;
    my @splet;
    my @splyt;
    my $annotation = "ta.ptt";
    my $wiggle = "tavsta.wiggle3";    
    my $wigglepos_remember=0;
    my $wigglename_remember=0;
    print "\nAnnotating insertion site flanking sequences\n";
    my $lengthline;
    #ugly hack. fix asap
    qx("cat" allinsertions.txt |"$cut" -f 2|"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > allinsertionsuniq.txt);
    qx("$cat" allinsertionsuniq.txt |"$tr" -d '-' >tavsta.wiggle2);
    qx("$paste" tavsta.wiggle2 allinsertionsuniq.txt |"$sort" -n -k 2 >tavsta.wiggle3);
    #qx("$cat" $ptt |"$sort" -k2 -n > genome.ptt.tmp);
    #qx("mv" genome.ptt.tmp $ptt);

    # rapidly annotate the insertion site flanking sequences (the wiggle track of the genome) with genes, which positions are annotated in genome.ptt.
    
    open (ANNOTATION, "> $annotation" )	or die("Can't open file: ".$annotation."");
    print ANNOTATION "counts\tnames\n";
    #open(LOCATIONS, 'genome.ptt') 	or die("Can't open file: genome.ptt"); ## why not $ptt (i.e. truncated.ptt)?
    open(LOCATIONS, $ptt) or die("Can't open file: ".$ptt.""); ## Makes more sense to me. Should only count same sites that would be counted in read files in doTAcounts subroutine.
    open(WIGGLE, "< $wiggle" ) 		or die("Can't open file: ".$wiggle."");
    while (<LOCATIONS>){
        $featureline = $_ ;
        chomp $featureline ;
        @splyt = split '\t', $featureline;
        $featurename =  $splyt[0];
        $featurestart = $splyt[1];
        $featureend = $splyt[2];
        $feature3 = $splyt[3]; 
        $feature4 = $splyt[4];
        $feature5 = $splyt[5];
        $feature6 = $splyt[6];
        $feature7 = $splyt[7];
        $featurecount=0;
        next if ($featurename eq "!Locus");
        while (<WIGGLE>){
            $wiggleline = $_ ;
            $lengthline = $_;
            chomp $wiggleline ;
            @splet = split '\t', $wiggleline;
            $wigglepos = $splet[1];
            $wiggledepth = $splet[0];
            $wigglename = $splet[3];
    	    
    	    if ($featurestart <= $wigglepos && $featureend >= $wigglepos) {
                print ANNOTATION "$wigglename\t$featurename\t$featurestart\t$featureend\t$feature3\t$feature4\t$feature5\t$feature6\t$feature7\n";
    	    }
    	    if ($featureend < $wigglepos) {last}
    	}
        seek(WIGGLE, -length($lengthline), 1);
    }
    close(WIGGLE);
    close(LOCATIONS);
    close(ANNOTATION);
} #END doTAannotation

#######################################################################################
# doStats
# This subroutine will perform the statistics.

sub doStats {
    print "\nPerforming statistical tests (edgeR) on conditionally essential genes\n";
    print "".$Rnew." --vanilla --args $normalization $adjustment $dispersion $smoothing $statmethod gene < ".$edger."\n" ;
    qx("$Rnew" --vanilla --args $normalization $adjustment $dispersion $smoothing $statmethod gene < "$edger" 2>>R.error );
    print "Performing statistical tests (edgeR) on essential genes\n";
    print "".$Rnew." --vanilla --args $normalization $adjustment $dispersion $smoothing $statmethod essentialgenes < ".$edger."\n";
    qx("$Rnew" --vanilla --args $normalization $adjustment $dispersion $smoothing $statmethod essentialgenes < "$edger" 2>>R.error );
    print "Performing statistical tests (edgeR) on insertion site flanking sequences\n";
    print "".$Rnew." --vanilla --args $normalization $adjustment $dispersion $smoothing $statmethod ta < ".$edger."\n"; 
    qx("$Rnew" --vanilla --args $normalization $adjustment $dispersion $smoothing $statmethod ta < "$edger" 2>>R.error );
}#END doStats

#######################################################################################
# truncategenes
# This subroutine will chop of a small part of the 3' end of a gene to filter out matches with transposon insertions that have no effect on function
# it takes the log(1.1) of 10% of the length of the gene and removes that from the end. 
# This empirical value was shown to perform quite well with both long and short genes.

sub truncategenes {
    use POSIX;

    my $pttline;
    my $locusptt;
    my $startptt;
    my $stopptt;
    my $strandptt;
    my $pttvar4;
    my $pttvar5;
    my $pttvar6;
    my $pttvar7;
    my $pttlength;
    my $newstartptt;
    my $newstopptt;
    my @spliit;

    print "Applying gene truncation filter\n";
    open(PTT, 'genome.ptt') or general::error("Can't open file: genome.ptt");
    open (TRUNCATED, "> truncated.ptt");
    while (<PTT>){
        $pttline = $_ ;
        chomp $pttline;
        @spliit = split '\t', $pttline;
        $locusptt = $spliit[0] ;
        $startptt = $spliit[1] ;
        $stopptt  = $spliit[2] ;
        $strandptt = $spliit[3] ;
        $pttvar4 =  $spliit[4] ;
        $pttvar5 = $spliit[5];
        $pttvar6 = $spliit[6];
        $pttvar7 = $spliit[7];
        if ($locusptt eq "!Locus") {
            print TRUNCATED "$locusptt\t$startptt\t$stopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
        }
        next if ($locusptt eq "!Locus");
        $pttlength=$stopptt-$startptt;
        if ($pttlength <20) {
            print TRUNCATED "$locusptt\t$startptt\t$stopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
        }
        if ($pttlength >=20) {
            if ($strandptt eq "-") {
            $newstartptt = $startptt;
            $newstartptt = floor($startptt + (log (0.05 * $pttlength))/(log(1.1)) );
            print TRUNCATED "$locusptt\t$newstartptt\t$stopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
            }
            if ($strandptt eq "+") {
            $newstopptt = $stopptt;
            $newstopptt = floor($stopptt - (log (0.05 * $pttlength))/(log(1.1)) );
            print TRUNCATED "$locusptt\t$startptt\t$newstopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
            }
        }
    }
    close PTT;
    close TRUNCATED;
}#END truncategenes

#subroutine to count number of sequences in fasta file
sub countlines {
    my $filename = shift @_;
    my $lines = 0;
    my $buffer = 0;
    open(FILE, $filename) or die "Can't open `$filename': $!";
    while (<FILE>){ 
	if ($_ =~ m/\>/) {
	    $lines=$lines+1;}
	}
    close FILE;
    print "The number of sequences in $filename is $lines.\n";
}#END countlines

#subroutine to count number of sequences/lines in raw text file
sub countlines2 {

    my $filename = shift @_;
    my $lines = 0;
    my $buffer = 0;
    open(FILE, $filename) or die "Can't open `$filename': $!";
    while (sysread FILE, $buffer, 4096) {
        $lines += ($buffer =~ tr/\n//);
    }
    close FILE;
    my $sequences = $lines ;
    print "The number of sequences in $filename is $sequences.\n";
}#END countlines2


#countwiggle
sub countwiggle {
    my $wigglefile = shift @_;
    my $splet2;
    my $totallines = 0;
    my $totaldepth = 0;
    open(WIGGLEFILE, "< $wigglefile" ) or general::error("Can't open file: ".$wigglefile."");
    while (<WIGGLEFILE>){
        my $wigglefileline = $_ ;
        chomp $wigglefileline ;
        my @splet2 = split '\t', $wigglefileline;
        my $wigglefiledepth = $splet2[0];
        $totaldepth = $totaldepth + $wigglefiledepth;
        $totallines = $totallines + 1;
    }
    print "$totaldepth reads aligned to $totallines insertion site flanking sequences\n";
}#END countwiggle

sub is_path {
    ## Subroutine based on StackOverflow post by Sinan Unur (https://stackoverflow.com/a/8243770)
    my $exe = shift;
    my @path = path;
    my @pathext = ( q{} );
    if ($^O eq 'MSWin32'){
        push @pathext, map { lc } split /;/, $ENV{PATHEXT};
    }
    for my $dir (@path){
        for my $ext (@pathext){
            my $f = catfile $dir, "$exe$ext";
            return ($f) if -x $f;
        }
    }
    return();
}
