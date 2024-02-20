# Essentials_local

### Analysis pipeline for transposon insertion sequencing experimental data.

This is a command line implementation that runs on your local computer rather than the original, now defunct [web-based version](http://bamics2.cmbi.ru.nl/websoftware/essentials/essentials_start.php
) developed by [Aldert Zomer et al](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043012). This version was developed using the Essentials source code available [here](https://trac.nbic.nl/essentials/browser).  

Original ESSENTIALS software was distrubted under [GNU Affero General Public License](https://trac.nbic.nl/essentials/wiki/license)

**Some differences from original Essentials software:**

* Command line only. No web browser implementation (yet?).
* This version does not perform read demultiplexing, trimming, or alignment. Inputs to essentials\_local are wiggle files listing sequence positions and read counts. To peform read processing from fastq read files to generate input for essentials\_local, you can / shoud use the scripts available at [TIS_tools](https://github.com/egonozer/TIS_tools) depending on what sequencing library protocol you used.
* Currently only TA insertion libraries can be used in this version of essentials\_local. Random insertion libraries are not yet implemented.
* We added some new figure outputs to the pipeline. In addition to the PCA plots, density plots, and fold-change vs signal plots output by original Essentials, this version also outputs MDS and volcano plots, both of which we've found useful in our analyses, as well as plots of read counts per insertion site (TA) and per gene 
before and after loess normalization. 
* Removed multi-library support from the Cox-Reid (CR) analysis option. It was confusingly implemented (to me anyway) and very prone to error if the libraries were not perfectly defined in the input. Would not be too tough to re-integrate if needed, but we have not needed it. 
* Updated to use most recent generation of EdgeR v3.x+ (currently v4.x)
* Removed dependency on 'pass' aligner for identifying unique insertion sites, replaced with Perl script.
* Other changes can be found in the CHANGELOG.txt


**Requirements:** 
   
* [R](https://cloud.r-project.org/) (I have version 4.1.1, but should work on v3 and up)  
* Perl  

**Installation:**  

* Install packages (edgeR, EnhancedVolcano, zoo) in R: 

```r
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')        
BiocManager::install('edgeR')
BiocManager::install('EnhancedVolcano')
install.packages('zoo')
```

* If you want the option to automatically compress the output folder, you will need to have the perl module [Archive::Zip](https://metacpan.org/pod/Archive::Zip) installed. This is (very) optional and the default setting for essentials\_local is no compression. 
 

**Usage:**

`perl essentials.pl -c <configfile.txt> -g <reference.gbk>` 

```
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
  -l or --libsize   Expected library size. Enter 0 to use all sites in each 
                    input wiggle file (default: 0)
  --rdleng          Expected length of reads used for alignment. Options:
                    24 = Boll protocol (default)
                    17 = Goodman protocol
  --insert          Transposon insertion site type. Options are:
                    'TA' (default) or 'random' (NONFUNCTIONAL)                    
  --full            Use full gene lengths (default: use 5' truncated genes)
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
  --adjust          p-value adjustment. Options are:
                    'BH' (default), 'holm', 'hochberg', 'hommel', 'bonferroni',
                    'BY', or 'none'
  --zip             Zip results files. Requires Archive::Zip perl module to be
                    installed. (default: results are not zipped)

  --settings        A settings file containing some or all of the above
                    settings as was used in previous versions of this program.
                    File settings will override any command line settings.
```                    

Questions? Open an issue on Github or send me an email.
