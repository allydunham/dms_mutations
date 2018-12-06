# Deep Mutagenesis Experiment File Type (given .dm extension)

This is the file type used to be able to work with the different studies consistently in this project.
It is not meant to be completely comprehensive or externally useful, but to aid automation of data processing
The basic format is a glorified tsv file with a series of meta data header lines

## Current notes/todos
* Might be good to rename ref\_seq to something indicating its the protein seq
    * make more general so that arbitatry multi-line fields are supported
    * Use this to add full gene ref seq to Araya et al. 2012

## Version 1.0
Initial version. Current rules:
* Meta lines marked by "#"
* Final meta line/table header started by "?" instead of "#"

* Meta lines for a key value pair with some information about the experiment
    * Key:value are separated by a :
    * Spaces are allowed in value but nowhere else
    * e.g. "#gene\_name:BRCA1"
    * The first meta line gives the file type version (although the structure means this coming first is only convention)

* The refference protein sequence is started by a special key - #ref\_seq: with the following lines containing sequence
    * Sequence lines start with "#+" and contain no spaces or ":"
    * Sequence can be split over multiple lines to aid readability
    * This structure could easily be extended for other multi-line values (DNA sequence for instance) using a different key

* The data header is followed by a tsv style table containing the actual data
    * Each line contains information on one tested variant protein
    * required columns are "variants", "score" and "raw\_score"
    * Other columns from the original data can optionally be included
    * The variants column contains a comma separated list of hgvs style variant labels, using single AA code (e.g. p.V1A)

### Version 1.1
* ref\_seq meta line now must contain a value giving the number of sequence lines to follow (e.g. #ref\_seq:3)
