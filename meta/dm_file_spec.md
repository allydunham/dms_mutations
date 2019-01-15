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

### Version 1.2
* Added 'pdb\_id' as a required meta data entry (will be added as NA when the file is read if not present)
* Changed 'pubmed\_id' meta entry to 'pmid' to allow all other IDs (uniprot, pdb, etc.) to be easily grouped as '\_id' fields while 'pmid' remains with the paper info.

## Version 2.0
* Renamed 'ref\_seq' to 'aa\_seq' to allow multiple sequence entries
* Added generic meta data entries for split strings (e.g. seqs) and lists (e.g. list of pdb IDs)
    * \#\+ indicates a split string entry, with the leading line giving the name and number of lines to read (name:N)
        * This is currently applied to fields ending in '\_seq', which as most likely to be long
    * \#\* indicates a list entry, with the leading line in the same format this is read into an array like structure rather than concatenated
* pdb\_id is a list entry of pdb(:chain) for each related pdb entry

#### Version 2.0.1
* Changed version number scheme so as not to be confused for a number
* corresponds to major.minor.clarification type changes
* Major changes are breaking

### Version 2.1
* Added optional third field to pdb\_id entries, giving the offset relative to the uniprot sequence, to allow better conversion
* New field looks like ID:Chain(:Offset)
