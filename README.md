
## Rustynetics

Rustynetics is a high-performance Rust library designed for bioinformatics applications, offering efficient and scalable handling of common genomic file formats. It supports reading and writing of widely used formats such as bigWig, bedGraph, BED, and GFF, making it an essential tool for genomic data processing pipelines.

The library excels in computing coverage tracks, summarizing sequence alignments or read counts across the genome, allowing users to generate coverage profiles over specified the genome. In addition, it offers advanced statistical features, such as the calculation of cross-correlations, which can be used to assess relationships between different genomic datasets, for example, in ChIP-seq or RNA-seq analysis.

One of the library's core strengths is its efficient handling of genomic ranges. It offers a highly optimized data structure for manipulating large genomic intervals, ensuring that operations like querying, merging, or intersecting genomic regions are performed with minimal overhead. Moreover, the library provides a pretty print feature for easily displaying genomic ranges in human-readable formats, facilitating better visualization and interpretation of complex data.

Designed with performance and usability in mind, this library is ideal for large-scale genomics projects requiring both speed and precision, whether for research in genomics, epigenetics, or other related fields.

## Documentation

Please find the API documentation [here](https://docs.rs/rustynetics/latest/rustynetics/).

## Examples

### Import genes from UCSC

```rust
use crate::genes::Genes;

// Import from local file
if let Ok(genes) = Genes::import_genes("data/hg19.knownGene.txt.gz") {

    println!("{}", genes);
}
// Retrieve from USCS server
if let Ok(genes) = Genes::import_genes_from_ucsc("hg19", "knownGene") {

    println!("{}", genes);
}
```
The result is:
```bash
       seqnames ranges                 strand |               names                  cds
     1 chr1     [    11868,     14409) +      | ENST00000456328.2_1       [11868, 11868)
     2 chr1     [    29553,     31097) +      | ENST00000473358.1_5       [29553, 29553)
     3 chr1     [    30266,     31109) +      | ENST00000469289.1_1       [30266, 30266)
     4 chr1     [    34553,     36081) -      | ENST00000417324.1_4       [34553, 34553)
     5 chr1     [    35244,     36073) -      | ENST00000461467.1_3       [35244, 35244)
       ...      ...                           |                 ...                  ...
254531 chrY     [ 59161253,  59162245) -      | ENST00000711258.1_1 [59161253, 59161253)
254532 chrY     [ 59208304,  59208554) +      | ENST00000711259.1_1 [59208304, 59208304)
254533 chrY     [ 59311662,  59311996) -      | ENST00000711266.1_1 [59311662, 59311662)
254534 chrY     [ 59318040,  59318920) -      | ENST00000711267.1_1 [59318040, 59318040)
254535 chrY     [ 59358334,  59360548) -      | ENST00000711270.1_1 [59358334, 59358334)
```

### Read a BAM file into a GRanges object

```rust
use crate::bam::BamReaderOptions;
use crate::granges::GRanges;

let mut options = BamReaderOptions::default();

options.read_cigar = true;
options.read_qual  = true;

if let Ok(granges) = GRanges::import_bam_single_end("tests/test_bam_2.bam", Some(options)) {
    println!("{}", granges);
}
```
The result is:
```bash
     seqnames     ranges                 strand | flag mapq cigar                                                qual
   1 chr15        [ 25969791,  25969842) +      |   99   60   51M @@@DFFFFHHHHGBHIBHHHGGGIHIEEHEIIIIIIGCHGHIGIGIIIIHH
   2 chr15        [ 25969837,  25969888) -      |  147   60   51M GJJIIIIHIHDIHIIHHEGEEGJIIHFHIHCIHHGEIDHHDDHFDFFD@C@
   3 chr1         [175925088, 175925139) -      |  153    0   51M IIIIIJJJIJJJJJJJJIJGIJIJHJJJJJJJIIJJJJHHHHHFFFFFBCC
   4 chrX         [ 71582197,  71582248) -      |   83   60   51M GGDDIIIGJIJJJJJJJJJHGEHGJJJJIHDEIIGIJJGHHFHFFFFFCC@
   5 chrX         [ 71581965,  71582016) +      |  163   60   51M @CCFFDFFHHDHHJJJIGCHGIGIGIGJJJIGCGCHBFGDBFGFGIJIJGC
     ...          ...                           |  ...  ...   ...                                                 ...
4887 chr11        [  9074777,   9074828) +      |  163   29   51M <@:B;DDDFH:CC>CFEAADFFFCDFHIEHIHJEGGEHIJJIIDGGIGHII
4888 chr7         [  3303179,   3303230) -      |   83   60   51M HIHH@GIIHGHGHCJHGJIIIIIJJJJIJJIIIIIIJJHHGHHFFFFFCCC
4889 chr7         [  3303050,   3303101) +      |  163   60   51M <@<DADADAAFFFC@>DGEHIICEGH@HCCEGHCCEBGGGFG:BFCGGGBB
4890 chr11        [  4737838,   4737889) -      |   83   60   51M DB9;HCD?D??:?:):)CCA<C2:@HFAHEEHF@<?<?:ACADB;:BB1@?
4891 chr11        [  4737786,   4737837) +      |  163   60   51M @@<DDBDDFD+C?A:1CFDHBFHC<?F9+CGGI:49CCGFACE99?DC990
```

### Reading BigWig files

BigWig files contain data in a binary format optimized for fast random access. In addition to the raw data, bigWig files typically contain several zoom levels for which the data has been summarized. The BigWigReader class allows to query data and it automatically selects an appropriate zoom level for the given binsize:
```rust
let seqname = "chrY"; // (can be a regular expression)
let from    = 1838100;
let to      = 1838600;
let binsize = 100;

// The reader accepts either a local file or a file
// hosted on a HTTP server
if let Ok(mut reader) = BigWigFile::new_reader("tests/test_bigwig_2.bw") {

    for item in reader.query(seqname, from, to, binsize) {
        println!("{}", item.unwrap());
    }

}
```
The result is:
```bash
(data=(chrom_id=5, from=1838100, to=1838200, statistics=(valid=1, min=1.0000, max=1.0000, sum=1.0000, sum_squares=1.0000)), type=3)
(data=(chrom_id=5, from=1838200, to=1838300, statistics=(valid=1, min=1.0000, max=1.0000, sum=1.0000, sum_squares=1.0000)), type=3)
(data=(chrom_id=5, from=1838300, to=1838400, statistics=(valid=1, min=0.0000, max=0.0000, sum=0.0000, sum_squares=0.0000)), type=3)
(data=(chrom_id=5, from=1838400, to=1838500, statistics=(valid=1, min=0.0000, max=0.0000, sum=0.0000, sum_squares=0.0000)), type=3)
(data=(chrom_id=5, from=1838500, to=1838600, statistics=(valid=1, min=0.0000, max=0.0000, sum=0.0000, sum_squares=0.0000)), type=3)
```