
## Rustynetics

Rustynetics is a high-performance Rust library designed for bioinformatics applications, offering efficient and scalable handling of common genomic file formats. It supports reading and writing of widely used formats such as bigWig, bedGraph, BED, and GFF, making it an essential tool for genomic data processing pipelines.

The library excels in computing coverage tracks, summarizing sequence alignments or read counts across the genome, allowing users to generate coverage profiles over specified the genome. In addition, it offers advanced statistical features, such as the calculation of cross-correlations, which can be used to assess relationships between different genomic datasets, for example, in ChIP-seq or RNA-seq analysis.

One of the library's core strengths is its efficient handling of genomic ranges. It offers a highly optimized data structure for manipulating large genomic intervals, ensuring that operations like querying, merging, or intersecting genomic regions are performed with minimal overhead. Moreover, the library provides a pretty print feature for easily displaying genomic ranges in human-readable formats, facilitating better visualization and interpretation of complex data.

Designed with performance and usability in mind, this library is ideal for large-scale genomics projects requiring both speed and precision, whether for research in genomics, epigenetics, or other related fields.

## Documentation

Please find the API documentation [here](https://docs.rs/rustynetics/latest/rustynetics/).

## Examples

### Read a BAM file into a GRanges object

```rust
use crate::bam::BamReaderOptions;
use crate::granges::GRanges;

let mut granges = GRanges::default();
let mut options = BamReaderOptions::default();

options.read_cigar = true;
options.read_qual  = true;

if Ok(_) = granges.import_bam_single_end("tests/test_bam_2.bam", Some(options)) {
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
