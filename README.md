
## Rustynetics

Rustynetics is a high-performance Rust library designed for bioinformatics applications, offering efficient and scalable handling of common genomic file formats. It supports reading and writing of widely used formats such as BAM, FASTQ, FASTA, bigWig, bedGraph, BED, and GFF, making it an essential tool for genomic data processing pipelines.

The library excels in computing coverage tracks, summarizing sequence alignments or read counts across the genome, allowing users to generate coverage profiles over specified the genome. In addition, it offers advanced statistical features, such as the calculation of cross-correlations, which can be used to assess relationships between different genomic datasets, for example, in ChIP-seq or RNA-seq analysis.

One of the library's core strengths is its efficient handling of genomic ranges. It offers a highly optimized data structure for manipulating large genomic intervals, ensuring that operations like querying, merging, or intersecting genomic regions are performed with minimal overhead. Moreover, the library provides sequence containers for FASTA data, motif/PWM utilities, k-mer counting, segmentation import/export, and pretty printing for displaying genomic ranges in human-readable formats.

Designed with performance and usability in mind, this library is ideal for large-scale genomics projects requiring both speed and precision, whether for research in genomics, epigenetics, or other related fields.

## Documentation

Please find the API documentation [here](https://docs.rs/crate/rustynetics/latest).

### Tools

The package contains the following command line tools:

| Tool                     | Description                                                              |
| ------------------------ | ------------------------------------------------------------------------ |
| bam-check-fastq          | check whether all BAM read names are present in a FASTQ file             |
| bam-check-bin            | check bin records of a bam file                                          |
| bam-genome               | print the genome (sequence table) of a bam file                          |
| bam-to-fastq             | reconstruct FASTQ records from a BAM file                                |
| bam-to-bigwig            | convert bam to bigWig (estimate fragment length if required)             |
| bam-view                 | print contents of a bam file                                             |
| bigwig-edit-chrom-names  | rewrite a bigWig with chromosome names transformed by a regex            |
| bigwig-extract           | extract bigWig data for BED regions as a table or bigWig                 |
| bigwig-extract-chroms    | write a new bigWig containing only selected chromosomes                  |
| bigwig-genome            | print the genome (sequence table) of a bigWig file                       |
| bigwig-histogram         | compute a histogram or cumulative histogram over track values            |
| bigwig-info              | print information about a bigWig file                                    |
| bigwig-map               | apply a shared-library mapping function across one or more bigWig tracks |
| bigwig-nil               | re-encode a bigWig track through the Rust implementation                 |
| bigwig-positive          | call joint positive regions across one or more bigWig tracks             |
| bigwig-quantile-normalize| quantile-normalize one bigWig track against a reference                  |
| bigwig-query             | retrieve data from a bigWig file                                         |
| bigwig-query-sequence    | retrieve sequences from a bigWig file                                    |
| bigwig-statistics        | print summary statistics for a bigWig track                              |
| chromhmm-tables-to-bigwig| convert ChromHMM per-chromosome tables to bigWig                         |
| count-kmers              | count or identify k-mers in FASTA sequences or BED regions               |
| draw-genomic-regions     | draw random genomic regions from a genome                                |
| fasta-extract            | extract FASTA subsequences for BED regions                               |
| fasta-unresolved-regions | report unresolved (`N`) FASTA intervals as BED                           |
| gtf-to-bed               | convert GTF records to BED6                                              |
| meme-extract             | extract PWM or PPM motif matrices from MEME or DREME XML                 |
| observed-over-expected-cpg | compute observed/expected CpG scores for regions or whole sequences    |
| pwm-scan-regions         | score genomic regions with one or more PWMs                              |
| pwm-scan-sequences       | scan FASTA sequences with a PWM and export a bigWig track                |
| segmentation-differential| merge and score differential chromatin states across segmentations       |


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

### Read GTF files

```rust

use crate::granges::GRanges;

let granges = GRanges::import_gtf("src/granges_gtf.gtf",
    vec!["gene_id", "gene_num"], // Names of optional fields
    vec!["str"    , "int"     ], // Types of optional fields
    vec![None     , Some("0") ], // Default values, can be an empty vector if omitted
).unwrap();
```
The result is:
```bash
  seqnames ranges         strand |                             source    feature gene_num         gene_id
1 1        [11869, 14409) +      | transcribed_unprocessed_pseudogene       gene        1 ENSG00000223972
2 1        [11870, 14410) +      |               processed_transcript transcript        0 ENSG00000223972
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

### Read and write FASTQ records

```rust
use std::io::{BufReader, BufWriter};
use std::fs::File;

use rustynetics::fastq::{FastqReader, FastqRecord, FastqWriter};

let input = File::open("reads.fastq").unwrap();
let mut reader = FastqReader::new(BufReader::new(input));

let output = File::create("copy.fastq").unwrap();
let mut writer = FastqWriter::new(BufWriter::new(output));

while let Some(record) = reader.read_record().unwrap() {
    writer.write_record(&record).unwrap();
}

let record = FastqRecord::new(
    "read1/1".to_string(),
    "ACGT".to_string(),
    "IIII".to_string(),
).unwrap();

writer.write_record(&record).unwrap();
writer.flush().unwrap();
```

### Reconstruct FASTQ from BAM

`bam-to-fastq` reconstructs FASTQ records from BAM alignments. Reverse-strand reads are emitted in original FASTQ orientation by reverse-complementing the sequence and reversing the qualities.

```bash
# Write a single FASTQ stream
bam-to-fastq reads.bam > reads.fastq

# Split paired-end output into separate files and keep singles
bam-to-fastq reads.bam \
  --output1 reads_R1.fastq \
  --output2 reads_R2.fastq \
  --output-single reads_single.fastq \
  --pair-suffixes
```

If a BAM file does not contain quality scores, use `--fill-missing-quality` to emit synthetic qualities:

```bash
bam-to-fastq reads.bam --fill-missing-quality I > reads.fastq
```

### Sequence, PWM, and k-mer tools

The package also contains utilities for FASTA extraction, motif scanning, and k-mer counting:

```bash
# Extract sequences for BED intervals
fasta-extract genome.fa peaks.bed > peaks.fa

# Report unresolved regions
fasta-unresolved-regions genome.fa > unresolved.bed

# Count k-mers in FASTA sequences
count-kmers 4 6 genome.fa counts.table

# Score genomic regions with one or more PWMs
pwm-scan-regions --input peaks.bed genome.fa motif1.table motif2.table > pwm.table

# Scan complete sequences and export the scores as a bigWig track
pwm-scan-sequences motif.table genome.fa motif.bw
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
(data=(chrom_id=chrY, from=1838100, to=1838200, statistics=(valid=1, min=1.0000, max=1.0000, sum=1.0000, sum_squares=1.0000)), type=fixed)
(data=(chrom_id=chrY, from=1838200, to=1838300, statistics=(valid=1, min=1.0000, max=1.0000, sum=1.0000, sum_squares=1.0000)), type=fixed)
(data=(chrom_id=chrY, from=1838300, to=1838400, statistics=(valid=1, min=0.0000, max=0.0000, sum=0.0000, sum_squares=0.0000)), type=fixed)
(data=(chrom_id=chrY, from=1838400, to=1838500, statistics=(valid=1, min=0.0000, max=0.0000, sum=0.0000, sum_squares=0.0000)), type=fixed)
(data=(chrom_id=chrY, from=1838500, to=1838600, statistics=(valid=1, min=0.0000, max=0.0000, sum=0.0000, sum_squares=0.0000)), type=fixed)
```

### BigWig utilities

```bash
# Apply a custom shared-library function across tracks
bigwig-map mapper.so result.bw track1.bw track2.bw

# Extract selected regions from a bigWig file as a table
bigwig-extract signal.bw regions.bed signal.table

# Keep only selected chromosomes
bigwig-extract-chroms chr1,chr2 signal.bw subset.bw

# Quantile-normalize one bigWig track against another
bigwig-quantile-normalize reference.bw input.bw normalized.bw

# Compute global statistics or a histogram
bigwig-statistics signal.bw
bigwig-histogram --bins 200 signal.bw

# Find regions where multiple tracks are jointly positive
bigwig-positive result.table track1.bw:2.0 track2.bw:2.0

# Preview chromosome renaming without modifying the file
bigwig-edit-chrom-names --dry-run signal.bw '^chr' ''
```

### Motif XML extraction

```bash
# Extract MEME motifs as log-odds matrices
meme-extract motifs.xml motif

# Extract DREME motifs as PPMs in JASPAR format
meme-extract --input-format dreme --output-type ppm --output-format jaspar dreme.xml dreme
```

### Compute coverage tracks from BAM files

We download a BAM file from a ChIP-seq experiment in Homo sapiens A549 with FOXS1 as target (ENCFF504WRM) in addition to the control data (ENCFF739ECZ):
```bash
wget https://www.encodeproject.org/files/ENCFF504WRM/@@download/ENCFF504WRM.bam
wget https://www.encodeproject.org/files/ENCFF739ECZ/@@download/ENCFF739ECZ.bam
```

```rust

use rustynetics::bam_coverage::bam_coverage;

let tracks_treatment = vec!["ENCFF504WRM.bam"];
let tracks_control   = vec!["ENCFF739ECZ.bam"];

// Set fragment length to 0, which means that fragments will not be extended.
// Setting this to None will trigger automatic fragment length estimation
let fraglen_treatment = vec![Some(0)];
let fraglen_control   = vec![Some(0)];

let (track, _treatment_fraglen_estimates, _control_fraglen_estimates) = bam_coverage(
    &tracks_treatment,
    &tracks_control,
    &fraglen_treatment,
    &fraglen_control,
    vec![]
).unwrap();

if let Err(e) = track.export_bigwig("track.bw", vec![]) {
    panic!("{}", e);
}
```
