use gxf2bed::{run, BedType, Config};
use indoc::indoc;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Writes a file to the temporary directory and returns its path.
fn write_temp_file(dir: &Path, name: &str, contents: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, contents).unwrap();
    path
}

/// Converts a small GTF to BED12 and validates coordinates and blocks.
#[test]
fn convert_gtf_to_bed12() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = indoc! {"
        chr1\tsrc\ttranscript\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; transcript_name \"tx1\";
        chr1\tsrc\texon\t100\t150\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; exon_number \"1\";
        chr1\tsrc\texon\t180\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; exon_number \"2\";
        chr2\tsrc\ttranscript\t1000\t1100\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\"; transcript_name \"tx2\";
        chr2\tsrc\texon\t1000\t1050\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\"; exon_number \"1\";
        chr2\tsrc\texon\t1070\t1100\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\"; exon_number \"2\";
    "};
    let input_path = write_temp_file(dir.path(), "input.gtf", gtf.trim());
    let output_path = dir.path().join("output.bed");

    let config = Config {
        input: input_path,
        output: output_path.clone(),
        threads: 2,
        parent_feature: None,
        child_features: None,
        parent_attribute: None,
        child_attribute: None,
        bed_type: BedType::Bed12,
        additional_fields: None,
        chunks: 1024,
        include_non_coding: false,
    };

    run(&config).unwrap();

    let output = std::fs::read_to_string(&output_path).unwrap();
    let mut by_name: HashMap<String, Vec<String>> = HashMap::new();
    for line in output.lines() {
        let fields = line.split('\t').map(|s| s.to_string()).collect::<Vec<_>>();
        assert_eq!(fields.len(), 12);
        by_name.insert(fields[3].clone(), fields);
    }

    let tx1 = by_name.get("tx1").unwrap();
    assert_eq!(tx1[0], "chr1");
    assert_eq!(tx1[1], "99");
    assert_eq!(tx1[2], "200");
    assert_eq!(tx1[5], "+");
    assert_eq!(tx1[9], "2");
    assert_eq!(tx1[10], "51,21,");
    assert_eq!(tx1[11], "0,80,");

    let tx2 = by_name.get("tx2").unwrap();
    assert_eq!(tx2[0], "chr2");
    assert_eq!(tx2[1], "999");
    assert_eq!(tx2[2], "1100");
    assert_eq!(tx2[5], "-");
    assert_eq!(tx2[9], "2");
    assert_eq!(tx2[10], "51,31,");
    assert_eq!(tx2[11], "0,70,");
}

/// Validates that --include-non-coding sets thickStart == thickEnd for non-coding
/// transcripts (those without CDS features) and preserves thick regions for coding ones.
#[test]
fn include_non_coding_sets_thin_rendering() {
    let dir = tempfile::tempdir().unwrap();
    // tx_coding has CDS → thick region should be preserved
    // tx_noncoding has no CDS → thick should become thickStart == thickEnd == chromStart
    let gtf = indoc! {"
        chr1\tsrc\ttranscript\t100\t300\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx_coding\"; transcript_name \"tx_coding\";
        chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx_coding\"; exon_number \"1\";
        chr1\tsrc\texon\t250\t300\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx_coding\"; exon_number \"2\";
        chr1\tsrc\tCDS\t120\t190\t.\t+\t0\tgene_id \"g1\"; transcript_id \"tx_coding\";
        chr1\tsrc\tCDS\t260\t290\t.\t+\t1\tgene_id \"g1\"; transcript_id \"tx_coding\";
        chr2\tsrc\ttranscript\t500\t800\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_noncoding\"; transcript_name \"tx_noncoding\";
        chr2\tsrc\texon\t500\t600\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_noncoding\"; exon_number \"1\";
        chr2\tsrc\texon\t700\t800\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_noncoding\"; exon_number \"2\";
    "};
    let input_path = write_temp_file(dir.path(), "mixed.gtf", gtf.trim());
    let output_path = dir.path().join("mixed_output.bed");

    let config = Config {
        input: input_path,
        output: output_path.clone(),
        threads: 2,
        parent_feature: None,
        child_features: None,
        parent_attribute: None,
        child_attribute: None,
        bed_type: BedType::Bed12,
        additional_fields: None,
        chunks: 1024,
        include_non_coding: true,
    };

    run(&config).unwrap();

    let output = std::fs::read_to_string(&output_path).unwrap();
    let mut by_name: HashMap<String, Vec<String>> = HashMap::new();
    for line in output.lines() {
        let fields = line.split('\t').map(|s| s.to_string()).collect::<Vec<_>>();
        assert_eq!(fields.len(), 12);
        by_name.insert(fields[3].clone(), fields);
    }

    // Coding transcript: thick region should span the CDS
    let coding = by_name.get("tx_coding").unwrap();
    let thick_start: u64 = coding[6].parse().unwrap();
    let thick_end: u64 = coding[7].parse().unwrap();
    assert!(thick_start < thick_end, "Coding transcript should have thick_start < thick_end");
    assert_eq!(thick_start, 119); // 120 - 1 (0-based)
    assert_eq!(thick_end, 290);

    // Non-coding transcript: thickStart == thickEnd (thin rendering)
    let noncoding = by_name.get("tx_noncoding").unwrap();
    let nc_thick_start: u64 = noncoding[6].parse().unwrap();
    let nc_thick_end: u64 = noncoding[7].parse().unwrap();
    assert_eq!(nc_thick_start, nc_thick_end, "Non-coding transcript should have thickStart == thickEnd");
    // thickStart should equal chromStart
    let nc_chrom_start: u64 = noncoding[1].parse().unwrap();
    assert_eq!(nc_thick_start, nc_chrom_start, "Non-coding thickStart should equal chromStart");
    // Should still have 2 exon blocks
    assert_eq!(noncoding[9], "2");
}

/// Validates that without --include-non-coding, non-coding transcripts get
/// thickStart == chromStart and thickEnd == chromEnd (fully thick, default genepred behavior).
#[test]
fn without_include_non_coding_default_thick() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = indoc! {"
        chr2\tsrc\ttranscript\t500\t800\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_nc\"; transcript_name \"tx_nc\";
        chr2\tsrc\texon\t500\t600\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_nc\"; exon_number \"1\";
        chr2\tsrc\texon\t700\t800\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_nc\"; exon_number \"2\";
    "};
    let input_path = write_temp_file(dir.path(), "nc_default.gtf", gtf.trim());
    let output_path = dir.path().join("nc_default_output.bed");

    let config = Config {
        input: input_path,
        output: output_path.clone(),
        threads: 2,
        parent_feature: None,
        child_features: None,
        parent_attribute: None,
        child_attribute: None,
        bed_type: BedType::Bed12,
        additional_fields: None,
        chunks: 1024,
        include_non_coding: false,
    };

    run(&config).unwrap();

    let output = std::fs::read_to_string(&output_path).unwrap();
    let line = output.lines().next().unwrap();
    let fields = line.split('\t').collect::<Vec<_>>();
    // Without flag: thickStart == chromStart, thickEnd == chromEnd (genepred default)
    assert_eq!(fields[6], fields[1], "Default: thickStart should equal chromStart");
    assert_eq!(fields[7], fields[2], "Default: thickEnd should equal chromEnd");
}
