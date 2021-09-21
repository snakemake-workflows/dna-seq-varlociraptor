#!/usr/bin/env rust-script
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! rust-htslib = "0.38"
//! ```

use std::error::Error;
use rust_htslib::bam::header::Header;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use std::collections::HashSet;


fn main() -> Result<(), Box<dyn Error>> {
    let mut input = bam::Reader::from_path(&snakemake.input[0])?;
    let header = input.header().to_owned();

    let mut primerless_writer = bam::Writer::from_path(
        &snakemake.output.primerless,
        &Header::from_template(input.header()),
        bam::Format::Bam,
    )?;
    let mut primer_writer = bam::Writer::from_path(
        &snakemake.output.primers,
        &Header::from_template(input.header()),
        bam::Format::Bam,
    )?;
    let mut primary_records = HashSet::new();
    for result in input.records() {
        let record = result?;
        if record.aux(b"ra").is_ok() {
            primary_records.insert((
                header.tid2name(record.tid() as u32),
                record.pos(),
                record.qname().to_vec(),
            ));
        }
    }
    let mut input = bam::Reader::from_path(&snakemake.input[0])?;
    for result in input.records() {
        let record = result?;
        if record.aux(b"ra").is_ok()
            || record.aux(b"ma").is_ok()
            || is_secondary_alignment(&record, &primary_records)?
        {
            primer_writer.write(&record)?
        } else {
            primerless_writer.write(&record)?
        }
    }
    Ok(())
}

fn is_secondary_alignment(
    record: &bam::Record,
    primary_records: &HashSet<(&[u8], i64, Vec<u8>)>,
) -> Result<bool, Box<dyn Error>> {
    if let Ok(Aux::String(sa_tag)) = record.aux(b"SA") {
        let split_tag = sa_tag.split(',').collect::<Vec<&str>>();
        let chrom = split_tag[0].as_bytes();
        let pos = split_tag[1].parse::<i64>()? -1;
        return Ok(primary_records.contains(&(chrom, pos, record.qname().to_vec())));
    }
    Ok(false)
}
