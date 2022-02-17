//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! cargo-features = ["edition2021"]
//! [dependencies]
//! indexmap = "1.8"
//! noodles = { version = "0.18.0", features = ["bam", "sam", "bgzf"] }
//! ```

use indexmap::IndexMap;
use noodles::bam::{Reader, Record, Writer};
use noodles::bgzf::writer;
use noodles::sam::{header::ReferenceSequence, record::data::field::Tag, Header};
use std::collections::HashSet;
use std::error::Error;
use std::ffi::CString;
use std::fs::File;
use std::i32;
use std::io::BufWriter;
use std::str::FromStr;


fn main() -> Result<(), Box<dyn Error>> {
    snakemake.redirect_stderr(&snakemake.log[0])?;
    let mut input = File::open(&snakemake.input[0]).map(Reader::new)?;
    let header = input.read_header()?.parse()?;
    let reference_sequences = input.read_reference_sequences()?;

    let mut primerless_writer = build_writer(&snakemake.output.primerless, &header, &reference_sequences)?;
    let mut primer_writer = build_writer(&snakemake.output.primers, &header, &reference_sequences)?;

    let ra_tag: Tag = Tag::from_str("ra")?;
    let mut primary_records = HashSet::new();
    for result in input.records() {
        let record = result?;

        let data = record.data();
        match data.get(ra_tag) {
            Some(Ok(_)) => {
                let idx = i32::from(record.reference_sequence_id().unwrap());
                let chr = reference_sequences
                    .get_index(idx as usize)
                    .unwrap()
                    .0
                    .as_bytes();
                let pos: i32 = i32::from(record.position().unwrap());
                primary_records.insert((chr, pos, record.read_name().unwrap().to_owned()));
            }
            _ => continue,
        }
    }
    let mut input = File::open(&snakemake.input[0]).map(Reader::new)?;
    input.read_header()?;
    input.read_reference_sequences()?;
    let ma_tag = Tag::from_str("ma")?;
    for result in input.records() {
        let record = result?;
        let data = record.data();
        if data.get(ra_tag).is_some()
            || data.get(ma_tag).is_some()
            || is_secondary_alignment(&record, &primary_records)?
        {
            primer_writer.write_record(&record)?;
        } else {
            primerless_writer.write_record(&record)?
        }
    }
    Ok(())
}

fn is_secondary_alignment(
    record: &Record,
    primary_records: &HashSet<(&[u8], i32, CString)>,
) -> Result<bool, Box<dyn Error>> {
    let data = record.data();
    match data.get(Tag::OtherAlignments) {
        Some(Ok(sa_entry)) => {
            let split_tag = sa_entry
                .value()
                .as_str()
                .unwrap()
                .split(',')
                .collect::<Vec<&str>>();
            let chrom = split_tag[0].as_bytes();
            let pos = split_tag[1].parse::<i32>()?;
            return Ok(primary_records.contains(&(
                chrom,
                pos,
                record.read_name().unwrap().to_owned(),
            )));
        }
        _ => Ok(false),
    }
}

fn build_writer(
    file_path: &str,
    header: &Header,
    reference_sequences: &IndexMap<String, ReferenceSequence>,
) -> Result<Writer<writer::Writer<BufWriter<File>>>, Box<dyn Error>> {
    let mut writer = std::fs::File::create(file_path)
        .map(BufWriter::new)
        .map(Writer::new)?;
    writer.write_header(header)?;
    writer.write_reference_sequences(reference_sequences)?;
    Ok(writer)
}