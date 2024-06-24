use alignment::SeqGraphAlignment;
use clap::Parser;
use consensus::Consensus;
use graph::POAGraph;
use std::{
    cmp::max,
    fmt::Write,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write as BufWrite},
};
mod alignment;
mod args;
mod consensus;
mod graph;

fn get_format(path: &str) -> char {
    if path.to_lowercase().ends_with(".csv") {
        ','
    } else if path.to_lowercase().ends_with(".tsv") {
        '\t'
    } else if path.to_lowercase().ends_with(".fasta") {
        '0'
    } else {
        '1'
    }
}

fn main() -> Result<(), String> {
    let args = args::POACli::parse();
    let mut records = vec![];
    let sep = get_format(&args.input);
    match sep {
        ',' | '\t' => {
            let file = File::open(&args.input)
                .map_err(|_| format!("Unable to open file: {}", &args.input))?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line.map_err(|_| "IO Error".to_string())?;
                let chunks: Vec<String> = line.split(sep).map(|s| s.to_string()).collect();
                records.push((chunks[0].clone(), chunks[1..].to_vec()));
            }
        }
        _ => {}
    }

    let mut poa = POAGraph::new(records[0].0.clone(), records[0].1.clone());

    for (label, seq) in &records[1..] {
        let aln = SeqGraphAlignment::align_seq_to_graph(label.clone(), seq.clone(), &poa.graph);
        if args.debug {
            let (width, s, m, g) = aln.get_string(&poa.graph);
            let width = width + 2;
            println!(
                "Graph    : {}",
                g.iter().fold(String::new(), |mut out, s| {
                    let _ = write!(out, "{s:^width$}");
                    out
                })
            );
            println!(
                "Match    : {}",
                m.iter().fold(String::new(), |mut out, s| {
                    let _ = write!(out, "{s:^width$}");
                    out
                })
            );
            println!(
                "Alignment: {}\n",
                s.iter().fold(String::new(), |mut out, s| {
                    let _ = write!(out, "{s:^width$}");
                    out
                })
            );
        }
        poa.add_alignment(aln);
    }

    if args.graph {
        let graph = poa.get_string();
        let file = File::create(format!("{}.graph", args.output))
            .map_err(|_| format!("Unable to create file: {}", &args.output))?;
        let mut writer = BufWriter::new(file);
        writer
            .write_all(graph.as_bytes())
            .map_err(|_| "IO Error".to_string())?;
    }

    let item_width = poa.width + 2;
    let label_width = poa.labels.iter().rfold(0usize, |acc, b| max(acc, b.len()));
    let con = Consensus::new(poa.graph, poa.start_indices, poa.labels);
    let file = File::create(&args.output)
        .map_err(|_| format!("Unable to create file: {}", &args.output))?;
    let mut writer = BufWriter::new(file);

    for (label, seq) in con.compute() {
        let padded_seq: String = seq.into_iter().fold(String::new(), |mut output, item| {
            let _ = write!(output, "{item:^item_width$}");
            output
        });
        writer
            .write_all(format!("{label:^label_width$} {padded_seq}\n").as_bytes())
            .map_err(|_| "IO Error".to_string())?;
    }

    Ok(())
}
