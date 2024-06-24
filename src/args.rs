use clap::Parser;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct POACli {
    /// Input file path
    #[arg(short, long)]
    pub input: String,

    /// Input file path
    #[arg(short, long)]
    pub output: String,

    /// Enable HTML output
    #[arg(long)]
    pub html: bool,

    /// Enable graph output
    #[arg(long)]
    pub graph: bool,

    /// Display intermediate alignments
    #[arg(long)]
    pub debug: bool,
}
