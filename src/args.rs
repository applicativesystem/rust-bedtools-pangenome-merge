use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]
pub struct PangenomeMergeArgs {
    /// please provide the path to the first alignment file
    pub bed1: String,
    /// please provide the reference fasta file
    pub bed2: String,
    /// please provide the path to the reference fasta file
    pub fasta: String,
}
