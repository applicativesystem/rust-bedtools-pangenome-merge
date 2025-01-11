mod args;
mod merge;
use crate::merge::Fasta;
use crate::merge::MergeBed;
use crate::merge::MergeBed1;
use crate::merge::MergeBed2;
use crate::merge::MergeFasta;
use crate::merge::UnmergedBed;
use crate::merge::UnmergedFasta;
use args::PangenomeMergeArgs;
use clap::Parser;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
*Author Gaurav Sablok
*Universitat Potsdam and SLB Potsdam
*Date 2025-1-11

 bedtools merge for the pangenome.By implementing the two struct approach,
 there is not need to sort the bed file for the genome. You can directly
 pass the bed files for the genome.

*/

fn main() {
    let args = PangenomeMergeArgs::parse();
    let output = pangenome_merge(&args.bed1, &args.bed2, &args.fasta).unwrap();
    println!("Results have been written:{}", output);
}

fn pangenome_merge(path1: &str, path2: &str, path3: &str) -> Result<String, Box<dyn Error>> {
    let bed1_open = File::open(path1).expect("file not present");
    let bed2_open = File::open(path2).expect("file not present");
    let bed1_read = BufReader::new(bed1_open);
    let bed2_read = BufReader::new(bed2_open);

    let mut bed1: Vec<MergeBed1> = Vec::new();
    let mut bed2: Vec<MergeBed2> = Vec::new();

    for i in bed1_read.lines() {
        let line = i.expect("line not present");
        let linevec = line.split("\t").collect::<Vec<_>>();
        bed1.push(MergeBed1 {
            bed1: linevec[0].to_string(),
            start1: linevec[1].parse::<usize>().unwrap(),
            end1: linevec[2].parse::<usize>().unwrap(),
            bedfragment1: linevec[3].to_string(),
        });
    }

    for i in bed2_read.lines() {
        let line = i.expect("line not present");
        let linevec = line.split("\t").collect::<Vec<_>>();
        bed2.push(MergeBed2 {
            bed2: linevec[0].to_string(),
            start2: linevec[1].parse::<usize>().unwrap(),
            end2: linevec[2].parse::<usize>().unwrap(),
            bedfragment2: linevec[3].to_string(),
        });
    }

    let fastaslice: Vec<Fasta> = fasta_estimate(path3).unwrap();

    let mut merge_bed: Vec<MergeBed> = Vec::new();
    let mut unmerged_bed1: Vec<UnmergedBed> = Vec::new();
    let mut unmerged_bed2: Vec<UnmergedBed> = Vec::new();

    for i in bed1.iter() {
        for j in bed2.iter() {
            if j.start2 > i.start1 && j.end2 > i.end1 && j.start2 < i.end1 {
                let mut bedholditer: Vec<(String, String)> = vec![];
                let mut fragmentmerge1: Vec<(usize, usize)> = vec![];
                let mut fragmentmerge2: Vec<(usize, usize)> = vec![];
                bedholditer.push((i.bedfragment1.clone(), j.bedfragment2.clone()));
                fragmentmerge1.push((i.start1, i.end1));
                fragmentmerge2.push((j.start2, j.end2));
                merge_bed.push(MergeBed {
                    bed: i.bed1.clone(),
                    start1: i.start1,
                    end1: j.end2,
                    bedfragment: bedholditer,
                    mergedfragment1: fragmentmerge1,
                    mergedfragement2: fragmentmerge2,
                })
            } else if j.start2 > i.start1 && j.end2 > i.end1 {
                unmerged_bed1.push(UnmergedBed {
                    bed: i.bed1.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    bedfragment: i.bedfragment1.clone(),
                });
                unmerged_bed2.push(UnmergedBed {
                    bed: j.bed2.clone(),
                    start1: j.start2,
                    end1: j.end2,
                    bedfragment: j.bedfragment2.clone(),
                })
            }
        }
    }

    let mut merged_bed_write = File::create("merged-bed.txt").expect("file not present");
    for i in merge_bed.iter() {
        writeln!(
            merged_bed_write,
            "{}\t{}\t{:?}",
            i.start1, i.end1, i.bedfragment
        )
        .expect("file not present");
    }
    for val in unmerged_bed1.iter() {
        writeln!(
            merged_bed_write,
            "{}\t{}\t{}",
            val.start1, val.end1, val.bedfragment
        )
        .expect("file not present");
    }
    for gen in unmerged_bed2.iter() {
        writeln!(
            merged_bed_write,
            "{}\t{}\t{}",
            gen.start1, gen.end1, gen.bedfragment
        )
        .expect("file not present");
    }

    let mut mergedfasta: Vec<MergeFasta> = Vec::new();
    let mut unmergedfasta1: Vec<UnmergedFasta> = Vec::new();
    let mut unmergedfasta2: Vec<UnmergedFasta> = Vec::new();

    let mut mergedfasta_write = File::create("mergedfasta").expect("file not present");
    let mut unmergedfasta_write1 = File::create("unmerged1.fasta").expect("file not present");
    let mut unmergedfasta_write2 = File::create("unmerged2.fasta").expect("file not present");

    for i in merge_bed.iter() {
        for j in fastaslice.iter() {
            if i.bed == j.header {
                mergedfasta.push(MergeFasta {
                    bed: i.bed.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    mergedfragment: j.sequence[i.start1..i.end1].to_string(),
                })
            }
        }
    }

    for i in unmerged_bed1.iter() {
        for j in fastaslice.iter() {
            if i.bed == j.header {
                unmergedfasta1.push(UnmergedFasta {
                    bed: i.bed.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    seq: j.sequence[i.start1..i.end1].to_string(),
                })
            }
        }
    }

    for i in unmerged_bed2.iter() {
        for j in fastaslice.iter() {
            if i.bed == j.header {
                unmergedfasta2.push(UnmergedFasta {
                    bed: i.bed.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    seq: j.sequence[i.start1..i.end1].to_string(),
                })
            }
        }
    }

    for i in mergedfasta.iter() {
        writeln!(
            mergedfasta_write,
            "{}\t{}\t{}\t{}",
            i.bed, i.start1, i.end1, i.mergedfragment
        )
        .expect("line not present");
    }

    for i in unmergedfasta1.iter() {
        writeln!(
            unmergedfasta_write1,
            "{}\t{}\t{}\t{}",
            i.bed, i.start1, i.end1, i.seq,
        )
        .expect("line not present");
    }

    for i in unmergedfasta2.iter() {
        writeln!(
            unmergedfasta_write2,
            "{}\t{}\t{}\t{}",
            i.bed, i.start1, i.end1, i.seq,
        )
        .expect("line not present");
    }

    Ok("pangenome intersect results have been written".to_string())
}

fn fasta_estimate(path: &str) -> Result<Vec<Fasta>, Box<dyn Error>> {
    let fastaopen = File::open(path).expect("file not present");
    let fastaread = BufReader::new(fastaopen);
    let mut fastaholder: Vec<Fasta> = Vec::new();
    let mut fastaheader: Vec<String> = Vec::new();
    let mut fastasequence: Vec<String> = Vec::new();
    for i in fastaread.lines() {
        let line = i.expect("line not present");
        if line.starts_with(">") {
            fastaheader.push(line.replace(">", ""));
        } else {
            fastasequence.push(line);
        }
    }

    for i in 0..fastaheader.len() {
        fastaholder.push(Fasta {
            header: fastaheader[i].clone(),
            sequence: fastasequence[i].clone(),
        })
    }

    Ok(fastaholder)
}
