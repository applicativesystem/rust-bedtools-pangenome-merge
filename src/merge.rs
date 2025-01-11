/*
 * holding all the structs in the separate files so that they
 * can be easily called as a reference call in the result.
 *
 *
 * */
#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct MergeBed1 {
    pub bed1: String,
    pub start1: usize,
    pub end1: usize,
    pub bedfragment1: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct MergeBed2 {
    pub bed2: String,
    pub start2: usize,
    pub end2: usize,
    pub bedfragment2: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct MergeBed {
    pub bed: String,
    pub start1: usize,
    pub end1: usize,
    pub mergedfragment1: Vec<(usize, usize)>,
    pub mergedfragement2: Vec<(usize, usize)>,
    pub bedfragment: Vec<(String, String)>,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct UnmergedBed {
    pub bed: String,
    pub start1: usize,
    pub end1: usize,
    pub bedfragment: String,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Fasta {
    pub header: String,
    pub sequence: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct MergeFasta {
    pub bed: String,
    pub start1: usize,
    pub end1: usize,
    pub mergedfragment: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct UnmergedFasta {
   pub bed: String,
   pub start1: usize,
   pub end1: usize,
   pub seq: String,

}
