# rust-bedtools-pangenome-merge
  - rust bedtools pangenome merge for multimerge.
  - anchor alignment approach to pangenome merge, sort the range and then intsert into different structs.
  - no ned to sort the file. 

 ```
 cargo build
 ```
 ```
  ./target/debug/rust-bedtools-pangenome-merge -h
  Usage: rust-bedtools-pangenome-merge <BED1> <BED2> <FASTA>

  Arguments:
  <BED1>   please provide the path to the first alignment file
  <BED2>   please provide the reference fasta file
  <FASTA>  please provide the path to the reference fasta file

  Options:
  -h, --help     Print help
  -V, --version  Print version
 ```

 Gaurav Sablok
