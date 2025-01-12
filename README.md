# rust-bedtools-pangenome-multi-merge
  - rust bedtools pangenome merge for multimerge.
  - anchor alignment approach to pangenome merge, sort the range and then intsert into different structs.
  - no need to sort the file.
  - general note: Incase of Golang and RUST, please see the last commit message and if it says compiled binary then it is completed or else still in development version. 

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
