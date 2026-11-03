fn main() {
    if let Err(error) = pathallele::cli::run() {
        eprintln!("PathAllele: {error}");
        std::process::exit(2);
    }
}
