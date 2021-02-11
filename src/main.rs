#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate phasst_lib;
extern crate disjoint_set;
extern crate statrs;
extern crate bio;

use bio::io::fasta;
use bio::utils::Text;
use bio::io::fasta::Record;
use std::path::Path;


use disjoint_set::DisjointSet;

use hashbrown::{HashMap, HashSet};
use phasst_lib::{
    get_reader, load_assembly_kmers, load_hic, load_hifi, load_linked_read_barcodes, Assembly,
    HicMols, HifiMols, Kmers, LinkedReadBarcodes,
};
use statrs::distribution::{Binomial};
use statrs::distribution::{Univariate};


use clap::App;

use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() {
    println!("Welcome to phasst scaff");
    let params = load_params();
    let kmers = Kmers::load_kmers(&params.het_kmers);
    eprintln!("loading assembly kmers");
    let assembly = load_assembly_kmers(&params.assembly_kmers, &params.assembly_fasta, &kmers);
    let (phasing, kmer_contigs) = load_phased_vcf(&params.phased_vcf, &kmers, &assembly);
    eprintln!("loading hic kmers");
    let hic_mols = load_hic(Some(&params.hic_mols), &kmers);
    eprintln!("building phasing consistency counts");
    let (phasing_consistency_counts, scaffolds) = phasing_consistency(&hic_mols, &phasing, &kmer_contigs, &assembly);
    output_modified_fasta(&scaffolds, &params, &assembly);
    let ordered_scaffolds = order_and_orient(scaffolds, &hic_mols);
}

fn output_modified_fasta(scaffolds: &Scaffold, params: &Params, assembly: &Assembly) {
    eprintln!("{}",params.assembly_fasta);
    let mut reader =  fasta::IndexedReader::from_file(&Path::new(&params.assembly_fasta)).expect("fasta not found");
    let mut scaffold_sizes: HashMap<usize, usize> = HashMap::new();
    for (scaffold, contigs) in scaffolds.chromosomes.iter() {
        let size = scaffold_sizes.entry(*scaffold).or_insert(0);
        for contig in contigs {
            *size += assembly.contig_sizes.get(contig).unwrap();
        }
    }
    let mut size_vec: Vec<(&usize, &usize)> = scaffold_sizes.iter().collect();
    size_vec.sort_by(|a, b| b.1.cmp(a.1));

    let mut writer = fasta::Writer::to_file(Path::new(&format!("{}/chromosomes.fa",params.output))).expect("cannot open fasta writer");
    for (scaffold_id, (scaffold, size)) in size_vec.iter().enumerate() {
        let mut contig_names: Vec<String> = Vec::new();
        for contig in scaffolds.chromosomes.get(scaffold).unwrap().iter() {
            let contig_name = &assembly.contig_names[*contig as usize];
            contig_names.push(contig_name.to_string());
        }
        let full_concat = contig_names.join("_");
        for contig in scaffolds.chromosomes.get(scaffold).unwrap().iter() {
            let contig_name = &assembly.contig_names[*contig as usize];
            reader.fetch_all(contig_name).expect("cannot find contig name in fasta");
            let mut seq: Text = Text::new();
            reader.read(&mut seq).expect("cannot read fasta");
            let id = format!("scaff_{}_contig_{}_{}_{}",scaffold_id+1, contig_name, size, full_concat);
            let record = Record::with_attrs(&id, None, &seq);
            writer.write_record(&record).expect("could not write record");
        }
    }
}

struct OrderedScaffolds {
    scaffolds: Vec<OrderedScaffold>,
}

struct OrderedScaffold {
    order: Vec<i32>,
    orientations: Vec<bool>,
    gaps: Vec<usize>, // gap after each contig, will append 0 on the end to make equal length with other vecs
}

fn order_and_orient(scaffolds: Scaffold, hic_mols: &HicMols) -> OrderedScaffolds {
    let mut scaffolds = OrderedScaffolds { scaffolds: Vec::new() };
    //let length_distribution: = get_empirical_hic_distribution(hic_mols); // index is number of bases apart, value is density
    scaffolds
}

struct LengthDistribution {
    length_distribution: Vec<f32>,
    tail_probability: f32,
}

impl LengthDistribution {
    fn density(&self, length: usize) -> f32 {
        if length/1000 >= self.length_distribution.len() {
            self.tail_probability
        } else {
            self.length_distribution[length/1000]
        }
    }
}

fn get_empirical_hic_distribution(hic_mols: &HicMols) -> LengthDistribution {
    let mut lengths: Vec<usize> = Vec::new();
    let mut to_return = LengthDistribution { length_distribution: Vec::new(), tail_probability: 0.0, };

    for hic in hic_mols.get_hic_molecules() {
        let contig: Option<i32> = None;
    }
    to_return
}

struct PhasingConsistencyCounts {
    counts: HashMap<(i32, i32), PhasingConsistency>,
}

impl PhasingConsistencyCounts {
    fn new() -> PhasingConsistencyCounts {
        PhasingConsistencyCounts {
            counts: HashMap::new(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct PhasingConsistency {
    cis1: u32,
    cis2: u32,
    trans1: u32,
    trans2: u32,
}

impl PhasingConsistency {
    fn new() -> PhasingConsistency {
        PhasingConsistency {
            cis1: 0,
            cis2: 0,
            trans1: 0,
            trans2: 0,
        }
    }
}

struct Scaffold {
    chromosomes: HashMap<usize, Vec<i32>>,
}

fn phasing_consistency(
    hic_mols: &HicMols,
    phasing: &Phasing,
    kmer_contigs: &KmerContigs, assembly: &Assembly,
) -> (PhasingConsistencyCounts,  Scaffold) {
    let mut phasing_consistency_counts = PhasingConsistencyCounts::new();
    let mut components: DisjointSet<i32> = DisjointSet::new();
    for (_name, id) in assembly.contig_ids.iter() {
        components.make_set(*id);
    }
    for hicmol in hic_mols.get_hic_molecules() {
        for vardex1 in 0..hicmol.len() {
            let var1 = hicmol[vardex1];
            if let Some(phase1) = phasing.get_phase(&var1.abs()) {
                if let Some(contigid1) = kmer_contigs.get_contig(&var1.abs()) {
                    for vardex2 in (vardex1 + 1)..hicmol.len() {
                        let var2 = hicmol[vardex2];
                        if let Some(contigid2) = kmer_contigs.get_contig(&var2.abs()){
                            if contigid1 == contigid2 {
                                continue;
                            }
                            if let Some(phase2) = phasing.get_phase(&var2.abs()) {
                                let min = contigid1.min(contigid2);
                                let max = contigid1.max(contigid2);
                                let counts = phasing_consistency_counts
                                    .counts
                                    .entry((*min, *max))
                                    .or_insert(PhasingConsistency::new());
                                if *phase1 && *phase2 {
                                    counts.cis1 += 1;
                                } else if *phase1 && !phase2 {
                                    counts.trans1 += 1;
                                } else if !phase1 && *phase2 {
                                    counts.trans2 += 1;
                                } else {
                                    counts.cis2 += 1;
                                }
                            } 
                        } 
                    }
                } 
            } 
        }
    }
    
    for ((contig1, contig2), counts) in phasing_consistency_counts.counts.iter() {
        let cis = (counts.cis1 + counts.cis2) as f32;
        let trans = (counts.trans1 + counts.trans2) as f32;
        let p_value = binomial_test(cis, trans);
        if counts.cis1 + counts.cis2 + counts.trans1 + counts.trans2 > 150 {
            //if cis > trans && cis/(cis + trans) > 0.9 {
            let min = (counts.cis1 + counts.cis2).min(counts.trans1 + counts.trans2) as f64;
            let max = (counts.cis1 + counts.cis2).max(counts.trans1 + counts.trans2) as f64;
            if p_value < 0.00001 && max/(min+max) > 0.7 { // ? change? keep? test. 
                let min = counts.cis1.min(counts.cis2) as f32;
                if min / cis > 0.25 {
                    components.union(*contig1, *contig2).expect("unable to merge, is this node in the set?");
                    eprintln!("match in cis {} -- {} = {:?}, p-value {}", contig1, contig2, counts, p_value);
                } else {
                    eprintln!("unrelated? {} -- {} = {:?}, p-value {} ", contig1, contig2, counts, p_value);
                }
            } else if trans > cis && trans/(cis + trans) > 0.9 {
                let min = counts.trans1.min(counts.trans2) as f32;
                if min / trans > 0.25 {
                    components.union(*contig1, *contig2).expect("unable to merge, is this node in the set?");
                    eprintln!("match in trans {} -- {} = {:?}, p-value {}", contig1, contig2, counts, p_value);
                } else {
                    eprintln!("unrelated? {} -- {} = {:?}, p-value {}", contig1, contig2, counts, p_value);
                }
            } else {
                eprintln!("unrelated, output for debug {} -- {} = {:?}, p-value {}", contig1, contig2, counts, p_value);
            }
        } else {
            eprintln!("unrelated, output for debug {} -- {} = {:?}, p-value {}", contig1, contig2, counts, p_value);
        }
    }

    let mut sizes: HashMap<usize, usize> = HashMap::new();
    let mut scaffolds: HashMap<usize, Vec<i32>> = HashMap::new();
    for (_name, id) in assembly.contig_ids.iter() {
        let comp = components.find(*id).unwrap_or(0);
        let compvec = scaffolds.entry(comp).or_insert(Vec::new());
        compvec.push(*id);
        let size = assembly.contig_sizes.get(id).unwrap();
        let compsize = sizes.entry(comp).or_insert(0);
        *compsize += *size;
    }
    let mut count_vec: Vec<(&usize, &usize)> = sizes.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    eprintln!("SCAFFOLD SIZES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    for (component, size) in count_vec.iter() {
        eprintln!("scaffold {} size {}", component, size);
    }

    (phasing_consistency_counts, Scaffold { chromosomes: scaffolds} )
}

fn binomial_test(cis: f32, trans: f32) -> f64 {
    let min = cis.min(trans) as f64;
    let max = cis.max(trans) as f64;
    let n = Binomial::new(0.5, (cis+trans) as u64).unwrap();
    let p_value = n.cdf(min) * 2.0;
    p_value
}

struct Phasing {
    phasing: HashMap<i32, bool>,
}

impl Phasing {
    fn get_phase(&self, kmer: &i32) -> Option<&bool> {
        self.phasing.get(kmer)
    }

    fn new() -> Phasing {
        Phasing {
            phasing: HashMap::new(),
        }
    }

    fn set_phase(&mut self, kmer: i32, phase: bool) {
        self.phasing.insert(kmer, phase);
    }
}

struct KmerContigs {
    contigs: HashMap<i32, i32>, // map from kmer id to contig id
}

impl KmerContigs {
    fn get_contig(&self, kmer: &i32) -> Option<&i32> {
        self.contigs.get(kmer)
    }

    fn set_contig(&mut self, kmer: i32, contig: i32) {
        self.contigs.insert(kmer, contig);
    }

    fn new() -> KmerContigs {
        KmerContigs {
            contigs: HashMap::new(),
        }
    }
}

fn load_phased_vcf(vcf: &String, kmers: &Kmers, assembly: &Assembly) -> (Phasing, KmerContigs) {
    let mut to_return: Phasing = Phasing::new();
    let mut kmer_contigs: KmerContigs = KmerContigs::new();
    let f = File::open(vcf.to_string()).expect("Unable to open file");
    let reader = BufReader::new(f);

    for line in reader.lines() {
        let line = line.expect("cannot read phased vcf");
        let toks: Vec<&str> = line.split("\t").collect();
        let reference = toks[3].to_string();
        let alternative = toks[4].to_string();
        let format = toks[9].to_string();
        let gt = format.split(":").collect::<Vec<&str>>()[0].to_string();
        let phasing = gt.split("|").collect::<Vec<&str>>();

        let hap1 = phasing[0].to_string();
        let hap2 = phasing[1].to_string();
        let kmer_id1 = kmers.kmer_ids.get(&reference).unwrap();
        let kmer_id2 = kmers.kmer_ids.get(&alternative).unwrap();
        let contig = toks[0].to_string();
        let contig_id = assembly.contig_ids.get(&contig).expect(&format!("cannot find contig name {}", contig));
        kmer_contigs.set_contig(*kmer_id1, *contig_id);
        kmer_contigs.set_contig(*kmer_id2, *contig_id);
        if hap1 == "0" && hap2 == "1" {
            to_return.set_phase(*kmer_id1, false);
            to_return.set_phase(*kmer_id2, true);
        } else if hap1 == "1" && hap2 == "0" {
            to_return.set_phase(*kmer_id1, true);
            to_return.set_phase(*kmer_id2, false);
        } // else dont add to phasing
    }

    (to_return, kmer_contigs)
}

#[derive(Clone)]
struct Params {
    het_kmers: String,
    linked_read_kmers: Vec<String>,
    hic_mols: Vec<String>,
    output: String,
    assembly_kmers: String,
    assembly_fasta: String,
    phased_vcf: String,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let het_kmers = params.value_of("het_kmers").unwrap();
    let output = params.value_of("output").unwrap();
    let txg_tmp = match params.values_of("txg_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut linked_read_kmers: Vec<String> = Vec::new();
    for x in txg_tmp {
        linked_read_kmers.push(x.to_string());
    }
    let hic_tmp = match params.values_of("hic_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut hic_mols: Vec<String> = Vec::new();
    for x in hic_tmp {
        hic_mols.push(x.to_string());
    }

    let assembly_kmers = params.value_of("assembly_kmers").unwrap();
    let assembly_fasta = params.value_of("assembly_fasta").unwrap();
    let phased_vcf = params.value_of("phased_vcf").unwrap();
    Params {
        het_kmers: het_kmers.to_string(),
        output: output.to_string(),
        linked_read_kmers: linked_read_kmers,
        hic_mols: hic_mols,
        assembly_kmers: assembly_kmers.to_string(),
        assembly_fasta: assembly_fasta.to_string(),
        phased_vcf: phased_vcf.to_string(),
    }
}
