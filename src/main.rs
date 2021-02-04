#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate phasst_lib;
extern crate disjoint_set;

use disjoint_set::DisjointSet;

use hashbrown::{HashMap, HashSet};
use phasst_lib::{
    get_reader, load_assembly_kmers, load_hic, load_hifi, load_linked_read_barcodes, Assembly,
    HicMols, HifiMols, Kmers, LinkedReadBarcodes,
};

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
    let hic_mols = load_hic(&Some(params.hic_mols), &kmers);
    eprintln!("building phasing consistency counts");
    let phasing_consistency_counts = phasing_consistency(&hic_mols, &phasing, &kmer_contigs, &assembly);
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
                                    eprintln!("adding, now {}", counts.cis1);
                                } else if *phase1 && !phase2 {
                                    counts.trans1 += 1;
                                    eprintln!("adding, now {}", counts.trans1);
                                } else if !phase1 && *phase2 {
                                    counts.trans2 += 1;
                                    eprintln!("adding, now {}", counts.trans2);
                                } else {
                                    counts.cis2 += 1;
                                    eprintln!("adding, now {}", counts.cis2);
                                }
                            } else { eprintln!("no phase2"); }
                        } else { eprintln!("no contig2"); }
                    }
                } else { eprintln!("no contig"); }
            } else { eprintln!("no phase"); }
        }
    }
    
    for ((contig1, contig2), counts) in phasing_consistency_counts.counts.iter() {
        if counts.cis1 + counts.cis2 + counts.trans1 + counts.trans2 > 150 {
            let cis = (counts.cis1 + counts.cis2) as f32;
            let trans = (counts.trans1 + counts.trans2) as f32;
            if cis > trans && cis/(cis + trans) > 0.9 {
                let min = counts.cis1.min(counts.cis2) as f32;
                if min / cis > 0.25 {
                    components.union(*contig1, *contig2).expect("unable to merge, is this node in the set?");
                    eprintln!("match in cis {} -- {} = {:?}", contig1, contig2, counts);
                }
            } else if trans > cis && trans/(cis + trans) > 0.9 {
                let min = counts.trans1.min(counts.trans2) as f32;
                if min / trans > 0.25 {
                    components.union(*contig1, *contig2).expect("unable to merge, is this node in the set?");
                    eprintln!("match in trans {} -- {} = {:?}", contig1, contig2, counts);
                }
            } else {
                eprintln!("unrelated, output for debug {} -- {} = {:?}", contig1, contig2, counts);
            }
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

    (phasing_consistency_counts, Scaffold {chromosomes: scaffolds} )
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
        let contig_id = assembly.contig_ids.get(&contig).unwrap();
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
