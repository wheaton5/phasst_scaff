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
    let hic_mols = load_hic(Some(&params.hic_mols), &kmers, false);
    eprintln!("building phasing consistency counts");
    let (phasing_consistency_counts, scaffolds) = phasing_consistency(&hic_mols, &phasing, &kmer_contigs, &assembly);
    let biggest_to_smallest_scaffolds = output_modified_fasta(&scaffolds,  &params, &assembly);
    //let ordered_scaffolds = order_and_orient(biggest_to_smallest_scaffolds, order_and_orientation_counts, &hic_mols, &assembly);
}

fn output_modified_fasta(scaffolds: &Scaffold,  params: &Params, assembly: &Assembly) -> Vec<Vec<i32>> {
    eprintln!("{}",params.assembly_fasta);
    let reader =  fasta::Reader::from_file(Path::new(&params.assembly_fasta)).expect("fasta not found");
    let path = Path::new(&params.assembly_fasta);
    let mut reader =  fasta::IndexedReader::from_file(&path).expect(&format!("fasta not found {}", path.to_str().unwrap()));
    let mut scaffold_sizes: HashMap<usize, usize> = HashMap::new();
    for (scaffold, contigs) in scaffolds.chromosomes.iter() {
        let size = scaffold_sizes.entry(*scaffold).or_insert(0);
        for contig in contigs {
            *size += assembly.contig_sizes.get(contig).unwrap();
        }
    }
    let mut size_vec: Vec<(&usize, &usize)> = scaffold_sizes.iter().collect();
    size_vec.sort_by(|a, b| b.1.cmp(a.1));
    let mut biggest_to_smallest_scaffolds: Vec<Vec<i32>> = Vec::new();
    for (scaffold, _size) in size_vec.iter() {
        let mut contigs: Vec<i32> = Vec::new();
        for contig in scaffolds.chromosomes.get(scaffold).unwrap().iter() {
            contigs.push(*contig);
        }
        biggest_to_smallest_scaffolds.push(contigs);
    }


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
    biggest_to_smallest_scaffolds
}

struct OrderedScaffolds {
    scaffolds: Vec<OrderedScaffold>,
}

struct OrderedScaffold {
    order: Vec<i32>,
    orientations: Vec<bool>,
    gaps: Vec<usize>, // gap after each contig, will append 0 on the end to make equal length with other vecs
}

fn order_and_orient(biggest_to_smallest_scaffolds: Vec<Vec<i32>>, order_orientation_counts: HashMap<(i32, i32), OrderOrient>, hic_mols: &HicMols, assembly: &Assembly) -> OrderedScaffolds {
    let mut scaffolds = OrderedScaffolds { scaffolds: Vec::new() };
    let length_distribution = get_empirical_hic_distribution(hic_mols, assembly); // index is number of bases apart, value is density
    for (scaffdex, scaffold) in biggest_to_smallest_scaffolds.iter().enumerate() {
        let mut order_and_orientation_counts: HashMap<(i32, i32), [u32; 4]> = HashMap::new();
        let mut order_and_orientation_totals: HashMap<(i32, i32), u32> = HashMap::new();
        eprintln!("Order and Orientation counts for scaffold {}", scaffdex);
        for contig1dex in 0..scaffold.len() {
            let contig1 = scaffold[contig1dex];
            for contig2dex in (contig1dex+1)..scaffold.len() {
                let contig2 = scaffold[contig2dex];
                let min = contig1.min(contig2);
                let max = contig1.max(contig2);
                let counts = order_and_orientation_counts.get(&(min, max)).unwrap();
                eprintln!("\tcontigs {} -- {} have counts {:?}", min, max, counts);
            }
        }
    }
    scaffolds
}

struct LengthDistribution {
    length_distribution: Vec<f32>,
    tail_probability: f32,
}

impl LengthDistribution {
    fn new() -> LengthDistribution {
        let mut length_distribution: Vec<f32> = Vec::new();
        for _index in 0..10000 {
            length_distribution.push(0.0);
        }
        LengthDistribution { length_distribution: length_distribution, tail_probability: 0.0, }
    }

    fn density(&self, length: usize) -> f32 {
        if length/1000 >= self.length_distribution.len() {
            self.tail_probability
        } else {
            self.length_distribution[length/1000]
        }
    }
}

fn get_empirical_hic_distribution(hic_mols: &HicMols, assembly: &Assembly) -> LengthDistribution {
    let mut to_return = LengthDistribution::new();
    let mut total_count = 0.0;
    for hic in hic_mols.get_hic_molecules() {
        let mut contig: Option<i32> = None;
        let mut min = usize::MAX;
        let mut max = 0;
        for var in hic {
            if let Some((contig_id, position)) = kmer_contig_position(*var, assembly) {
                if contig == None {
                    contig = Some(contig_id);
                    min = min.min(position);
                    max = max.max(position);
                } else if contig == Some(contig_id) {
                    min = min.min(position);
                    max = max.max(position);
                }
            }
        }
        if min != usize::MAX && max - min > 0 {
            if (max-min)/1000 < to_return.length_distribution.len() {
                to_return.length_distribution[(max-min)/1000] += 1.0;
                total_count += 1.0;
            }
        }
    }
    for index in 0..to_return.length_distribution.len() {
        to_return.length_distribution[index] /= total_count;
    }
    to_return
}

fn kmer_contig_position(kmer: i32, assembly: &Assembly) -> Option<(i32, usize)> {
    if let Some((contig_id, number_seen, _order, position)) = assembly.variants.get(&kmer.abs()) {
        if *number_seen == 1 {
            return Some((*contig_id, *position));
        }
    } else if let Some((contig_id, number_seen, _order, position)) = assembly.variants.get(&Kmers::pair(kmer.abs())) {
        if *number_seen == 1 {
            return Some((*contig_id, *position));
        }
    }
    None
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

struct OrderOrientPhase {
    cis: OrderOrient,
    trans: OrderOrient,
}

#[derive(Debug, Clone, Copy)]
struct OrderOrient {
    start1_start2: u32,
    start1_end2: u32,
    end1_start2: u32,
    end1_end2: u32,
}

impl OrderOrient {
    fn new() -> OrderOrient {
        OrderOrient{start1_start2: 0, start1_end2: 0, end1_start2: 0, end1_end2: 0}
    }
}
impl OrderOrientPhase {
    fn new() -> OrderOrientPhase {
        OrderOrientPhase{ cis: OrderOrient::new(), trans: OrderOrient::new() }
    }
}

fn phasing_consistency(
    hic_mols: &HicMols,
    phasing: &Phasing,
    kmer_contigs: &KmerContigs, assembly: &Assembly,
) -> (PhasingConsistencyCounts,  Scaffold) {
    let mut phasing_consistency_counts = PhasingConsistencyCounts::new();
    let mut components: DisjointSet<i32> = DisjointSet::new();
    let mut overly_complex_data_structure: HashMap<(i32, i32), OrderOrientPhase> = HashMap::new();
    //let mut order_and_oriention_counts: HashMap<(i32, i32), OrderOrient> = HashMap::new();

    let mut num_assembly_kmers_start_end: HashMap<(i32, bool),f32> = HashMap::new();
    for (_kmer, (contig_id, number_seen, _order, position)) in assembly.variants.iter() {
        if *number_seen == 1 {
            let contig_size = assembly.contig_sizes.get(contig_id).unwrap();
            if *position < 100000 {
               
                let count = num_assembly_kmers_start_end.entry((*contig_id, true)).or_insert(0.0);
                *count += 1.0;
                 eprintln!("adding count {}",count);
            } else if contig_size - position < 100000 {
                let count = num_assembly_kmers_start_end.entry((*contig_id, false)).or_insert(0.0);
                *count += 1.0;
                 eprintln!("adding count {}",count);
            }
        }
    }

    for (_name, id) in assembly.contig_ids.iter() {
        components.make_set(*id);
    }
    for hicmol in hic_mols.get_hic_molecules() {
        for vardex1 in 0..hicmol.len() {
            let var1 = hicmol[vardex1];
            if let Some(phase1) = phasing.get_phase(&var1.abs()) {
                let mut phase1 = phase1;
                if let Some(contigid1) = kmer_contigs.get_contig(&var1.abs()) {
                    //let (_contig_id, _number_seen, _order, position1) = assembly.variants.get(&var1.abs()).unwrap();
                    let (_contig_id, position1) = kmer_contig_position(var1.abs(), assembly).expect("why");
                    let contig1_size = assembly.contig_sizes.get(contigid1).unwrap();
                    let mut contig1_start_or_end = None;
                    if *contig1_size - position1 < 100000 {
                        contig1_start_or_end = Some(false);
                    } else if position1 < 100000 {
                        contig1_start_or_end = Some(true);
                    }
                    for vardex2 in (vardex1 + 1)..hicmol.len() {
                        let var2 = hicmol[vardex2];
                        if let Some(contigid2) = kmer_contigs.get_contig(&var2.abs()){
                            if contigid1 == contigid2 {
                                continue;
                            }
                            if let Some(phase2) = phasing.get_phase(&var2.abs()) {
                                let mut phase2 = phase2;
                                //let (_contig_id, _number_seen, _order, position2) = assembly.variants.get(&var2.abs()).unwrap();
                                let (_contig_id, position2) = kmer_contig_position(var2.abs(), assembly).expect("why");
                                let contig2_size = assembly.contig_sizes.get(contigid2).unwrap();
                                let mut contig2_start_or_end = None;
                                //if (position2 as f32)/(*contig2_size as f32) > 0.5 {
                                if *contig2_size - position2 < 100000 {
                                    contig2_start_or_end = Some(false);
                                } else if position2 < 100000 {
                                    contig2_start_or_end = Some(true);
                                }
                                if contigid1 > contigid2 {
                                    let tmp = contig1_start_or_end;
                                    contig1_start_or_end = contig2_start_or_end;
                                    contig2_start_or_end = tmp;
                                    let tmp = phase1;
                                    phase1 = phase2;
                                    phase2 = tmp;
                                }
                                let min = contigid1.min(contigid2);
                                let max = contigid1.max(contigid2);
                                let counts = phasing_consistency_counts
                                    .counts
                                    .entry((*min, *max))
                                    .or_insert(PhasingConsistency::new());
                                let mut cis_or_trans = true;
                                if *phase1 && *phase2 {
                                    counts.cis1 += 1;
                                } else if *phase1 && !phase2 {
                                    counts.trans1 += 1;
                                    cis_or_trans = false;
                                } else if !phase1 && *phase2 {
                                    counts.trans2 += 1;
                                    cis_or_trans = false;
                                } else {
                                    counts.cis2 += 1;
                                }
                                let order_orient_phase = overly_complex_data_structure.entry((*min,*max)).or_insert(OrderOrientPhase::new());
                                if let Some(contig1_se) = contig1_start_or_end {
                                    if let Some(contig2_se) = contig2_start_or_end {
                                        if cis_or_trans {
                                            if contig1_se && contig2_se {
                                                order_orient_phase.cis.start1_start2 += 1;
                                            } else if contig1_se && !contig2_se {
                                                order_orient_phase.cis.start1_end2 += 1;
                                            } else if !contig1_se && contig2_se {
                                                order_orient_phase.cis.end1_start2 += 1;
                                            } else {
                                                order_orient_phase.cis.end1_end2 += 1;
                                            }   
                                        } else {
                                            if contig1_se && contig2_se {
                                                order_orient_phase.trans.start1_start2 += 1;
                                            } else if contig1_se && !contig2_se {
                                                order_orient_phase.trans.start1_end2 += 1;
                                            } else if !contig1_se && contig2_se {
                                                order_orient_phase.trans.end1_start2 += 1;
                                            } else {
                                                order_orient_phase.trans.end1_end2 += 1;
                                            }
                                        } // end cis or trans
                                    } // end let Some(contig1_se)
                                } // end let Some(contig2_se)
                            } // end let Some(phase2)
                        } // end contig2_id 
                    } // end var2dex
                } // end contig1_id
            } // end phase1
        } // end var1dex
    } // end hic read loop
    
    let mut kmer_coverages: Vec<f64> = Vec::new();
    for ((contig1, contig2), counts) in phasing_consistency_counts.counts.iter() {
        let cis = (counts.cis1 + counts.cis2) as f32;
        let trans = (counts.trans1 + counts.trans2) as f32;
        //let min = (cis).min(trans) as f64;
        let max = (cis).max(trans) as f64;
        let contig1_kmers = assembly.molecules.get(contig1).expect(&format!("couldnt load contig {}", assembly.contig_names[*contig1 as usize])).len() as f64;
        let contig2_kmers = assembly.molecules.get(contig2).expect(&format!("couldnt load contig {}", assembly.contig_names[*contig2 as usize])).len() as f64;
        let dominant_kmers = contig1_kmers.min(contig2_kmers);
        let coverage = max/dominant_kmers;
        kmer_coverages.push(coverage);
    }
    kmer_coverages.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let coverage_threshold = 10.0 * kmer_coverages[(0.5 * kmer_coverages.len() as f64) as usize]; // threshold is 10 * median

    for ((contig1, contig2), counts) in phasing_consistency_counts.counts.iter() {
        let cis = (counts.cis1 + counts.cis2) as f32;
        let trans = (counts.trans1 + counts.trans2) as f32;
        let p_value = binomial_test(cis, trans);
        let contig1_kmers = assembly.molecules.get(contig1).unwrap().len() as f64;
        let contig2_kmers = assembly.molecules.get(contig2).unwrap().len() as f64;
        let dominant_kmers = contig1_kmers.min(contig2_kmers);

        let min = (counts.cis1 + counts.cis2).min(counts.trans1 + counts.trans2) as f64;
        let max = (counts.cis1 + counts.cis2).max(counts.trans1 + counts.trans2) as f64;
        let coverage = max/dominant_kmers;

        if counts.cis1 + counts.cis2 + counts.trans1 + counts.trans2 > 100 {
            if cis > trans && p_value < 0.000001 && max/(min+max) > 0.8 { // ? change? keep? test. 
                let min = counts.cis1.min(counts.cis2) as f32;
                if min / cis > 0.25 && coverage > coverage_threshold{
                    components.union(*contig1, *contig2).expect("unable to merge, is this node in the set?");
                    eprintln!("match in cis {} -- {} = {:?}, kmer coverage {}, p-value {}", assembly.contig_names[*contig1 as usize], 
                        assembly.contig_names[*contig2 as usize], counts, coverage, p_value);
                    //let counts = order_and_oriention_counts.entry((*contig1, *contig2)).or_insert(OrderOrient::new());
                    let mut counts: [f32;4] = [0.0;4];
                    let getcounts = overly_complex_data_structure.get(&(*contig1, *contig2)).unwrap();
                    //let start1 = *num_assembly_kmers_start_end.get(&(*contig1, true)).unwrap();
                    //let end1 = *num_assembly_kmers_start_end.get(&(*contig1, false)).unwrap();
                    //let start2 = *num_assembly_kmers_start_end.get(&(*contig2, true)).unwrap();
                    //let end2 = *num_assembly_kmers_start_end.get(&(*contig2, false)).unwrap();
                    counts[0] = getcounts.cis.start1_start2 as f32;//(start1.min(start2));
                    counts[1] = getcounts.cis.start1_end2 as f32;//(start1.min(end2));
                    counts[2] = getcounts.cis.end1_start2 as f32;//(end1.min(start2));
                    counts[3] = getcounts.cis.end1_end2 as f32;//(end1.min(end2));
                    eprintln!("\torder and orientation counts {:?}", counts);
                } else {
                    eprintln!("unrelated . . {} -- {} = {:?}, kmer coverage {}, p-value {} ", 
                        assembly.contig_names[*contig1 as usize], assembly.contig_names[*contig2 as usize], 
                        counts, coverage, p_value);
                }
            } else if p_value < 0.000001 && max/(min+max) > 0.8  {
                let min = counts.trans1.min(counts.trans2) as f32;
                if min / trans > 0.25  && coverage > coverage_threshold {
                    components.union(*contig1, *contig2).expect("unable to merge, is this node in the set?");
                    eprintln!("match in trans {} -- {} = {:?}, kmer coverage {}, p-value {}", 
                        assembly.contig_names[*contig1 as usize], assembly.contig_names[*contig2 as usize], counts, coverage, p_value);
                    //let counts = order_and_oriention_counts.entry((*contig1, *contig2)).or_insert(OrderOrient::new());
                    let mut counts: [f32;4] = [0.0;4];
                    let getcounts = overly_complex_data_structure.get(&(*contig1, *contig2)).unwrap();
                    //let start1 = *num_assembly_kmers_start_end.get(&(*contig1, true)).unwrap();
                    //let end1 = *num_assembly_kmers_start_end.get(&(*contig1, false)).unwrap();
                    //let start2 = *num_assembly_kmers_start_end.get(&(*contig2, true)).unwrap();
                    //let end2 = *num_assembly_kmers_start_end.get(&(*contig2, false)).unwrap();
                    counts[0] = getcounts.trans.start1_start2 as f32;//(start1.min(start2));

                    counts[1] = getcounts.trans.start1_end2 as f32;//(start1.min(end2));
                    counts[2] = getcounts.trans.end1_start2 as f32;//(end1.min(start2));
                    counts[3] = getcounts.trans.end1_end2 as f32;//(end1.min(end2));
                    eprintln!("\torder and orientation counts {:?}", counts);
                } else {
                    eprintln!("unrelated . . {} -- {} = {:?}, kmer coverage {}, p-value {}", assembly.contig_names[*contig1 as usize], 
                    assembly.contig_names[*contig2 as usize], counts, coverage, p_value);
                }
            } else {
                eprintln!("unrelated . . {} -- {} = {:?}, kmer coverage {}, p-value {}", 
                assembly.contig_names[*contig1 as usize], assembly.contig_names[*contig2 as usize], counts, coverage, p_value);
            }
        } else {
            eprintln!("unrelated . . {} -- {} = {:?}, kmer coverage {}, p-value {}", 
            assembly.contig_names[*contig1 as usize], assembly.contig_names[*contig2 as usize], counts, coverage, p_value);
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
        if phasing.len() == 1 {
            continue;
        }
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
    let mut txg_mols: Vec<String> = Vec::new();
    match params.value_of("linked_read_barcodes") {
        Some(txg_fofn) => {
            let f = File::open(txg_fofn).expect("Unable to open txg fofn");
            let f = BufReader::new(f);

            for line in f.lines() {
                let line = line.expect("Unable to read txg fofn line");
                txg_mols.push(line.to_string());
            }
        },
        None => (),
    }


    let mut hic_mols: Vec<String> = Vec::new();
    match params.value_of("hic_mols") {
        Some(hic_fofn) => {
            let f = File::open(hic_fofn).expect("Unable to open hic fofn");
            let f = BufReader::new(f);

            for line in f.lines() {
                let line = line.expect("Unable to read txg fofn line");
                hic_mols.push(line.to_string());
            }
        },
        None => (),
    }

    let assembly_kmers = params.value_of("assembly_kmers").unwrap();
    let assembly_fasta = params.value_of("assembly_fasta").unwrap();
    let phased_vcf = params.value_of("phased_vcf").unwrap();
    Params {
        het_kmers: het_kmers.to_string(),
        output: output.to_string(),
        linked_read_kmers: txg_mols,
        hic_mols: hic_mols,
        assembly_kmers: assembly_kmers.to_string(),
        assembly_fasta: assembly_fasta.to_string(),
        phased_vcf: phased_vcf.to_string(),
    }
}
