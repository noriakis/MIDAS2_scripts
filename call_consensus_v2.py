#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## Modified to read MIDAS2 output, and work standalone (2022/12/06)

import argparse, sys, os, gzip, numpy as np, random, csv
import lz4.frame
import io
__version__ = '1.3.0'

class Sample:
	""" Base class for sample """
	def __init__(self, info):
		self.id = info['sample_name']
		self.info = info
		self.mean_depth = float(self.info['mean_coverage'])
		self.fract_cov = float(self.info['fraction_covered'])
		self.consensus = "" # consensus nucleotide sequence

	def filter(self, mean_depth, fract_cov):
		if self.fract_cov < fract_cov:
			return True
		elif self.mean_depth < mean_depth:
			return True
		else:
			return False

class Species:
    """ Base class for species in MIDAS2 """

    def __init__(self, dir):
        if "/" in dir:
            self.spid = dir.split("/")[-1]
        else:
            self.spid = dir
        self.dir = dir
        self.id = os.path.basename(self.dir)
        self.init_paths()
        self.init_files()
        self.init_samples()

    def init_paths(self):
        self.paths = {}
        for type in ['freqs', 'depth', 'info']:
            self.paths[type] = '%s/%s.snps_%s.tsv.lz4' % (self.dir, self.spid, type)
    def init_files(self):
        self.files = {}
        for type in ['freqs', 'depth', 'info']:
            ## Read LZ4 files
            file = lz4.frame.open(self.paths[type]).read().decode('utf-8')
            if type in ['info']:
                self.files[type] = csv.DictReader(io.StringIO(file), delimiter='\t')
            else:
                self.files[type] = csv.reader(io.StringIO(file), delimiter='\t')
        print(self.spid)
        ## All SNV summary is in snps in MIDAS2
        prev_path = os.path.abspath(os.path.join(self.dir, os.pardir))
        file = open(prev_path+'/snps_summary.tsv')
        reader = csv.DictReader(file, delimiter='\t')
        self.files["summary"]=[row for row in reader if self.spid in row["species_id"]]
    def init_samples(self):
        self.sample_ids = None
        for file in ['freqs', 'depth']:
            self.sample_ids = next(self.files[file])[1:]


class GenomicSite:
	""" Base class for genomic sites """
	def __init__(self, species, samples):
		try:
			# fetch site info
			self.info = next(species.files['info'])
			self.id = self.info['site_id']
			# self.ref_allele = self.info['ref_allele']
			self.minor_allele = self.info['minor_allele']
			self.major_allele = self.info['major_allele']
			self.gene_id = self.info['gene_id']
			self.locus_type = self.info['locus_type']
			self.site_type = self.info['site_type']
						
			# copy samples
			self.samples = samples
			
			# fetch site data from freq and depth matrixes
			#	self.samples[sample.id].freq
			#	self.samples[sample.id].depth
			self.fetch_row(species)

		except StopIteration:
			self.id = None
		
	def fetch_row(self, species):
		""" Fetch next row from freq and depth matrices
		    Store in sample objects: sample.freq, sample.depth """
		freqs = next(species.files['freqs'])[1:]
		depths = next(species.files['depth'])[1:]
		for sample in self.samples.values():
			self.samples[sample.id].freq = float(freqs[sample.index])
			self.samples[sample.id].depth = int(depths[sample.index])
				
	def flag_samples(self, site_depth, site_ratio, allele_support):
		""" Filter samples at site based on coverage
			Sets flag: sample.keep = [True/False] 
		
			site_depth: minimum per-sample site depth
			site_ratio: maximum per-sample ratio of site depth to genome depth; used to filter sites with abnormal depth 
			allele_support: minimum fraction of reads supporting consensus allele
		"""
		for sample in self.samples.values():
			sample.flags = []
			sample.keep = True
			if sample.depth < site_depth:
				 sample.keep = False
				 sample.flags.append('site-depth')
			if sample.depth/sample.mean_depth > site_ratio:
				sample.keep = False
				sample.flags.append('depth-ratio')
			if max(sample.freq, 1-sample.freq) < allele_support:
				sample.keep = False
				sample.flags.append('allele-support')

	def filter(self, site_prev=None, site_maf=None, locus_type=None, site_type=None):
		""" determine if site passes quality control """
		self.flags = []
		self.keep = True
		# if self.ref_allele not in ['A','T','C','G']:
		# 	self.flags.append('ref-allele')
		# 	self.keep = False
		if site_prev and self.prevalence < max(1e-6, site_prev):
			self.flags.append('site-prev')
			self.keep = False
		if site_maf and self.pooled_maf < site_maf:
			self.flags.append('site-maf')
			self.keep = False
		if locus_type and self.locus_type != locus_type:
			self.flags.append('locus-type')
			self.keep = False
		if site_type and self.site_type != site_type:
			self.flags.append('site-type')
			self.keep = False

	def compute_prevalence(self):
		self.count_samples = sum([1 for s in self.samples.values() if s.keep])
		return float(self.count_samples)/len(self.samples)

	def summary_stats(self, weight):
		""" compute summary stats for site across samples """
		self.prevalence = self.compute_prevalence()
		self.pooled_maf = self.compute_pooled_maf(weight)

	def compute_pooled_maf(self, weight=False):
		""" compute average frequency of reference allele with optional weighting of samples """
		if self.count_samples == 0:
			return 0.0
		elif weight:
			depth = sum([s.depth for s in self.samples.values() if s.keep])
			maf = sum([s.depth * s.freq for s in self.samples.values() if s.keep])
			return maf/depth
		else:
			return np.mean([s.freq for s in self.samples.values() if s.keep])

	def resample_reads(self, rand_reads, replace_reads):
		""" resample random number of reads per sample """
		index = 0
		for sample in self.samples.values():
			sample.depth = rand_reads
			if sample.freq > 0 and sample.freq < 1:
				count_minor = int(round(sample.freq * sample.depth))
				count_major = sample.depth - count_minor
				alleles = np.random.choice([1]*count_minor+[0]*count_major, rand_reads, replace=replace_reads)
				sample.freq = np.mean(alleles)
			index += 1
			
	def call_consensus(self):
		""" call consensus allele at each position """
		for sample in self.samples.values():
			sample.freq = round(sample.freq)
			
	def fetch_consensus(self, sample):
		if not sample.keep:
			return '-'
		elif sample.depth == 0:
			return '-'
		elif sample.freq >= 0.5:
			return self.minor_allele
		else:
			return self.major_allele
			
def fetch_samples(species, mean_depth=0, fract_cov=0, max_samples=float('inf'),
                  keep_samples=None, exclude_samples=None, rand_samples=None):
	""" Select samples from input
		mean_depth: filter samples based on average genome/gene depth
		fract_cov: filter samples based on the number of genomic sites covered by >=1 read
		max_samples: stop when this numbe of samples reached
		keep_samples: only keep samples in this list
		exclude_samples: exclude any sample in this list
		rand_samples: select random subset of samples
	"""
	samples = {}
	# select samples that pass filters
	for index, info in enumerate(species.files['summary']):
		sample = Sample(info)
		sample.index = index
		if sample.filter(mean_depth, fract_cov):
			continue
		if keep_samples and sample.id not in keep_samples:
			continue
		if exclude_samples and sample.id in exclude_samples:
			continue
		if len(samples) >= max_samples:
			continue
		samples[sample.id] = sample
	# check there's at least 1 sample
	if len(samples) == 0:
		error = "\nError: no samples satisfied your selection criteria.\n"
		error += "Try running again with more lenient parameters\n"
		sys.exit(error)
	# select random subset of samples
	if rand_samples:
		if rand_samples > len(samples):
			error = "\nError: --rand_samples cannot exceed the number of samples\n"
			sys.exit(error)
		ids = np.random.choice(samples.keys(), rand_samples, replace=False)
		for id in samples.keys():
			if id not in ids:
				del samples[id]
	return samples

def fetch_sites(species, samples):
	""" yield genomic sites from species across samples """
	index = 0
	while True:
		site = GenomicSite(species, samples)
		if not site.id:
			break
		else:
			index += 1
			yield site

def print_copyright(log=None):
	lines = []
	lines.append("")
	lines.append("MIDAS: Metagenomic Intra-species Diversity Analysis System")
	lines.append("version %s; github.com/snayfach/MIDAS" % __version__)
	lines.append("Copyright (C) 2015-2016 Stephen Nayfach")
	lines.append("Freely distributed under the GNU General Public License (GPLv3)")
	lines.append("")
	if log is not None: log.write('\n'.join(lines)+'\n')
	sys.stdout.write('\n'.join(lines)+'\n')

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description:
Build FASTA file of consensus sequences for a species per sample
Useful for building phylogenetic trees
Before running this script, you'll need to have run `merge_midas.py snps`

Usage: call_consensus_v2.py indir [options]
""",
		epilog="""
Examples:
1) Build multi-FASTA of core-genome sequences (recommended)
-core-genome sites defined as >=5 reads in >=90% of samples
-use only variable positions (>=1% minor allele frequency across samples)
-only include samples with sufficient data (>=10x mean-depth, >=40% of sites with >=1 mapped read)
-exclude sites with abnormal depth (>5x mean-depth or <1/5 mean-depth)

call_consensus_v2.py /path/to/snps --out /path/to/seqs --site_maf 0.01 --site_depth 5 --site_prev 0.90 --sample_depth 10 --sample_cov 0.40 --site_ratio 5.0

2) Build multi-FASTA using defaults
call_consensus_v2.py /path/to/snps --out /path/to/seqs

3) Run a quick test
call_consensus_v2.py /path/to/snps --out /path/to/output --max_sites 10000

""")
	parser.add_argument('indir', metavar='PATH', type=str,
		help="""path to output from `merge_midas.py snps` for one species
directory should be named according to a species_id and contains files 'snps_*.txt')""")
	parser.add_argument('--out', metavar='PATH', type=str, default="/dev/stdout",
		help="""path to output file""")

	sample = parser.add_argument_group("Sample filters (select subset of samples from INDIR)")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=0.0, metavar='FLOAT',
		help="""minimum average read depth per sample (0.0)""")
	sample.add_argument('--sample_cov', dest='fract_cov', type=float, default=0.0, metavar='FLOAT',
		help="""fraction of reference sites covered by at least 1 read (0.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT', default=float('Inf'),
		help="""maximum number of samples to process.
useful for quick tests (use all)""")
	sample.add_argument('--keep_samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to include
samples will still be subject to other filters""")
	sample.add_argument('--exclude_samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to exclude.
samples will still be subject to other filters""")

	snps = parser.add_argument_group("Site filters (select subset of genomic sites from INDIR)")
	snps.add_argument('--site_list', metavar='PATH',type=str,
		help="""path to list of sites to include; other filters still apply""")
	snps.add_argument('--site_depth', type=int, default=2, metavar='INT',
		help="""minimum number of mapped reads per site (2)""")
	snps.add_argument('--site_prev', type=float, default=0.0, metavar='FLOAT',
		help="""site has at least <site_depth> coverage in at least <site_prev> proportion of samples (0.0)
a value of 1.0 will select sites that have sufficent coverage in all samples.
a value of 0.0 will select all sites, including those with low coverage in many samples 
NAs recorded for included sites with less than <site_depth> in a sample """)
	snps.add_argument('--site_maf', type=float, default=0.0, metavar='FLOAT',
		help="""minimum average-minor-allele-frequency of site across samples (0.0)
setting this above zero (e.g. 0.01, 0.02, 0.05) will only retain variable sites
by default invariant sites are also retained.""")
	snps.add_argument('--site_ratio', type=float, default=float('Inf'), metavar='FLOAT',
		help="""maximum ratio of site-depth to mean-genome-depth (None)
a value of 10 will filter genomic sites with 10x greater coverage than the genomic background""")
	snps.add_argument('--allele_support', type=float, default=0.5, metavar='FLOAT',
		help="minimum fraction of reads supporting consensus allele (0.5)")	
	snps.add_argument('--locus_type', choices=['CDS', 'RNA', 'IGR'],
		help="""use genomic sites that intersect: 'CDS': coding genes, 'RNA': rRNA and tRNA genes, 'IGS': intergenic regions""")
	snps.add_argument('--site_type', choices=['1D','2D','3D','4D'],
		help="""if locus_type == 'CDS', use genomic sites with specified degeneracy: 4D indicates synonymous and 1D non-synonymous sites""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""maximum number of sites to include in output (use all)
useful for quick tests""")
	args = vars(parser.parse_args())
	keep = args['keep_samples'].rstrip(',').split(',') if args['keep_samples'] else None
	exclude = args['exclude_samples'].rstrip(',').split(',') if args['exclude_samples'] else None
	return args

def print_args(args):
	lines = []
	
def print_args(args):
	lines = []
	lines.append("Command: %s" % ' '.join(sys.argv))
	lines.append("Script: call_consensus_v2.py")
	lines.append("Input directory: %s" % args['indir'])
	lines.append("Output file: %s" % args['out'])
	lines.append("Sample filters:")
	lines.append("  sample_depth: %s" % args['sample_depth'])
	lines.append("  fract_cov: %s" % args['fract_cov'])
	lines.append("  max_samples: %s" % args['max_samples'])
	lines.append("  keep_samples: %s" % args['keep_samples'])
	lines.append("  exclude_samples: %s" % args['exclude_samples'])
	lines.append("Site filters:")
	lines.append("  site_list: %s" % args['site_list'])
	lines.append("  site_depth: %s" % args['site_depth'])
	lines.append("  site_prev: %s" % args['site_prev'])
	lines.append("  site_maf: %s" % args['site_maf'])
	lines.append("  site_ratio: %s" % args['site_ratio'])
	lines.append("  allele_support: %s" % args['allele_support'])
	lines.append("  locus_type: %s" % args['locus_type'])
	lines.append("  site_type: %s" % args['locus_type'])	
	lines.append("  max_sites: %s" % args['max_sites'])
	sys.stdout.write('\n'.join(lines)+'\n')

def check_args(args):
	if not os.path.isdir(args['indir']):
		sys.exit("\nError: Specified input directory '%s' does not exist\n" % args['indir'])
	if args['site_depth'] < 1:
		sys.exit("\nError: --site_depth must be >=1\n")
	if args['max_sites'] < 1:
		sys.exit("\nError: --max_sites must be >= 1 to calculate nucleotide variation\n")
	if args['max_samples'] < 1:
		sys.exit("\nError: --max_samples must be >= 1 to calculate nucleotide variation\n")
	if args['site_ratio'] < 0:
		sys.exit("\nError: --site_ratio cannot be a negative number\n")
	if args['site_depth'] < 0:
		sys.exit("\nError: --site_depth cannot be a negative number\n")
	if args['sample_depth'] < 0:
		sys.exit("\nError: --sample_depth cannot be a negative number\n")
	if not 0 <= args['site_maf'] <= 1:
		sys.exit("\nError: --site_maf must be between 0 and 1\n")
	if not 0 <= args['site_prev'] <= 1:
		sys.exit("\nError: --site_prev must be between 0 and 1\n")
	if not 0 <= args['fract_cov'] <= 1:
		sys.exit("\nError: --fract_cov must be between 0 and 1\n")

def percent_missing(seq):
	if len(seq) > 0:
		return round(100*seq.count('-')/float(len(seq)),2)
	else:
		return 'NA'

def sequence_description(sample):
	desc = {}
	desc['length'] = len(sample.consensus)
	desc['percent_missing'] = percent_missing(sample.consensus)
	desc['mean_depth'] = round(sample.mean_depth, 2)
	return desc

def write_consensus(args, samples):
	outfile=open(args['out'], 'w')
	for sample_id in sorted(samples):
		sample = samples[sample_id]
		desc = sequence_description(sample)
		outfile.write('>'+sample.id+'\t')
		outfile.write(' '.join(['%s=%s' % (key, value) for key, value in desc.items()])+'\n')
		outfile.write(sample.consensus+'\n')
	outfile.close()

def format_site_type(site_type):
	if site_type == 'ALL':
		return ['NC','1D','2D','3D','4D']
	elif site_type == 'CDS':
		return ['1D','2D','3D','4D']
	else:
		return [args['site_type']]
		
if __name__ == '__main__':
	
	# setup args
	args = parse_arguments()
	check_args(args)
	print_copyright()
	print_args(args)
	
	# init species, samples, sequences, site list
	species = Species(args['indir'])
	samples = fetch_samples(species, args['sample_depth'], args['fract_cov'], args['max_samples'])
	if args['site_list']: 
		site_list = set([_.rstrip() for _ in open(args['site_list'])])
	
	# loop over genomic sites
	sites = fetch_sites(species, samples)
	retained_sites = 0
	for index, site in enumerate(sites):
	
		# stop early
		if retained_sites >= args['max_sites']: break
			
		# prune low quality samples for site:
		#   site.samples['sample'].keep = [True|False]
		site.flag_samples(args['site_depth'], args['site_ratio'], args['allele_support'])
		
		# compute site summary stats
		#   site.prevalence
		#   site.pooled_maf
		site.summary_stats(weight=False)
		
		# filter genomic site
		#   site.keep = [True|False]
		if not args['site_list']:
			site.filter(args['site_prev'], args['site_maf'], args['locus_type'], args['site_type']) 
		elif site.id in site_list:
			site.keep = True
		else:
			site.keep = False
		
		# store consensus
		if site.keep:
			retained_sites += 1		
			for sample in site.samples.values():
				samples[sample.id].consensus += site.fetch_consensus(sample)

	# write consensus
	write_consensus(args, samples)

