import vcf
import sys
from collections import Counter, defaultdict

def analyze_vcf(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        # Skip non-biallelic sites
        if len(record.ALT) != 1:
            continue  # skip multiallelic sites

        genotypes = {}
        allele_to_samples = defaultdict(set)
        valid = True

        for sample in record.samples:
            sample_name = sample.sample.replace(".DmagnaLRV01.filtered.sorted.nd.bam", "")
            gt = sample['GT']
            if gt is None or '.' in gt:
                valid = False
                break

            alleles = gt.replace('|', '/').split('/')

            # Skip if alleles other than 0 or 1 are present
            if any(a not in {'0', '1'} for a in alleles):
                valid = False
                break

            genotypes[sample_name] = alleles
            for allele in set(alleles):
                allele_to_samples[allele].add(sample_name)

        if not valid:
            continue

        # Check if exactly one sample carries a unique allele
        unique_allele_sample = None
        if len(allele_to_samples['1']) == 1 and len(allele_to_samples['0']) > 1:
            unique_allele = '1'
            unique_allele_sample = list(allele_to_samples['1'])[0]
            most_frequent_allele = '0'
        elif len(allele_to_samples['0']) == 1 and len(allele_to_samples['1']) > 1:
            unique_allele = '0'
            unique_allele_sample = list(allele_to_samples['0'])[0]
            most_frequent_allele = '1'
        else:
            continue  # No uniquely carried allele

        # Get allele bases
        alleles = record.alleles  # [REF, ALT]
        sample_alleles = genotypes[unique_allele_sample]
        genotype_nucleotides = '/'.join([str(alleles[int(a)]) for a in sample_alleles])
        most_freq_nt = str(alleles[int(most_frequent_allele)])
        least_freq_nt = str(alleles[int(unique_allele)])

        print(f"{unique_allele_sample}\t{record.CHROM}\t{record.POS}\t{genotype_nucleotides}\t{most_freq_nt}\t{least_freq_nt}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_vcf.py <input.vcf>")
        sys.exit(1)
    analyze_vcf(sys.argv[1])
