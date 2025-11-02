import csv
import gzip
import sys

def parse_sample_file(sample_file):
    sample_positions = []
    with open(sample_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            sample = row[0]
            chromosome = row[1]
            position = int(row[2])
            ref_alt = row[3]
            most_freq_nt = row[4]
            least_freq_nt = row[5]
            sample_positions.append((sample, chromosome, position, ref_alt, most_freq_nt, least_freq_nt))
    return sample_positions

def parse_vcf(vcf_file, sample_positions):
    position_dict = {
        (chrom, pos): (sample, ref_alt, most_nt, least_nt)
        for sample, chrom, pos, ref_alt, most_nt, least_nt in sample_positions
    }

    positions_to_find = set(position_dict.keys())
    found_positions = set()
    header_printed = False

    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    sample_name = header[9]  # assumes single-sample VCF
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])

            key = (chrom, pos)
            if key not in positions_to_find:
                continue

            found_positions.add(key)
            ref = fields[3]
            alt_alleles = [a for a in fields[4].split(',') if a != '<*>']
            all_alleles = [ref] + alt_alleles  # index 0 = ref, index 1+ = alts

            format_fields = fields[8].split(':')
            sample_data = fields[9].split(':')
            ad_index = format_fields.index('AD') if 'AD' in format_fields else None

            sample, ref_alt, most_nt, least_nt = position_dict[key]
            most_freq = least_freq = "NA"

            if ad_index is not None:
                ad_field = sample_data[ad_index]
                try:
                    ad_values = list(map(int, ad_field.split(',')))
                    total = sum(ad_values)

                    if total > 0:
                        # Frequency of least frequent nucleotide (unique allele)
                        if least_nt in all_alleles:
                            j = all_alleles.index(least_nt)
                            least_freq = ad_values[j] / total if j < len(ad_values) else 0.0
                        else:
                            least_freq = 0.0

                        # Frequency of most frequent nucleotide (common allele)
                        if most_nt in all_alleles:
                            i = all_alleles.index(most_nt)
                            most_freq = ad_values[i] / total if i < len(ad_values) else 0.0
                        else:
                            most_freq = 0.0
                except (ValueError, IndexError):
                    most_freq = least_freq = "NA"

            alleles_in_ref_alt = ref_alt.split('/')
            hom_or_het = 1 if len(alleles_in_ref_alt) == 2 and alleles_in_ref_alt[0] == alleles_in_ref_alt[1] else 2
            vcf_genotype = f"{most_nt}/{least_nt}"

            if not header_printed:
                print("Sample\tChromosome\tPosition\tVCF_Genotype\tPositions_Genotype\tRef_Allele_Frequency\tAlt_Allele_Frequency\tHomozygous_or_Heterozygous")
                header_printed = True

            if isinstance(most_freq, float) and isinstance(least_freq, float):
                print(f"{sample}\t{chrom}\t{pos}\t{vcf_genotype}\t{ref_alt}\t{most_freq:.4f}\t{least_freq:.4f}\t{hom_or_het}")
            else:
                print(f"{sample}\t{chrom}\t{pos}\t{vcf_genotype}\t{ref_alt}\t{most_freq}\t{least_freq}\t{hom_or_het}")

            if found_positions == positions_to_find:
                break

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_frequencies.py <sample_positions.txt> <variants.vcf[.gz]>")
        sys.exit(1)

    sample_file = sys.argv[1]
    vcf_file = sys.argv[2]

    sample_positions = parse_sample_file(sample_file)
    parse_vcf(vcf_file, sample_positions)
