
## Filter the Caendr vcf file for founders strains
bcftools view -S founders_list.txt WI.20250625.hard-filter.isotype.vcf.gz -O v -o WI.20250625.hard-filter.isotype.CemeeFounders.vcf.gz

## Keep only the polymorphic sites
bcftools view -c 1 WI.20250625.hard-filter.isotype.CemeeFounders.vcf.gz -O v -o WI.20250625.hard-filter.isotype.CemeeFounders_poly.vcf.gz