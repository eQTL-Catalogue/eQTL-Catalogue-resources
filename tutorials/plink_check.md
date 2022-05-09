
# Checking raw genotype files prior to imputation

## VCF to PLINK conversion

If your genotypes are in VCF format, then you need to first convert them to the binary PLINK format (.bed/.bim/.fam):

```bash
plink --vcf <path_to_vcf_file> --make-bed --out <plink_file_prefix>

```

## Checking the X chromosome name

The eQTL-Catalogue/genimpute workflow assumes that the input genotypes are in the binary plink format (.bed/.bim/.fam),
use GRCh37 coordinates and the name of the X chromosome is 'X' with PAR and non-PAR regions merged. You can use to following command to check the chromosome names in a plink file:

```bash
(base) [a72094@rocket plink_raw]$ cut -f 1 fetal_geno.bim | uniq -c
    625 0
 138820 1
 146994 2
 124540 3
 116400 4
 110321 5
 117265 6
  97959 7
  96179 8
  79584 9
  91288 10
  87854 11
  85541 12
  64298 13
  58756 14
  55913 15
  58850 16
  50339 17
  53909 18
  36098 19
  43963 20
  24952 21
  25849 22
  31538 23
    494 24
   1861 25
     46 26
```

In this case, the X chromsome genotypes are encoed as 23 (Non-PAR) and 25 (PAR). Thus, we need to merge the PAR and non-PAR regions and rename both to 'X":

```bash
plink -bfile fetal_geno --merge-x --make-bed --output-chr MT --out merged_x
```

And now the chromosome names are correct:

```bash
(base) [a72094@rocket plink_raw]$ cut -f 1 merged_x.bim | uniq -c
    625 0
 138820 1
 146994 2
 124540 3
 116400 4
 110321 5
 117265 6
  97959 7
  96179 8
  79584 9
  91288 10
  87854 11
  85541 12
  64298 13
  58756 14
  55913 15
  58850 16
  50339 17
  53909 18
  36098 19
  43963 20
  24952 21
  25849 22
  33399 X
    494 Y
     46 MT
```

# Checking for individuals with high levels of missingness

Optionally, you can also immediately check if there are some individual with many missing genotypes (see manual QC steps below). Individual samples with high levels of missingness (e.g. > 5%) should be excluded, because their presence can cause the imputation workflow to fail.

```bash
plink --bfile <plink_file_prefix> --missing
```
