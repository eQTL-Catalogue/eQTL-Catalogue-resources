
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

## Checking for individuals with high levels of missingness

Optionally, you can also immediately check if there are some individual with many missing genotypes (see manual QC steps below). Individual samples with high levels of missingness (e.g. > 5%) should be excluded, because their presence can cause the imputation workflow to fail.

```bash
plink --bfile <plink_file_prefix> --missing
```

We can then sort the the output file by the 6th column.

```bash
(base) [a72094@rocket plink_raw]$ sort -k6n plink.imiss | tail -n 10
 565     566          Y    17028   759993  0.02241
 553     554          Y    19595   759993  0.02578
 241     242          Y    22634   752646  0.03007
 155     155          Y    22967   759993  0.03022
1086   155_2          Y    24584   759993  0.03235
 931     932          Y    43064   752646  0.05722
 919     920          Y    45395   752646  0.06031
 913     914          Y    73313   752646  0.09741
 925     926          Y    90780   752646   0.1206
 844     845          Y   176952   752646   0.2351
 ```
 In this case, we can see that there are five samples with more than 5% of the genotypes missing.
 
 To remove those samples, we first need to make text file with corresponding sample ids:

```bash
(base) [a72094@rocket plink_raw]$ cat remove_list.txt
931 932
919 920
913 914
925 926
844 845
```

We can then use the --remove option to remove those samples from the plink file

```bash
plink -bfile OneK1K_AllChr --remove remove_list.txt --make-bed --out OneK1K_nonmissing
```

## Avoding numerical sample ids

While numerical sample ids are allowed both in plink and VCF formats, they can cause potential downstream problems when loading genotype matrices into R. Thus it's a good practice to make sure that genotype ids start with a character. We also do not need family ids (fid) in our analysis.

If the original .fam file looks something like this:
```bash
(base) [a72094@rocket plink_raw]$ head -n5 OneK1K_AllChr.fam
1 1 0 0 2 -9
2 2 0 0 2 -9
3 3 0 0 2 -9
4 4 0 0 2 -9
5 5 0 0 2 -9
```

We could change it to something like this:

```bash
(base) [a72094@rocket plink_raw]$ head -n5 OneK1K_nonmissing.fam
0 OneK1K_1 0 0 2 -9
0 OneK1K_2 0 0 2 -9
0 OneK1K_3 0 0 2 -9
0 OneK1K_4 0 0 2 -9
0 OneK1K_5 0 0 2 -9
```
