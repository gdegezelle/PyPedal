# References

These are the published works used by the methods this manual describes.
Paper PDFs are copyrighted and are not committed in the repository.

## Methods

Ballou, J. D. 1997. Ancestral inbreeding only minimally affects
inbreeding depression in mammalian populations. *Journal of Heredity*
88:169–178.

Baumung, R., J. Farkas, D. Boichard, G. Mészáros, J. Sölkner, and
I. Curik. 2015. GRAIN: a computer program to calculate ancestral and
partial inbreeding coefficients using a gene dropping approach.
*Journal of Animal Breeding and Genetics* 132:100–108.
DOI 10.1111/jbg.12145.

Boichard, D., L. Maignel, and É. Verrier. 1997. The value of using
probabilities of gene origin to measure genetic variability in a
population. *Genetics Selection Evolution* 29:5–23.
DOI 10.1186/1297-9686-29-1-5.

Lacy, R. C. 1989. Analysis of founder representation in pedigrees:
founder equivalents and founder genome equivalents. *Zoo Biology*
8:111–123.

MacCluer, J. W., J. L. VandeBerg, B. Read, and O. A. Ryder. 1986.
Pedigree analysis by computer simulation. *Zoo Biology* 5:147–160.

Meuwissen, T. H. E., and Z. Luo. 1992. Computing inbreeding coefficients
in large populations. *Genetics Selection Evolution* 24:305–313.

Mrode, R. A. 2005. *Linear Models for the Prediction of Animal Breeding
Values*. 2nd ed. Wallingford, UK: CAB International. Table 2.1 is the
six-animal worked pedigree used in this manual. Appendix B.2 presents
the modified Meuwissen–Luo algorithm.

Pattie, W. 1965. Selection for weaning weight in Merino sheep.
*Australian Journal of Experimental Agriculture and Animal Husbandry*
5:353–360. Cited only as historical context for the stored `gencoeff` /
pedformat `p` field. PyPedal 4.0 does not compute Pattie coefficients.

Quaas, R. L. 1995. *Fx* algorithms, as presented in Mrode (2005)
Appendix B.2 (`method="mod_meu_luo"`).

Suwanlee, S., R. Baumung, J. Sölkner, and I. Curik. 2007. Evaluation of
ancestral inbreeding coefficients: Ballou’s formula versus gene
dropping. *Conservation Genetics* 8:489–495.
DOI 10.1007/s10592-006-9187-9.

VanRaden, P. M. 1992. Accounting for inbreeding and crossbreeding in
genetic evaluation of large populations. *Journal of Dairy Science*
75:3136–3144. Pedigree algorithm behind
`pyp_nrm.inbreeding(method="vanraden")`.

VanRaden, P. M. 2008. Efficient methods to compute genomic predictions.
*Journal of Dairy Science* 91:4414–4423. Method 1 of the genomic
relationship matrix. Not the 1992 pedigree method.

Wright, S. 1922. Coefficients of inbreeding and relationship.
*The American Naturalist* 56:330–338.

## PyPedal and project citations

Cole, J. B. 2007. PyPedal: a computer program for pedigree analysis.
*Computers and Electronics in Agriculture* 57:107–113.

Cole, J. B. 2012. A Manual for use of PyPedal: A software package for
pedigree analysis. Animal Improvement Programs Laboratory, Agricultural
Research Service, United States Department of Agriculture.

Cole, J. B., D. E. Franke, and E. A. Leighton. 2004. Population structure
of a colony of dog guides. *Journal of Animal Science* 82:2906–2912.

This PyPedal 4 manual is a newly written successor for the Python 3
line. It is not a reprint or relicensed edition of the 2012 USDA manual.
Scientific facts, method names, and citations above are used as
reference information. See [Notices](notices.md).
