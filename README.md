## **vcf_from_snvphyl**

This tool is essentially glue between the SNVPhyl workflow and *tbvcfreport*. It
reads the variant report from SNVPhyl and generates annotated (ala. SnpEff) VCF
files, using the COMBAT-TB eXplorer database for annotation. In the process
it computes the effect of variants on genes and proteins, thus the
including `snptools.vcf` module might be more broadly useful.
