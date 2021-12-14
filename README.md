# OpenCustomDB
Introduction
OpenCMS is a tool for the construction of custom (sample specific) protein sequence databases based on RNA-seq data. It uses the OpenProt annotations and the OpenVar variant annotator, thus includes canonical or reference proteins, novel isoforms, and alternative proteins encoded in non-canonical open reading frames larger than 29 codons.
Databases generated with OpenCMS can directly be used for protein database search in mass spectrometry (MS)-based proteomics experiments with any search engine.

By default, the 100,000 proteins encoded in the most abundant transcripts in the sample are annotated in the generated database. This makes it possible not to obtain a database that is too large. This value can be changed. Rather than limiting the total number of proteins, TPM values can be used instead. 

Inputs files:
OpenCustomDB needs two inputs:
Input 1: A VCF file generated with an identifier of genetic variants. Several variant detectors are available, including SAMtools, the GATK Best Practices pipeline, CTAT, FreeBayes, MuTect2, Strelka2, and VarScan2
Input 2: A quantification of transcript abundance generated by Kallisto-quant.
Options:
protnumber: number of proteins added in the database (value by default :100,000)
tpmnumber: use a tpm threshold instead of a protein number.
ip-iiban: altProts and novel isoforms are not added into the database (default = No)
trxexclude: requires a list of transcripts. Proteins encoded in these transcripts will not be added to the database. Format as shown below.

```
ENST00000378024.8
ENST00000322586.2
ENST00000346458.1
ENST00000313147.7
ENST00000300569.5
```

trsave: requires a list of transcripts. Proteins encoded in these transcripts will be added to the database even if the total number of proteins exceeds the threshold. Format as shown below.

```
ENST00000378024.8
ENST00000322586.2
ENST00000346458.1
ENST00000313147.7
ENST00000300569.5
```

Output:
One custom protein sequence database in fasta format.
One tab-file resuming all variant proteins found by the annotator: 

```
Prot	Transcrit	HGVS_P	HGVS_C	Potential_Error
IP_3379289	FBgn0013687	p.Lys11Asn	c.33A>T	WARNING_TRANSCRIPT_NO_START_CODON
```

One summary file in tab format:
```
effective_threshold	altProt (IP)	novel isoform (II)	refProt	altProt_variants (IP)	novel isoform_variants (II)	refProt_variants
0.19579	            63497	        1655	              18088	  16180	                580	                        0
```


Running OpenCustomDB with Python
Clone github repository https://github.com/NoeGuilloy/OpenCustomDatabase.

Launch from Python as shown below :

```
from OpenCMS.OpenCMS import OpenCMS                        #Load OpenCMS
OPCMS = OpenCMS(vcf_path=vcf,                              #Vcf Path
                        expname = name,                    #experience name
                        input_kallisto=kallisto,
                   ip_iiban =None,                         #yes = non-canonical proteins are excluded, default = no 
                   protnumber =None,                       #default = 100000
                   tpmnumber = None,                       #default = 0
                   trxexclude = None,                      #list of transcript to exclude
                   trsave = None)                          #list of transcript to rescue
    OPCMS.run()

```


 

