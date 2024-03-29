# SNPs3DProfile [Packaging under progress]

## 1. Introduction
SNPs3DProfile scores the entropic effects of single nucleotide polymorphisms (SNPs) to predict their deleterious effects  
For more details about the package, contact John Moult at [jmoult@tunc.org].

## 2. Liscense
SNPs3DProfile is open source and is released under the [insert] liscence

## 3. Prerequisites
**Operating System:** UNIX/Linux (64 bit) <br/>
**Software:** Perl v5.10.1 <br/>
**3rd Party Applications:** SVM-Light v6.02, BLAST+ 2.10.1 (or later) <br/>
**3rd Party Database:** BLAST database 2.9.0 <br/>

## 4. Installation
Download link for SNPs3DProfile - [https://github.com/MoultLab/SNPs3DProfile](https://github.com/MoultLab/SNPs3DProfile)
### 4.1 Third party tools
The third party applications and database can be found at ```/moulthome/choe/third_party_apps.tar.gz``` and ```/moulthome/choe/profile_data.tar.gz```.  The two applications should be placed in ```/path/to/SNPs3DProfile/third_party_apps``` while the database should be placed in ```/path/to/SNPs3DProfile/profile_data```

## 5. Running Script
SNPs3DProfile can now be run as follows: <br/>
**1.** Go to the SNPs3DProfile 'bin' directory ``` $ cd /path/to/SNPs3DProfile/bin``` <br/>
**2.** Run script with appropriate arguments ``` $ perl profile_pipeline.pl --options <arg>``` <br/>

All arguments are given as strings unless otherwise stated.  Arguments include: <br/>
```--dbserver          Database server``` <br/>
```--mysqldb           SQL database.  ```<br/>
```--dbuser            Database username``` <br/>
```--dbpass            Database password``` <br/>
```--path              File path to script.  Does not extend to the 'bin' directory (i.e. should end with /SNPs3DProfile)```<br/>
```--list              File with mutations.  Each line should include the mutation name, sequence NP number, and mutation separated by a tab.  See example file for format ```<br/>
  Optional arguments include:  <br/>
```--table             SQL table.  Set to "SNPs3D_Profile_2017_precomputed_psimtx_500_1_0_JAN2018" by default ```<br/>
```--blastbin          BLAST bin.  Default assumes third party bin was downloaded as described in Installation. ```<br/>
```--blastdb           BLAST database.  See --blastbin.  ```<br/>
```--rpb               Run PSI-BLAST.  Default integer value set to 1 to run on missing sequences; set to 0 otherwise   ```<br/>
```--output            Output file.  Set to "pred_HWUL2.txt" by default ```<br/>
```--log               Log file.  Set to "log.profile_pipeline" by default ```<br/>
An example set of 300 sequences and its subsequent output are given under the testdata directory.
       
## 6. Contributions
Ezra Cho - tested pipeline, packaged dependencies, shipped code to GitHub, generated testset, and benchmarked code

## 7. Publication
Yue, P., &amp; Moult, J. (2006). Identification and Analysis of Deleterious Human SNPs. Journal of Molecular Biology, 356(5), 1263-1274. [doi:10.1016/j.jmb.2005.12.025](https://pubmed.ncbi.nlm.nih.gov/16412461/) <br/>
Yue, P., Melamud, E., & Moult, J. (2006).  SNPs3D: candidate gene and SNP selection for association studies.  BMC Bioinformatics, 7, 166. 
 [doi:doi.org/10.1186/1471-2105-7-166](https://pubmed.ncbi.nlm.nih.gov/16551372/)

