# SNPs3DProfile

## 1. Introduction
SNPs3DProfile scores the entropic effects of single nucleotide polymorphisms (SNPs) to predict their deleterious effects  

## 2. Liscense
SNPs3DProfile is [open source?] and is released under the [lisence?] liscence

## 3. Prerequisites
**Operating System:** UNIX/Linux (64 bit) <br/>
**Software:** Perl v5.10.1 <br/>
**3rd Party Applications:** SVM-Light v6.02, BLAST+ 2.2.30 <br/>
**3rd Party Database:** BLAST database <version?> <br/>

## 4. Installation
Download link for SNPs3DProfile - [https://github.com/MoultLab/SNPs3DProfile](https://github.com/MoultLab/SNPs3DProfile)
### 4.1 Third party tools
The third party applications and database can be found at <link>.  The two applications should be placed in ```/path/to/SNPs3DProfile/third_party_apps``` while the database should be placed in ```/path/to/SNPs3DProfile/profile_data```

## 5. Running Script
SNPs3DProfile can now be run as follows: <br/>
**1.** Go to the SNPs3DProfile 'bin' directory ``` $ cd /path/to/SNPs3DProfile/bin``` <br/>
**2.** Run script with appropriate arguments ``` $ perl profile_pipeline.pl --options <arg>``` <br/>

All arguments are given as strings unless otherwise stated.  Arguments include: <br/>
```--dbserver          Database server``` <br/>
```--dbuser            Database user``` <br/>
```--dbpass            Database password``` <br/>
```--path              File path to script.  Does not extend to the 'bin' directory (i.e. should end with SNPs3DProfile)```<br/>
```--list              Text file with mutations.  Each line should include the mutation name, sequence NP number, and mutation separated by a tab.  See example file for format ```<br/>
  Optional arguments include:  <br/>
```--mysqldb           SQL database.  Set to "SNPs3D_Dev" by default ```<br/>
```--table             SQL table.  Set to "SNPs3D_Profile_2017_precomputed_psimtx_500_1_0_JAN2018" by default ```<br/>
```--blastbin          BLAST bin.  Default assumes third party bin was downloaded as described in Installation. ```<br/>
```--blastdb           BLAST database.  See --blastbin.  ```<br/>
```--rpb               Run PSI-BLAST.  Default integer value set to 1 to run on missing sequences; set to 0 otherwise   ```<br/>
```--output            Output file.  Set to "pred_HWUL2.txt" by default ```<br/>
```--log               Log file.  Set to "log.profile_pipeline" by default ```<br/>
       
## 6. Contributions
[add contributions]
