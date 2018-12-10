# FastEPRR working version
Program for estimating the population recombination rate written by Feng Gao.

See http://www.picb.ac.cn/evolgen/softwares/FastEPRR.html for more detailed information.

Linux/Unix
There are two ways to install FastEPRR package in Linux/Unix. We recommend the first way,
since it is convenient and it avoids opening the R console.
[1] Command Line
Once the R environment has been installed, the user can type the following command in the
command line to install FastEPRR package. That is,
$ R CM D INSTALL "package_path/FastEPRR_1.0.tar.gz"
[2] Within R console
Typing “R” in the command line can start the R console. That is,
$ R
The results of the above command show the basic information about R version, platform and
introductions and so on. Then in the R console, the user can type the following command at the R
prompt (the > symbol) to install FastEPRR package.
>install.packages("package_path /FastEPRR_1.0.tar.gz", repos = NULL, type="source")
The package_path is the file path of FastEPRR_1.0.tar.gz
