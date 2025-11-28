This folder contains reference data from literature kindly shared with the pyMBE-dev team by its original authors.
We use the data to periodically benchmark pyMBE to ensure than future releases are still able to reproduce the existing results.
If you intend to use this data for your research, please acknowledge its original authors.
A list of files and its corresponding sources in the literature follows.

## Files in parent folder
These files have been standarized to facilitate accessing the data within the pyMBE framework. 
Minor changes in the data have also been done in some cases to correct for small errors in the original data reported by its authors.
The changes done in the original data can be checked in `mantainer/standarize_data.py`, which is the same script used to generate the data in this folder.

List of files:
* [^1] Blanco2020a.csv       
* [^2] Landsgesell2020a.csv
* [^3] Landsgesell2022a.csv  
* [^4] Lunkad2021a.csv
* [^5] Lunkad2021b.csv
* [^6] Torres2017.csv   
* [^7] Torres2022.csv
                 
## Files in src
The files in the `src` subfolder have not been modified as compared to what was shared with the pyMBE-dev team.
This means that the data, formatting and filename is kept unaltered.
The folder is chaotic due to its very nature, and therfore the files in this subfolder can be hard to navigate without prior knowledge of its contents.
To ease such endevour, we provide below a map between the files in the parent folder and those in `src`.

List of files:
* 1beb-10mM-torres.dat: Torres2017.csv
* 1f6s-10mM-torres.dat: Torres2022.csv  
* data_landsgesell.csv: Landsgesell2020a.csv           
* equilibrium_values_gel_MD.csv: Landsgesell2022a.csv  
* Glu-HisMSDE.csv: Lunkad2021b.csv           
* histatin5_SoftMatter.txt: Blanco2020a.csv  
* Lys-AspMSDE.csv: Lunkad2021a.csv
* weak-gel_total_data.csv: Landsgesell2022a.csv

[^1]: P. M. Blanco, S. Madurga, J. L. Garcés, F. Mas, R. S. Dias. Soft Matter(2021), 17(3), 655-669. doi: [10.1039/D0SM01475C](https://doi.org/10.1039/D0SM01475C).
[^2]: J. Landsgesell, P. Hebbeker, O. Rud, R. Lunkad, P. Kosovan, C. Holm.  Macromolecules (2020), 53(8), 3007-3020. doi: [10.1021/acs.macromol.0c00260](https://doi.org/10.1021/acs.macromol.0c00260).
[^3]: J. Landsgesell, D. Beyer, P. Hebbeker, P. Košovan, C. & Holm. Macromolecules (2022), 55(8), 3176-3188. doi: [10.1021/acs.macromol.1c02489](https://doi.org/10.1021/acs.macromol.1c02489).
[^4]: R. Lunkad, A. Murmiliuk, P. Hebbeker, M. Boublík, Z. Tošner, M. Štěpánek, P. Košovan. Molecular Systems Design & Engineering (2021), 6(2), 122-131. doi: [10.1039/D0ME00147C](https://doi.org/10.1039/D0ME00147C).
[^5]: R. Lunkad, A. Murmiliuk, Z. Tošner, M. Štěpánek, P. Košovan. Polymers (2021), 13(2), 214. doi: [10.3390/polym13020214](https://doi.org/10.3390/polym13020214).
[^6]: P. B. Torres, L. Bojanich, F. Sanchez-Varretti, A. J. Ramirez-Pastor, E. Quiroga, V. Boeris, C. F. Narambuena. Colloids and Surfaces B: Biointerfaces (2017), 160, 161-168. doi: [10.1016/j.colsurfb.2017.09.018](https://doi.org/10.1016/j.colsurfb.2017.09.018).
[^7]: P. B. Torres, P. M. Blanco, J. L. Garcés, C. F. Narambuena.  The Journal of Chemical Physics (2022), 157(20). doi: [10.1063/5.0122275](https://doi.org/10.1063/5.0122275).