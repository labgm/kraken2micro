# kraken2micro
Code to automatically create microanalyst input tables from kraken outputs

## How to use

python kraken2micro.py --rank 'S' --organism 'Bacteria' --files file1_mpa.txt file2_mpa.txt	file3_mpa.txt

## Parameters
--rank   Last rank for tables construction (default= 'S') (options = 'S', 'G', 'F', 'O', 'C', 'P', 'K')

--organism  Filter organisms that are not the one provided (default= 'Bacteria') (options = 'Bacteria', 'Viruses', 'Archaea', 'Eukaryota')

--files list of kraken2 mpa files

