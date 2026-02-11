# kraken2micro
Code to automatically create microanalyst input tables from kraken outputs

## How to use

python kraken2micro.py --rank 'S' --organism 'Bacteria' --files_dir input_dir/ --output_dir example_micro

## Parameters
--rank   Last rank for tables construction (default= 'S') (options = 'S', 'G', 'F', 'O', 'C', 'P', 'K')

--organism  Filter organisms that are not the one provided (default= 'Bacteria') (options = 'Bacteria', 'Viruses', 'Archaea', 'Fungi')

--files_dir directory with kraken2 mpa files (Must be with '.txt' extension)

--output_dir Location of the microbiome analysist format files
