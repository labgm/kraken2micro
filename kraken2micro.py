import pandas as pd
import argparse

"""
autor: TMDavi

data: 05/06/2024
objective: parse kraken2 mpa files to generate a otu table and a taxonomy table in the microbiome analyst format

input: list of mpa files
output: microanalyst_otu_table.txt and microanalyst_tax_table.txt

How to use

python kraken2micro.py --rank 'S' --organism 'Bacteria' --files MP1_A_mpa.txt MP1_B_mpa.txt	MP1_C_mpa.txt MP2_A_mpa.txt MP2_B_mpa.txt MP2_C_mpa.txt	MP3_A_mpa.txt MP3_B_mpa.txt MP3_C_mpa.txt MP4_A_mpa.txt MP4_B_mpa.txt MP4_C_mpa.txt MP5_A_mpa.txt MP5_B_mpa.txt MP5_C_mpa.txt
"""


def parseMPA(file_mpa):
    columns=['taxonomy','reads']
    df = pd.read_csv(file_mpa,sep='\t',header=None,names=columns)
    return df

def make_otu_list(df, rank):

    tax_dic = {'D': 'd__','K':'k__','P': 'p__', 'C': 'c__', 'O': 'o__', 'F': 'f__', 'G': 'g__', 'S': 's__'}
    
    lines = []
    index=0
    for l in df['taxonomy']:
        l = l.split('|')
        if l[-1][:3] == tax_dic[rank]: #Caso os tres primeiros caracteres do ultimo elemento da lista sejam iguais a rank correspondente, reversar a linha
            reads = df.loc[index,'reads']
            lines.append([l, reads])
        index+=1

    return lines

def parse_otu_taxonomy(otu_table):

    lines = []
    for l in otu_table['#NAME']:
        l = l.split(';')
        lines.append(l)
    return lines

def make_otu_table(kraken_reports, rank):
    otu_table = pd.DataFrame()
    temp_table = pd.DataFrame()
    for report in kraken_reports:
        df_mpa = parseMPA(report)
        otu_list = make_otu_list(df_mpa, rank)

        temp_table = pd.DataFrame({
            '#NAME':[';'.join(rank[0]) for rank in otu_list],
            f"{report.replace('.txt','')}":[rank[1] for rank in otu_list]
        })
        if otu_table.empty:
            otu_table = temp_table
        else:
            otu_table = pd.merge(otu_table, temp_table, on='#NAME', how='outer')
    otu_table = otu_table.fillna(0)
    return otu_table

def trimotutable(otu_table):
    lines = parse_otu_taxonomy(otu_table)
    otu_table['#NAME'] = [i[-1] for i in lines]

    return otu_table

def make_tax_table(otu_table, rank):

    otu_list = parse_otu_taxonomy(otu_table)
    
    table = []
    rank_lists = {
        'S': ["d__", "k__","p__", "c__", "o__", "f__", "g__", "s__"],
        'G': ["d__", "k__","p__", "c__", "o__", "f__", "g__"],
        'F': ["d__", "k__","p__", "c__", "o__", "f__"],
        'O': ["d__", "k__","p__", "c__", "o__"],
        'C': ["d__", "k__","p__", "c__"],
        'P': ["d__", "k__","p__"],
        'K': ["d__","k__",],
        'D': ["d__"]
    }
    tax_rank = rank_lists[rank][:]
    
    for line in otu_list:
        for r in line:
            for id in range(0, len(tax_rank)):
                if r[:3] == tax_rank[id]:
                    tax_rank[id] = r
        table.append(tax_rank[:])
        tax_rank = rank_lists[rank][:]

    tax_final = [[linha[i] for linha in table] for i in range(len(tax_rank))]

    possible_dicts = {
        'S': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[0], "Kingdom":tax_final[1],"Phylum": tax_final[2], "Class": tax_final[3],
              "Order": tax_final[4], "Family": tax_final[5], "Genus": tax_final[6],
              "Species": tax_final[7]},
        'G': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[0], "Kingdom":tax_final[1],"Phylum": tax_final[2], "Class": tax_final[3],
              "Order": tax_final[4], "Family": tax_final[5], "Genus": tax_final[6]},
        'F': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[0], "Kingdom":tax_final[1], "Phylum": tax_final[2], "Class": tax_final[3],
              "Order": tax_final[4], "Family": tax_final[5]},
        'O': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[0], "Kingdom":tax_final[1], "Phylum": tax_final[2], "Class": tax_final[3],
              "Order": tax_final[4]},
        'C': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[0], "Kingdom":tax_final[1], "Phylum": tax_final[2], "Class": tax_final[3]},
        'P': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[4], "Phylum": tax_final[5]},
        'K': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[0], "Kingdom":tax_final[1]},
        'D': {'#TAXONOMY': [rank[-1] for rank in otu_list],
              "Domain": tax_final[0]}
    }

    tax_table = pd.DataFrame.from_dict(possible_dicts[rank])

    return tax_table

def selectTaxa(otu_table, tax_table, org='Bacteria'):
    organisms = {'Bacteria':'d__Bacteria', 'Viruses':'d__Viruses','Archaea':'d__Archaea','Eukaryota':'d__Eukaryota','Fungi':'k__Fungi'}
    otu_table = trimotutable(otu_table)

    
    if org=='Fungi':
        tax_table = tax_table[tax_table['Kingdom'] == organisms[org]]
    else:
        tax_table = tax_table[tax_table['Domain'] == organisms[org]]

    otu_table = otu_table[otu_table['#NAME'].isin(tax_table['#TAXONOMY'])]

    return otu_table, tax_table


def main():
    parser = argparse.ArgumentParser(description="Process Kraken reports and generate OTU and taxonomy tables.")
    parser.add_argument('--rank', type=str, choices=['S', 'G', 'F', 'O', 'C', 'P', 'K'], required=True, help='Taxonomic rank to use.')
    parser.add_argument('--organism', type=str, choices=['Bacteria', 'Viruses', 'Archaea', 'Fungi'], required=True, help='Organism to filter.')
    parser.add_argument('--files', type=str, nargs='+', required=True, help='List of Kraken report files to parse.')

    args = parser.parse_args()

    otu_table = make_otu_table(args.files, args.rank)
    tax_table = make_tax_table(otu_table, args.rank)
    otu_table, tax_table = selectTaxa(otu_table, tax_table, org=args.organism)

    if args.organism != 'Fungi':
        tax_table = tax_table.drop(columns='Kingdom')
        tax_table = tax_table.rename(columns={'Domain': 'Kingdom'})
        
    if args.organism == 'Fungi':
        tax_table = tax_table.drop(columns='Domain')
    

    otu_table.to_csv('microanalyst_otu_table.txt', sep='\t',index=False)
    tax_table.to_csv('microanalyst_taxonomy_table.txt', sep='\t',index=False)


if __name__ == "__main__":
    main()


