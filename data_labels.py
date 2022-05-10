import pandas as pd
import re
import blosum as bl
import csv


def create_dataframe(input_parq, output_filename):
    # in dataset kolom class -> 1 = pathogeen (ziek), 0 = benign
    df = pd.read_parquet(input_parq, engine='pyarrow')
    df_rownames = df.index.values

    with open(output_filename, "a") as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')

        for i in range(len(df_rownames)):
            aa = re.split("[0-9]+", df_rownames[i])[1:]
            ref = three_to_one(aa[0])
            alt = three_to_one(aa[1])
            label = df.iloc[i, -1]
            score = bl.BLOSUM(62)[ref+alt]
            tsv_writer.writerow([ref, alt, label, score])

    out_file.close()


def three_to_one(three_letter_code):
    aa_list = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
          'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return aa_list[three_letter_code]


def create_tsv(input_parq):
    output_filename = input_parq.replace("parq", "tsv")

    with open(output_filename, 'w') as out_file:
        output_tsv = csv.writer(out_file, delimiter='\t')
        output_tsv.writerow(["REF", "ALT", "Label", "Score"])
    out_file.close()

    return output_filename


def main():
    input_files = ["test_data_bio_prodict.parq", "train_data_bio_prodict.parq", "valid_data_bio_prodict.parq"]
    for input_parq in input_files:
        output_filename = create_tsv(input_parq)
        create_dataframe(input_parq, output_filename)


main()
