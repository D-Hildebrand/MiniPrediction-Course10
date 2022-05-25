from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from sklearn import metrics
import numpy as np
import re
import blosum as bl


def file_reader(file):
    """
    Converts a TSV file to a pandas dataframe
    :param file: Input .TSV file
    :return df: Pandas dataframe
    """
    return pd.read_table(file)


def three_to_one(three_letter_code):
    aa_list = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return aa_list[three_letter_code]


def create_parq_df(input_parq):
    df = pd.read_parquet(input_parq, engine='pyarrow')
    df_rownames = df.index.values
    for i, row in df.iterrows():
        aa = re.split("[0-9]+", i)
        ref = three_to_one(aa[1])
        alt = three_to_one(aa[2])
        score = bl.BLOSUM(62)[ref + alt]
        df.at[i, 'score'] = score

    return df


def parq_model(df):
    parq = df[['conservationPro', 'conservationAla', 'conservationHis',
               'conservationThr', 'conservationGln', 'conservationTyr',
               'conservationGly', 'conservationArg', 'conservationVal',
               'consWildType', 'conservationGlu', 'conservationMet',
               'conservationLys', 'conservationIle', 'conservationPhe',
               'conservationLeu', 'conservationAsn', 'conservationSer',
               'conservationAsp', 'conservationCys', 'consVariant',
               'conservationTrp', 'score']]
    pathogenicity = df['class']

    x_train, x_test, y_train, y_test = train_test_split(parq, pathogenicity,
                                                        test_size=0.25)

    clf = RandomForestClassifier(n_estimators=100)

    clf.fit(x_train, y_train)

    y_pred = clf.predict(x_test)

    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))


def blosum_classifier(df):
    blosum = np.vstack(df['Score'])
    pathogenicity = df['Label']

    x_train, x_test, y_train, y_test = train_test_split(blosum, pathogenicity,
                                                        test_size=0.25)

    clf = RandomForestClassifier(n_estimators=100)

    clf.fit(x_train, y_train)

    y_pred = clf.predict(x_test)

    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))


def main():
    # input_files = ["train_data_bio_prodict.parq",
    #                "valid_data_bio_prodict.parq"]
    # for input_parq in input_files:
    df = create_parq_df("train_data_bio_prodict.parq")
    parq_model(df)

    # file = "train_data_bio_prodict.tsv"
    # df = file_reader(file)
    # blosum_classifier(df)



main()
