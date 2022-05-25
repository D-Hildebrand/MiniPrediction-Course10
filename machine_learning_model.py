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


def parq_model(df_train, df_valid):
    x_train = df_train[['conservationPro', 'conservationAla', 'conservationHis',
               'conservationThr', 'conservationGln', 'conservationTyr',
               'conservationGly', 'conservationArg', 'conservationVal',
               'consWildType', 'conservationGlu', 'conservationMet',
               'conservationLys', 'conservationIle', 'conservationPhe',
               'conservationLeu', 'conservationAsn', 'conservationSer',
               'conservationAsp', 'conservationCys', 'consVariant',
               'conservationTrp', 'score']]
    y_train = df_train['class']

    x_valid = df_valid[['conservationPro', 'conservationAla', 'conservationHis',
               'conservationThr', 'conservationGln', 'conservationTyr',
               'conservationGly', 'conservationArg', 'conservationVal',
               'consWildType', 'conservationGlu', 'conservationMet',
               'conservationLys', 'conservationIle', 'conservationPhe',
               'conservationLeu', 'conservationAsn', 'conservationSer',
               'conservationAsp', 'conservationCys', 'consVariant',
               'conservationTrp', 'score']]
    y_valid = df_valid['class']

    clf = RandomForestClassifier(n_estimators=100)

    clf.fit(x_train, y_train)

    y_pred = clf.predict(x_valid)

    print("Accuracy:", metrics.accuracy_score(y_valid, y_pred))
    return clf


def clf_test(df_test, clf):
    x_test = df_test[
        ['conservationPro', 'conservationAla', 'conservationHis',
         'conservationThr', 'conservationGln', 'conservationTyr',
         'conservationGly', 'conservationArg', 'conservationVal',
         'consWildType', 'conservationGlu', 'conservationMet',
         'conservationLys', 'conservationIle', 'conservationPhe',
         'conservationLeu', 'conservationAsn', 'conservationSer',
         'conservationAsp', 'conservationCys', 'consVariant',
         'conservationTrp', 'score']]
    y_pred = clf.predict(x_test)
    df_test['class'] = y_pred.tolist()
    df_test.to_csv('test_data_predictions.tsv', sep="\t")


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
    df_train = create_parq_df("train_data_bio_prodict.parq")
    df_valid = create_parq_df("valid_data_bio_prodict.parq")
    df_test = create_parq_df("test_data_bio_prodict.parq")
    clf = parq_model(df_train, df_valid)
    clf_test(df_test, clf)

    # file = "train_data_bio_prodict.tsv"
    # df = file_reader(file)
    # blosum_classifier(df)



main()
