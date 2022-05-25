from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from sklearn import metrics
import re
import blosum as bl


def three_to_one(three_letter_code):
    """
    Converts the 3 letter protein codes to the corresponding 1 letter code
    :param three_letter_code: 3 lettered protein code
    :return: 1 letter protein code
    """
    aa_list = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return aa_list[three_letter_code]


def create_parq_df(input_parq):
    """
    Creates a pandas dataframe out of the parq file
    :param input_parq: Parq files provided by BioProdict
    :return df: Pandas dataframe containing information from the parq files
    """
    df = pd.read_parquet(input_parq, engine='pyarrow')
    for i, row in df.iterrows():
        aa = re.split("[0-9]+", i)
        ref = three_to_one(aa[1])
        alt = three_to_one(aa[2])
        score = bl.BLOSUM(62)[ref + alt]
        df.at[i, 'score'] = score

    return df


def parq_model(df_train, df_valid):
    """
    Creates the random forest classifier
    :param df_train: Dataframe with training parq data + blosum scores
    :param df_valid: Dataframe with validation parq data + blosum scores
    :return clf: Random forest classifier created by using the training and
                 validation sets
    """
    df_columns = ['conservationPro', 'conservationAla', 'conservationHis',
         'conservationThr', 'conservationGln', 'conservationTyr',
         'conservationGly', 'conservationArg', 'conservationVal',
         'consWildType', 'conservationGlu', 'conservationMet',
         'conservationLys', 'conservationIle', 'conservationPhe',
         'conservationLeu', 'conservationAsn', 'conservationSer',
         'conservationAsp', 'conservationCys', 'consVariant',
         'conservationTrp', 'score']

    # Retrieves the train and validation information and labels
    x_train = df_train[df_columns]
    y_train = df_train['class']
    x_valid = df_valid[df_columns]
    y_valid = df_valid['class']

    # Creates the random forest classifier
    clf = RandomForestClassifier(n_estimators=100)

    # Fits the model
    clf.fit(x_train, y_train)

    # Predicts the pathogenicity
    y_pred = clf.predict(x_valid)

    print("Accuracy:", metrics.accuracy_score(y_valid, y_pred))
    return clf


def clf_test(df_test, clf):
    """
    Classifier that predicts the pathogenicity of the test data
    :param df_test: Dataframe with test parq data + blosum scores
    :param clf: Random forest classifier created by using the training and
           validation sets
    :return test_data_prediction.tsv: File with predictions added to the test
            data
    """
    x_test = df_test[
        ['conservationPro', 'conservationAla', 'conservationHis',
         'conservationThr', 'conservationGln', 'conservationTyr',
         'conservationGly', 'conservationArg', 'conservationVal',
         'consWildType', 'conservationGlu', 'conservationMet',
         'conservationLys', 'conservationIle', 'conservationPhe',
         'conservationLeu', 'conservationAsn', 'conservationSer',
         'conservationAsp', 'conservationCys', 'consVariant',
         'conservationTrp', 'score']]

    # Predicts the pathogenicity
    y_pred = clf.predict(x_test)

    # Adds the prediction to the output file
    y_pred.tolist().to_csv('test_data_predictions.tsv', sep="\t")


def main():
    # Creates the train, validation and test pandas dataframes
    df_train = create_parq_df("train_data_bio_prodict.parq")
    df_valid = create_parq_df("valid_data_bio_prodict.parq")
    df_test = create_parq_df("test_data_bio_prodict.parq")

    # Creates the random forest classifier
    clf = parq_model(df_train, df_valid)

    # Uses the random forest classifier to predict the data from the
    # test set
    clf_test(df_test, clf)


main()
