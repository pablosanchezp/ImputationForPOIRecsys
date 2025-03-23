import pandas as pd
import argparse
from model_imputation import EASE


def main():
    parser = argparse.ArgumentParser(
        description="launch ease")
    '''
    parser.add_argument("--training_imputation", type=str, help="Path of the imputed training file.", default="LastfmHetrecRandomGlobal8020TrainImputedExample.txt")
    parser.add_argument("--original_training", type=str, help="Path of the original training file.", default="LastfmHetrecRandomGlobal8020Train.txt")
    parser.add_argument("--test", type=str, help="Path of the test file.", default="LastfmHetrecRandomGlobal8020Test.txt")
    parser.add_argument("--implicit", type=bool, help="Boolean for implicit.", default="True")
    parser.add_argument("--result", type=str, help="Path of the result file.", default="salida.txt")
    parser.add_argument("--lamb", type=float, help="lamda.", default=0.05)
    parser.add_argument("--nI", type=int, help="numberItems.", default="50")
    '''

    parser.add_argument("--training_imputation", type=str, help="Path of the imputed training file.")
    parser.add_argument("--original_training", type=str, help="Path of the original training file (NO IMPUTATION).")
    parser.add_argument("--test", type=str, help="Path of the test file.")
    parser.add_argument("--implicit", type=bool, help="Boolean for implicit.", default=False)
    parser.add_argument("--result", type=str, help="Path of the result file.")
    parser.add_argument("--lamb", type=float, help="lamda.")
    parser.add_argument("--nI", type=int, help="numberItems.")


    args = parser.parse_args()

    training_imputation_df = pd.read_csv(args.training_imputation, header=None, sep="\t")
    columns = training_imputation_df.columns.tolist()
    training_imputation_df.columns = ['user_id', 'item_id', 'rating'] + columns[3:]

    original_training_df = pd.read_csv(args.original_training, header=None, sep="\t")
    columns_original_training = original_training_df.columns.tolist()
    original_training_df.columns = ['user_id', 'item_id', 'rating'] + columns_original_training[3:]


    test_df = pd.read_csv(args.test, header=None, sep="\t")
    unique_users = test_df[0].unique()

    ease_rec = EASE()
    ease_rec.fit(training_imputation_df, args.lamb, args.implicit)
    df_result = ease_rec.predict(training_imputation_df, original_training_df, unique_users, training_imputation_df['item_id'].unique(), args.nI)
    df_result.to_csv(args.result, index=False, header=None, sep="\t")


if __name__ == "__main__":
    main()
