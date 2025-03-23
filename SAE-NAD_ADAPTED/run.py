from argparse import ArgumentParser
import heapq
import numpy as np
import scipy
import pandas as pd

import eval_metrics
import dataset
from model import AutoEncoder

import torch
from torch.autograd import Variable

from dataset import Dataset
import numpy as np
import scipy
from sklearn.metrics.pairwise import rbf_kernel
from scipy.sparse import csr_matrix

if torch.cuda.is_available():
    import torch.cuda as T
else:
    import torch as T





def get_mini_batch(X, weight_mask, batch_user_index):
    batch_item_index = []
    for user_id in batch_user_index:
        batch_item_index.append(X.getrow(user_id).indices)
    return X[batch_user_index].toarray(), weight_mask[batch_user_index].toarray(), batch_item_index


def log_surplus_confidence_matrix(B, alpha, epsilon):
    # To construct the surplus confidence matrix, we need to operate only on the nonzero elements.
    # This is not possible: S = alpha * np.log(1 + B / epsilon)
    S = B.copy()
    S.data = alpha * np.log(1 + S.data / epsilon)
    return S


def train_autoencoder(train_matrix, place_coords:list, data: Dataset, test_set: str, original_training_set: str, result_file:str, number_items: int):
    num_users, num_items = train_matrix.shape
    weight_matrix = log_surplus_confidence_matrix(train_matrix, alpha=args.alpha, epsilon=args.epsilon)
    train_matrix[train_matrix > 0] = 1.0
    place_correlation = cal_place_pairwise_dist(place_coords, args.gamma)

    assert num_items == place_correlation.shape[0]
    print(train_matrix.shape)

    # Construct the model by instantiating the class defined in model.py
    model = AutoEncoder(num_items, args.inner_layers, num_items, da=args.num_attention, dropout_rate=args.dropout_rate)
    if torch.cuda.is_available():
        print("CUDA")
        model.cuda()

    criterion = torch.nn.MSELoss(size_average=False, reduce=False)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay)

    batch_size = args.batch_size
    user_indexes = np.arange(num_users)

    model.train()
    for t in range(args.epoch):
        print("epoch:{}".format(t))
        np.random.shuffle(user_indexes)
        avg_cost = 0.
        for batchID in range(int(num_users / batch_size)):
            start = batchID * batch_size
            end = start + batch_size

            batch_user_index = user_indexes[start:end]

            batch_x, batch_x_weight, batch_item_index = get_mini_batch(train_matrix, weight_matrix, batch_user_index)
            batch_x_weight += 1
            batch_x = Variable(torch.from_numpy(batch_x).type(T.FloatTensor), requires_grad=False)

            y_pred = model(batch_item_index, place_correlation)

            # Compute and print loss
            batch_x_weight = Variable(torch.from_numpy(batch_x_weight).type(T.FloatTensor), requires_grad=False)
            loss = (batch_x_weight * criterion(y_pred, batch_x)).sum() / batch_size

            print(batchID, loss.data)

            # Zero gradients, perform a backward pass, and update the weights.
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            avg_cost += loss / num_users * batch_size

        print("Avg loss:{}".format(avg_cost))

        # print the prediction score for the user 0
        print(model(
            [train_matrix.getrow(0).indices],
            place_correlation)[:, T.LongTensor(train_matrix.getrow(0).indices.astype(np.int32))])
        print(model([train_matrix.getrow(0).indices], place_correlation))

    # Evaluation
    model.eval()

    test_df = pd.read_csv(test_set, header=None, names=["user_id", "poi_id", "score", "timestamp"], sep="\t")
    user_test_ids = set(test_df["user_id"])

    training_no_impu = pd.read_csv(original_training_set, header=None, names=["user_id", "poi_id", "score", "timestamp"], sep="\t")
    rec_list = open(result_file, 'w')
    print(f"Users in test set: {len(user_test_ids)}")
    for uidx in range(num_users):
        real_user_id = data.uidx_to_uid[uidx]
        if real_user_id in user_test_ids:

            user_rating_vector = train_matrix.getrow(uidx).toarray()
            pred_rating_vector = model([train_matrix.getrow(uidx).indices], place_correlation)
            pred_rating_vector = pred_rating_vector.cpu().data.numpy()
            user_rating_vector = user_rating_vector[0]
            pred_rating_vector = pred_rating_vector[0]
            # pred_rating_vector[user_rating_vector > 0] = 0 # Consumed in training

            item_recommended_dict = dict()
            for item_inner_id, score in enumerate(pred_rating_vector):
                item_recommended_dict[item_inner_id] = score

            sorted_items_list = sorted(item_recommended_dict.items(), key=lambda x: x[1], reverse=True)

            counter = 0
            poi_ids_consumed = set(training_no_impu.loc[training_no_impu["user_id"] == real_user_id, "poi_id"])
            for key, value in sorted_items_list:
                real_poi = data.pidx_to_pid[key]
                if real_poi not in poi_ids_consumed:
                    rec_list.write(f"{real_user_id}\t{real_poi}\t{value}\n")
                    counter +=1
                    if counter >= number_items: # equivalent to top k
                        break


    rec_list.close()





    '''
    precision, recall, MAP = [], [], []
    for k in [5, 10, 15, 20]:
        precision.append(eval_metrics.precision_at_k(test_set, recommended_list, k))
        recall.append(eval_metrics.recall_at_k(test_set, recommended_list, k))
        MAP.append(eval_metrics.mapk(test_set, recommended_list, k))

    print(precision)
    print(recall)
    print(MAP)
    '''

def cal_place_pairwise_dist(place_coordinates, gamma):
    # this method calculates the pair-wise rbf distance
    place_correlation = rbf_kernel(place_coordinates, gamma=gamma)
    np.fill_diagonal(place_correlation, 0)
    place_correlation[place_correlation < 0.1] = 0
    place_correlation = csr_matrix(place_correlation)
    return place_correlation

if __name__ == "__main__":
    parser = ArgumentParser(description="SAE-NAD")
    parser.add_argument('--epoch', type=int, default=10, help='number of epochs for GAT')
    parser.add_argument('--batch_size', type=int, default=256, help='batch size for training')
    parser.add_argument('--alpha', type=float, default=2.0, help='the parameter of the weighting function')
    parser.add_argument('--epsilon', type=float, default=1e-5, help='the parameter of the weighting function')
    parser.add_argument('--learning_rate', type=float, default=1e-3, help='learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-3, help='weight decay')
    parser.add_argument('--num_attention', type=int, default=40, help='the number of dimension of attention')
    parser.add_argument('--inner_layers', nargs='+', type=int, default=[200, 50, 200], help='the number of latent factors')
    parser.add_argument('--dropout_rate', type=float, default=0.5, help='the dropout probability')
    # parser.add_argument('--seed', type=int, default=0, help='random state to split the data')
    parser.add_argument('--training_file', type=str, default="GB_London_K2_AgTLAST_APSUM_TTrain.dat",help='training_file')
    parser.add_argument('--original_training_set', type=str, default="GB_London_K2_AgTLAST_APSUM_TTrain.dat", help='original training non imputed)')
    parser.add_argument('--coord_file', type=str, default="GB_LondonPOIS_Coords.txt", help='coord_file')
    parser.add_argument('--test_set', type=str, default="GB_London_K2_AgTLAST_APSUM_TTest.dat", help='test_set')
    parser.add_argument("--result_file", type=str, default="salida.txt", help="Path of the result file.")
    parser.add_argument("--nI", type=int, default=100, help="number items to recommend.")
    parser.add_argument("--gamma", type=int, default=90, help="gamma")
    args = parser.parse_args()
    for key, value in vars(args).items():
        print(f"{key}: {value}")
        
    # try attention model
    data_set = Dataset(args.training_file, args.coord_file)
    train_matrix, test_matrix, place_coords = data_set.generate_data()
    train_autoencoder(train_matrix, place_coords, data_set, args.test_set, args.original_training_set, args.result_file, args.nI)
