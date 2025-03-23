from sklearn.model_selection import train_test_split
import scipy.sparse as sparse


class Dataset(object):
    def __init__(self, path_file:str, coord_file:str):
        self.user_num = 0
        self.poi_num = 0
        self.training_file = path_file

        self.uidx_to_uid, self.uid_to_uidx = {}, {}
        self.pidx_to_pid, self.pid_to_pidx = {}, {}
        self.coord_file = coord_file


    def read_raw_data(self):
        # FIRST PASS - INDEXES
        with open(self.training_file, 'r') as f:
            for line in f:
                uid, pid, *_ = line.strip().split()
                uid = int(uid)
                pid = int(pid)
                if uid not in self.uid_to_uidx:
                    self.uid_to_uidx[uid] = self.user_num
                    self.uidx_to_uid[self.user_num] = uid
                    self.user_num += 1
                
                if pid not in self.pid_to_pidx:
                    self.pid_to_pidx[pid] = self.poi_num  
                    self.pidx_to_pid[self.poi_num] = pid
                    self.poi_num +=1

        # Second pass create the matrix
        sparse_raw_matrix = sparse.dok_matrix((self.user_num, self.poi_num))
        all_data = open(self.training_file, 'r').readlines()
        for eachline in all_data:
            parts = eachline.strip().split()  # Dividir la línea en columnas
            uid, lid, freq = parts[:3]  # Tomar solo las primeras tres columnas
            uid = int(uid)
            lid = int(lid)
            uidx, pidx, freq = self.uid_to_uidx[uid], self.pid_to_pidx[lid], int(float(freq))
            sparse_raw_matrix[uidx, pidx] = sparse_raw_matrix[uidx, pidx] + freq

        return sparse_raw_matrix.tocsr()

    '''
    def split_data(self, raw_matrix, random_seed=0):
        train_matrix = sparse.dok_matrix((self.user_num, self.poi_num))
        test_set = []
        for user_id in range(self.user_num):
            place_list = raw_matrix.getrow(user_id).indices
            freq_list = raw_matrix.getrow(user_id).data
            train_place, test_place, train_freq, test_freq = train_test_split(place_list, freq_list, test_size=0, random_state=random_seed)

            for i in range(len(train_place)):
                train_matrix[user_id, train_place[i]] = train_freq[i]
            test_set.append(test_place.tolist())

        return train_matrix.tocsr(), test_set
    '''
    def split_data(self, raw_matrix):
        train_matrix = sparse.dok_matrix((self.user_num, self.poi_num))
        test_set = []  

        for user_id in range(self.user_num):
            place_list = raw_matrix.getrow(user_id).indices
            freq_list = raw_matrix.getrow(user_id).data

            for i in range(len(place_list)):
                train_matrix[user_id, place_list[i]] = freq_list[i]

            test_set.append([])  # No habrá datos en el test_set

        return train_matrix.tocsr(), test_set

    def read_poi_coos(self):

        poi_coos = {}
        poi_data = open(self.coord_file, 'r').readlines()
        for eachline in poi_data:
            lid, lat, lng = eachline.strip().split()
            lid = int(lid)
            iidx, lat, lng = self.pid_to_pidx[lid], float(lat), float(lng)
            poi_coos[iidx] = (lat, lng)

        #Generate a list starting from 0 to NPOIs-1, so that each element of the list has the coordinates of iidx 
        place_coords = []
        sorted_keys = sorted(poi_coos.keys())
        for k in sorted_keys:
            lat, long = poi_coos[k]
            place_coords.append([lat, long])

        return place_coords

    def generate_data(self, random_seed=0):
        raw_matrix = self.read_raw_data()
        train_matrix, test_set = self.split_data(raw_matrix)
        place_coords =self.read_poi_coos()
        return train_matrix, test_set, place_coords


