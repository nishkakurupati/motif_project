import create_dataset_utils
import random
import parameters


# Keep a deterministic seed to allow for easy debug
random.seed(0)

# Run through combinations and create dataset
for k in range (10):
    for i in range (len(parameters.params_list)):
        ml = parameters.params_list[i][0]
        sc = parameters.params_list[i][1]
        sl = parameters.params_list[i][2]

        create_dataset_utils.create_files(ml, sc, sl, k)

