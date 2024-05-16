## \file correction_module.py
#
# \brief Given IC, find velocity correction according to neural network.

import numpy as np
from sklearn.preprocessing import MinMaxScaler
import joblib   # To save/load scaler
import tensorflow as tf
import keras


DIM = 6     # dimension of RTBP system

## Given IC, find correction maneuver according to neural network.
#
# @param[in]  q 	init. cond. in position-velocity (x,y,z,dx,dy,dz)
# @param[out] dv 	modulus of correction maneuver (in the direction of the
# STable direction

def correction_NN(q):

    # Scale all the features
    scaler = joblib.load('scaler.gz')
    scaler_target = joblib.load('scaler_target.gz')

    # Reshape q with 1 row and 6 columns
    q_reshaped = q.reshape(1,-1)
    #print("q_reshaped: ", q_reshaped)

    q_scaled = scaler.transform(q_reshaped)
    #print("q_scaled: ", q_scaled)

    # Load back the NN model
    model = keras.models.load_model('best_model.keras')

    # Use model to predict dv
    dv_scaled = model.predict(q_scaled, verbose=0)

    dv = scaler_target.inverse_transform(dv_scaled)
    return dv
