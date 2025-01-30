## \file correction_module.py
#
# \brief Given IC, find velocity correction according to trained predictor.

import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PolynomialFeatures
import joblib   # To save/load scaler
import tensorflow as tf
import keras


DIM = 6     # dimension of RTBP system

## Given IC, find correction maneuver according to polynomial regression.
#
# @param[in]  q 	init. cond. in position-velocity (x,y,z,dx,dy,dz)
# @param[out] dv 	modulus of correction maneuver (in the direction of the
# STable direction

def correction_regression(q):

    # Scale all the features
    scaler = joblib.load('scaler_supervised.gz')
    scaler_target = joblib.load('scaler_target_supervised.gz')

    # Reshape q with 1 row and 6 columns
    q_reshaped = q.reshape(1,-1)
    #print("q_reshaped: ", q_reshaped)

    q_scaled = scaler.transform(q_reshaped)
    #print("q_scaled: ", q_scaled)

    # Load back the regression model
    # Here you can replace pickle with joblib or cloudpickle
    from pickle import load
    with open("model_supervised.pkl", "rb") as f:
        lin_reg = load(f)

    # Use model to predict dv
    poly_features = PolynomialFeatures(degree=2, include_bias=False)
    q_poly = poly_features.fit_transform(q_scaled)
    dv_scaled = lin_reg.predict(q_poly)

    dv = scaler_target.inverse_transform(dv_scaled)
    return dv


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
