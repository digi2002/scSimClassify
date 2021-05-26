from keras.layers import Dense,Dropout
from keras.models import Sequential,Model,Input
from keras import initializers
from keras.initializers import Orthogonal
from collections import defaultdict
from keras.models import model_from_yaml
import string
import operator,os
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
from keras.layers import BatchNormalization
from keras.models import load_model
import random as rn
import numpy as np
from tensorflow import set_random_seed


def convertPrediction(labelset, predictions):
    label_pred = []
    for pred in predictions:
        scoredict = {}
        for idx, label in zip(range(len(labelset)), labelset):
            scoredict[label] = pred[idx]  # label, score
        inverse = [(value, key) for key, value in scoredict.items()]
        label_pred.append(max(inverse)[1])
    return np.asarray(label_pred)


def mlpmodel(modeldir,modelname,labelDict, layersizes, dropput,X_train, X_test, y_train, y_test, label_test):
    SEED = 0
    os.environ['PYTHONHASHSEED'] = str(SEED)
    np.random.seed(SEED)
    set_random_seed(SEED)
    rn.seed(SEED)

    initializer = Orthogonal()
    model = Sequential()
    if len(layersizes)==1:
        model.add(Dense(layersizes[0], kernel_initializer=initializer,activation='relu'))
        model.add(Dropout(dropput))

    elif len(layersizes)==2:
        model.add(Dense(layersizes[0], kernel_initializer=initializer,activation='relu'))
        model.add(Dropout(dropput))

        model.add(BatchNormalization())
        model.add(Dense(layersizes[1], kernel_initializer=initializer,activation='relu'))
        model.add(Dropout(dropput))

    elif len(layersizes)==3:
        model.add(Dense(layersizes[0], kernel_initializer=initializer,activation='relu'))
        model.add(Dropout(dropput))

        model.add(BatchNormalization())
        model.add(Dense(layersizes[1], kernel_initializer=initializer,activation='relu'))
        model.add(Dropout(dropput))

        model.add(BatchNormalization())
        model.add(Dense(layersizes[2],  kernel_initializer=initializer,activation='relu'))
        model.add(Dropout(dropput))

    model.add(BatchNormalization())
    model.add(Dense(len(set(labelDict.values())), kernel_initializer=initializer,activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

    model.fit(X_train, y_train, batch_size=len(X_train), epochs=200, verbose=2, validation_data=(X_test, y_test))
    print(model.summary())
    model.save(modeldir + modelname)

    predictions = model.predict(X_test)
    label_pred = convertPrediction(set(labelDict.values()), predictions)
    return label_pred
