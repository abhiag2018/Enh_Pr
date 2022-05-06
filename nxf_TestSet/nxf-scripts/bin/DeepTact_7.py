#!/usr/bin/env python
#keras version: keras-1.2.0

import sys
from pathlib import Path
import os, re
import glob
import time
import random
import itertools
import datetime
import string
import numpy as np
import hickle as hkl
from sklearn import metrics
import pandas as pd
import h5py
from natsort import natsorted

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Convolution2D, MaxPooling2D, Flatten
from tensorflow.keras.layers import LSTM, Bidirectional
from tensorflow.keras.layers import Dense, Activation

from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import BatchNormalization
# from tensorflow.keras.layers import Reshape, Merge, Permute
from tensorflow.keras.layers import Reshape, Concatenate, Permute
from tensorflow.keras import optimizers
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
import tensorflow.keras.backend as K
from tensorflow.keras.models import load_model
# from tensorflow.keras.engine.topology import Layer, InputSpec
from tensorflow.keras.layers import Layer
# from tensorflow.keras import initializations

import graphtools as gt
"""
DeepTACT.py

@author: abhiag
"""

########################### Training ##########################
# Attention GRU network
class AttLayer(Layer):
    def __init__(self, **kwargs):
        # self.init = initializations.get('normal')
        #self.input_spec = [InputSpec(ndim=3)]
        super(AttLayer, self).__init__(**kwargs)

    def build(self, input_shape):
        assert len(input_shape)==3
        self.W = self.add_weight(shape=(input_shape[-1],), initializer="random_normal", trainable=True)
        super(AttLayer, self).build(input_shape)  # be sure you call this somewhere!

    def call(self, x):
        M = K.tanh(x)
        alpha = K.dot(M,K.expand_dims(self.W, axis=-1))
        alpha = K.squeeze(alpha,axis=-1)#.dimshuffle(0,2,1)

        ai = K.exp(alpha)
        weights = ai/K.expand_dims(K.sum(ai, axis=1),axis=-1)
        weighted_input = x*K.expand_dims(weights,axis=-1)
        return K.tanh(K.sum(weighted_input,axis=1))

    def get_output_shape_for(self, input_shape):
        return (input_shape[0], input_shape[-1])

def model_def():
    inp_region1_seq = Input(shape=(1, NUM_SEQ, RESIZED_LEN),name='enh_seq')
    inp_region2_seq = Input(shape=(1, NUM_SEQ, 1000),name='pr_seq')
    inp_region1_expr = Input(shape=(1, NUM_REP, RESIZED_LEN),name='enh_dnase')
    inp_region2_expr = Input(shape=(1, NUM_REP, 1000),name='pr_dnase')
    drop_rate = 0.5 
    conv_enhancer_seq = Sequential()
    conv_enhancer_seq.add(Convolution2D(1024, (NUM_SEQ, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_SEQ, RESIZED_LEN)))
    conv_enhancer_seq.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_enhancer_seq.add(Reshape((1024, (RESIZED_LEN-40+1)//20)))
    out_enh_seq = conv_enhancer_seq(inp_region1_seq)

    conv_promoter_seq = Sequential()
    conv_promoter_seq.add(Convolution2D(1024, (NUM_SEQ, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_SEQ, 1000)))
    conv_promoter_seq.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_promoter_seq.add(Reshape((1024, 48)))
    out_pr_seq = conv_promoter_seq(inp_region2_seq)

    merged_seq = Concatenate()([out_enh_seq, out_pr_seq])

    conv_enhancer_DNase = Sequential()
    conv_enhancer_DNase.add(Convolution2D(1024, (NUM_REP, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_REP, RESIZED_LEN)))
    conv_enhancer_DNase.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_enhancer_DNase.add(Reshape((1024, (RESIZED_LEN-40+1)//20)))
    out_enh_DNase = conv_enhancer_DNase(inp_region1_expr)

    conv_promoter_DNase = Sequential()
    conv_promoter_DNase.add(Convolution2D(1024, (NUM_REP, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_REP, 1000)))
    conv_promoter_DNase.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_promoter_DNase.add(Reshape((1024, 48)))
    out_pr_DNase = conv_promoter_DNase(inp_region2_expr)

    merged_DNase = Concatenate()([out_enh_DNase, out_pr_DNase])
    #   
    merged = Concatenate(axis=-2)([merged_seq, merged_DNase]) 

    merged = Permute((2, 1))(merged)
    merged = BatchNormalization()(merged)
    merged = Dropout(drop_rate)(merged)
    merged = Bidirectional(LSTM(100, return_sequences=True), merge_mode="concat")(merged)
    merged = AttLayer()(merged)
    merged = BatchNormalization()(merged)
    merged = Dropout(drop_rate)(merged)
    merged = Dense(925)(merged)
    merged = BatchNormalization()(merged)
    merged = Activation('relu')(merged)
    merged = Dropout(drop_rate)(merged)
    merged = Dense(1, activation = 'sigmoid')(merged)

    model = Model(inputs=[inp_region1_seq, inp_region2_seq, inp_region1_expr, inp_region2_expr], outputs=merged)
    return model

def f1(y_true, y_pred):
    cast = lambda x:tf.keras.backend.cast(x,dtype='float64')
    TP = K.sum(cast(K.equal(y_true, 1) & K.equal(K.round(y_pred), 1)))
    FP = K.sum(cast(K.equal(y_true, 0) & K.equal(K.round(y_pred), 1)))
    FN = K.sum(cast(K.equal(y_true, 1) & K.equal(K.round(y_pred), 0)))
    TN = K.sum(cast(K.equal(y_true, 0) & K.equal(K.round(y_pred), 0)))
    P = TP / (TP + FP + K.epsilon())
    R = TP / (TP + FN + K.epsilon())
    F1 = 2 * P * R / (P + R + K.epsilon())
    return F1

def data_gen(path, randomize=False):
    with h5py.File(path,'r') as h5f:
        # for i in np.random.permutation(h5f['label'].shape[0]):
        for i in range(h5f['label'].shape[0]):
            yield {'enh_dnase': h5f['enh_dnase'][i], 'pr_dnase': h5f['pr_dnase'][i], 'enh_seq': h5f['enh_seq'][i], 'pr_seq': h5f['pr_seq'][i]}, h5f['label'][i]

def data_gen_1(path, randomize=True):
    for dataset in natsorted(glob.glob(path)):

        DF = pd.read_csv(dataset)


        # print('\n','/'.join(path.split(os.path.sep)[-2:]))
        DF[DF.columns.drop('label')] = DF[DF.columns.drop('label')].applymap(eval).applymap(list)

        DF['enh_dnase'] = DF.apply(lambda row: np.array([row[f'enh_dnase_{x}'] for x in range(NUM_REP)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)
        DF['pr_dnase'] = DF.apply(lambda row: np.array([row[f'pr_dnase_{x}'] for x in range(NUM_REP)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)
        DF['enh_seq'] = DF.apply(lambda row: np.array([row[f'enh_seq_{x}'] for x in range(NUM_SEQ)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)
        DF['pr_seq'] = DF.apply(lambda row: np.array([row[f'pr_seq_{x}'] for x in range(NUM_SEQ)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)

        DF = DF[['enh_dnase', 'pr_dnase', 'enh_seq', 'pr_seq','label']]
        if randomize:
            DF = DF.sample(frac=1, axis=0).reset_index(drop=True)

        DF[['enh_dnase', 'pr_dnase', 'enh_seq', 'pr_seq']] = DF[['enh_dnase', 'pr_dnase', 'enh_seq', 'pr_seq']].applymap(tf.convert_to_tensor)
        for index, row in DF.iterrows():
            label = row.pop('label')
            yield dict(row),label

def train(cell_data, train_csv, val_csv, tensorboard_logdir):
    ### CHECK GPU usage
    print(tf.config.list_physical_devices('GPU'),flush=True)
    print(tf.test.gpu_device_name(),flush=True)

    output_types = {'enh_seq':tf.float64, 'pr_seq':tf.float64, 'enh_dnase':tf.float64, 'pr_dnase':tf.float64}
    output_shapes = {'enh_seq':[1,NUM_SEQ,RESIZED_LEN], 'pr_seq':[1,NUM_SEQ,1000], 'enh_dnase':[1,NUM_REP,RESIZED_LEN], 'pr_dnase':[1,NUM_REP,1000]}

    shuffle_buffer_size = 256

    train_set = tf.data.Dataset.from_generator(lambda:data_gen(f"{cell_data}/{train_csv}"), output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
    train_set = train_set.shuffle(shuffle_buffer_size, reshuffle_each_iteration=False).batch(BATCH_SIZE)

    val_set = tf.data.Dataset.from_generator(lambda:data_gen(f"{cell_data}/{val_csv}"), output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
    val_set = val_set.shuffle(shuffle_buffer_size, reshuffle_each_iteration=False).batch(BATCH_SIZE)

    model = model_def()
    print('compiling...')
    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers.Adam(lr = 0.00001),
                  metrics = ['acc', f1])

    modelCheckpoint_1 = ModelCheckpoint(f"{tensorboard_logdir}/best_model_acc.h5", monitor = 'val_acc', save_best_only = True, mode = 'max')
    modelCheckpoint_2 = ModelCheckpoint(f"{tensorboard_logdir}/best_model_loss.h5", monitor = 'val_loss', save_best_only = True, mode = 'min')
    modelCheckpoint_3 = ModelCheckpoint(f"{tensorboard_logdir}/"+"weights.{epoch:02d}.h5")

    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=tensorboard_logdir, histogram_freq=1)

    print('fitting...')
    model.fit(train_set, validation_data=val_set, epochs = 40, callbacks = [modelCheckpoint_1, modelCheckpoint_2, modelCheckpoint_3,tensorboard_callback])


########################### Evaluation ##########################


def evaluate(model_paths, eval_data, limit_data=None, outfile="out.npz"):
    output_types = {'enh_seq':tf.float64, 'pr_seq':tf.float64, 'enh_dnase':tf.float64, 'pr_dnase':tf.float64}
    output_shapes = {'enh_seq':[1,NUM_SEQ,RESIZED_LEN], 'pr_seq':[1,NUM_SEQ,1000], 'enh_dnase':[1,NUM_REP,RESIZED_LEN], 'pr_dnase':[1,NUM_REP,1000]}


    model = model_def()
    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers.Adam(lr = 0.00001),
                  metrics = ['acc', f1])

    # test_labels = np.array(pd.concat([pd.read_csv(f).label for f in natsorted(glob.glob(eval_data))], axis = 0))
    with h5py.File(eval_data,'r') as h5f:
        test_labels = np.array(h5f['label'])
    avg_score = np.zeros((test_labels.shape[0],1))


    for model_p in model_paths:
        print("evaluating ", model_p, flush=True)
        print()
        print()
        model.load_weights(model_p)

        test_set = tf.data.Dataset.from_generator(lambda:data_gen(eval_data, randomize=False), output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
        test_set = test_set.batch(BATCH_SIZE)

        # from keras_tqdm import TQDMCallback
        score = model.predict(test_set, verbose=1)
        avg_score = avg_score + score

    avg_score = avg_score/len(model_paths)    

    np.savez( outfile, pred_score=avg_score, true_labels=test_labels)

    return


########################### Input #############################
if __name__=="__main__":
    job = sys.argv[1]
    if len(sys.argv)<2:
        print('[USAGE] python DeepTACT.py job <args>')
        print('For example, python DeepTACT.py split <data-dir>')
        sys.exit()
    
    RESIZED_LEN = 2000 # enhancer input: number of base pairs

    NUM_SEQ = 4 # number of base pairs DNA acids A,T,C,G = 4
    # NUM_ENSEMBL = int(sys.argv[5])

    BATCH_SIZE=32

    ############################ MAIN ###############################
    if job == "train":
        CELL = sys.argv[2]
        NUM_REP = int(sys.argv[3])
        time_append = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+'-'+''.join(random.choices(string.ascii_uppercase + string.digits,k=3))
        LOG_DIR = f"{CELL}/logs_{time_append}"
        os.makedirs(LOG_DIR, exist_ok=True)

        cell_data = f"{CELL}/data"
        train_csv = "dataset-train.h5" #training dataset
        val_csv = "dataset-val.h5" #val dataset

        train(cell_data, train_csv, val_csv, LOG_DIR)
    elif job == "test":
        eval_data = sys.argv[2] # {EVAL_CELL}/data/dataset-train.csv

        CELL = sys.argv[3]
        NUM_REP = int(sys.argv[4])

        append_str = sys.argv[5]
        cell_bootstrap_suffix = np.array([sys.argv[6:]]).flatten()

        model_paths = [CELL+f"/logs_{time_append}/best_model_loss.h5" for time_append in cell_bootstrap_suffix]


        outfile = f"eval_out_{append_str}.npz"


        evaluate(model_paths, eval_data, limit_data=None, outfile=outfile)

        print("output file:", outfile)














