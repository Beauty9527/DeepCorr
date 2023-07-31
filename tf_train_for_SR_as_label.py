import argparse
import os
import tensorflow as tf

from tensorflow.keras import models, layers
from tensorflow.keras.utils import plot_model
#from keras import models,layers

import data_generator
import numpy as np
import tensorflow.keras.datasets as Dataset


# tf.config.list_physical_devices('GPU')
# os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"  # 使用Gpu


def build_parser():
    """Setup option parsing for sample."""
    parser = argparse.ArgumentParser(description='Encoder the split mpileup files into split hdf file.')
    parser.add_argument('--hdf-files', type=str, help='The complete mpileup file.')
    parser.add_argument('--check-point-path', type=str, help='Path to save checkpoint file.')
    parser.add_argument('--model_hdf5_path', type=str, help='Path to save model.')


    args = parser.parse_args()

    return args

parsed_args = build_parser()


num_labels = 5
hdf_files = parsed_args.hdf_files
check_point_path = parsed_args.check_point_path
model_hdf5_path = parsed_args.model_hdf5_path
tensor_dtypes = [np.float32, np.int64]  # should be homologous with tensor_keys
tensor_keys = ["features", "labels"]  # what do you want from hdf file with a generator
input_feature_size = 11
chunk_len=100



model = models.Sequential([
    # layers.Bidirectional(layers.GRU(128,return_sequences=True),merge_mode='concat'),# 双向GRU实现
    layers.Bidirectional(layers.GRU(128, return_sequences=True,), input_shape=(1000,input_feature_size), name='Bi-GRU'),
    layers.Dense(units=num_labels, activation='softmax',name='Dense'),
])

model.build() # 还得看输入数据的格式 #input_shape=(None,1000,input_feature_size)
model.summary()
plot_model(model, to_file='./NFM_model.png', show_shapes=True)



model.compile(optimizer=tf.optimizers.Adam(1e-3),
              loss=tf.losses.sparse_categorical_crossentropy,
              metrics=['accuracy'],

)

# model.fit(train_data,
#           batch_size=128,
#           epochs=50,
#           validation_data=val_data,
#           callbacks=tf.keras.callbacks.ModelCheckpoint(filepath=check_point_path,
#                                                        save_best_only=True,
#                                                        save_weights_only=True,
#                                                        monitor='val_loss',
#                                                        verbose=1)
# )

iter_train_data = data_generator.load_hdf_iter_percent(hdf_files, keys=tensor_keys, mode="train")
iter_val_data = data_generator.load_hdf_iter_percent(hdf_files, keys=tensor_keys, mode="val")
iter_test_data = data_generator.load_hdf_iter_percent(hdf_files, keys=tensor_keys, mode="test")

callback_wrs = [tf.keras.callbacks.ModelCheckpoint(filepath=check_point_path,
                                                       save_best_only=True,
                                                       save_weights_only=True,
                                                       monitor='val_loss',
                                                       verbose=1,
                                                       ),
                tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=8, mode='min', min_delta=0.01)]

history = model.fit(iter_train_data,
          # steps_per_epoch=1,  # Sequence 实例的这个值好像不用指定
          epochs=2, #it is simple and data enough,so epochs 2 is good200
          validation_data=iter_val_data,
          callbacks=callback_wrs,
          )
scores = model.evaluate(iter_test_data,verbose=1)

import matplotlib.pyplot as plt
# 绘图函数
def print_history(history):
    # 绘制训练 & 验证的准确率值
    plt.plot(history.history['accuracy'])
    plt.plot(history.history['loss'])
    plt.title('Model accuracy&loss')
    plt.xlabel('Epoch')
    plt.legend(['Train_acc', 'Train_loss'])
    plt.show()


print(history.params)
# print_history(history)  # 调用绘图函数

#scores = model.evaluate(test_data, verbose=1)  # 这句虽然简单但是就这样用的

print(model.metrics_names)
print("the scores of model.evaluate(test_data):")
print('test loss', scores[0])
print('test accuracy', scores[1])

model.save(model_hdf5_path)
print("model_train success")




# plt.plot(history.history['val_accuracy'])
# plt.plot(history.history['val_loss'])
# plt.title('Model validation accuracy & loss')
# plt.xlabel('Epoch')
# plt.legend(['Val_acc', 'Val_loss'])
# #plt.show()
# # 保存图片到文件
# plt.savefig('/itmslhppc/itmsl0212/validation/deephc_sr/fruit_fly/ont/val_loss.png')
#
# # 关闭 Matplotlib 绘图
# plt.close()
# README
# We tested the model on the original train dataset,the acc and loss are the same as test data
# reasons come from three point: dataset too small,not repeat enough,model is too simple,encoder is too simple
# check the reason after finishing infer program
# 调试参数：--hdf-files ./covered_fasta.hdf --check-point-path ./weight/cp.ckpt --model_hdf5_path ./weight/model.h5