import numpy as np
import random
import matplotlib.pyplot as plt
import tensorflow as tf
import os

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense


# Target function f(x) = m*x + b
m_true = 2
b_true = 1

# Set seed for reproducibility
seed = 0
random.seed(seed)
np.random.seed(seed)
tf.random.set_seed(seed)

# Dati di validazione fissi per tutti gli esperimenti
x_valid = np.random.uniform(-1, 1, 100)
x_valid.sort()

# Funzione per creare, addestrare e plottare
def train_and_plot(ntrain, nepochs, sigma, filename_suffix):
    x_train = np.random.uniform(-1, 1, ntrain)
    y_target = m_true * x_valid + b_true
    y_train = np.random.normal(m_true * x_train + b_true, sigma)
    y_valid = np.random.normal(m_true * x_valid + b_true, sigma)

    model = Sequential()
    model.add(Dense(1,input_shape=(1,)))#essendo la funzione lineare, bastano 2 neuroni, uno di input e uno di output. il neurone di output 
    #ha attivazione lineare, cioè come la funzione stessa
    
    # Compilazione del modello: scelgo ottimizzatore, funzione di loss e metriche
    model.compile(optimizer='sgd', loss='mse' 
              #,metrics=['mse']
              )
    #loss è la funzione da minimizzare durante l'addestramento. Metrics invece 
    #sono valori calcolati per valutare le prestazioni del modello, ma non influenzano direttamente l’addestramento.

    #addestramento del modello
    history = model.fit(x_train, y_train,
                        batch_size=32,
                        epochs=nepochs,
                        shuffle=True, #mischio i dati ad ogni epoca
                        validation_data=(x_valid, y_valid),
                        verbose=0) #verbose = 0 non stampa a schero i dati dell'addestramento

    weights, bias = model.get_weights()
    m_fit = weights[0][0] #peso
    b_fit = bias[0] #intercetta/bias

    #predizione della rete neurale
    y_pred = model.predict(x_valid)

    # Grafico dei dati di validazione, della funzione target e della funzione predetta
    plt.figure()
    plt.plot(x_valid, y_target, label='Target (m=2, b=1)')
    plt.scatter(x_valid, y_valid, color='r', label='Validation data', s=10)
    plt.plot(x_valid, y_pred, label=f'Predicted (m={m_fit:.2f}, b={b_fit:.2f})')
    plt.legend()
    plt.grid(True)
    plt.title(f"ntrain={ntrain}, sigma={sigma}, epochs={nepochs}")
    plt.savefig(f"grafici/funzione_{filename_suffix}.png")
    plt.close()

    # Grafico della loss per i dati di validazione e per i dati di allenamento
    plt.figure()
    plt.plot(history.history['loss'], label='Training Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    plt.title(f"Loss - ntrain={ntrain}, sigma={sigma}, epochs={nepochs}")
    plt.xlabel("Epoch")
    plt.ylabel("MSE Loss")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"grafici/loss_{filename_suffix}.png")
    plt.close()


# ciclo 1: vario n_epochs con sigma=0, ntrain=1000 
for nepochs in [5, 30, 100]:
    train_and_plot(ntrain=1000, nepochs=nepochs, sigma=0, filename_suffix=f"epochs{nepochs}")

# ciclo 2: vario n_train con sigma=0, epochs=30 
for ntrain in [10, 1000, 10000]:
    train_and_plot(ntrain=ntrain, nepochs=30, sigma=0, filename_suffix=f"ntrain{ntrain}")

# ciclo 3: vario sigma con ntrain=1000, epochs=30 
for sigma in [0.0, 0.2, 0.6]:
    train_and_plot(ntrain=1000, nepochs=30, sigma=sigma, filename_suffix=f"sigma{int(sigma*10)}")