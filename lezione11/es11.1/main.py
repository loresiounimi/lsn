import numpy as np
import random
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation


# Parametri target della funzione f(x) = m*x + b
m = 2  # coefficiente angolare (pendenza)
b = 1  # intercetta

# imposto i seed per replicabilità dei dati
seed = 0
random.seed(seed)
np.random.seed(seed)
tf.random.set_seed(seed)

#genero i dati di input
x_train = np.random.uniform(-1, 1, 1000)
x_valid = np.random.uniform(-1, 1, 100)
x_valid.sort()
y_target = m * x_valid + b # target ideale, funzione lineare

sigma = 0.0 # rumore sulle y
y_train = np.random.normal(m * x_train + b, sigma) # valori su cui eseguire l'addestramento
y_valid = np.random.normal(m * x_valid + b, sigma) #valori su cui eseguire la validazione

#creazione del modello
model = tf.keras.Sequential()
model = tf.keras.Sequential()
model.add(Dense(15, input_shape=(1,), activation='relu'))  # hidden layer 1
model.add(Dense(18, activation='relu'))                    # hidden layer 2
model.add(Dense(12, activation='relu'))                    # hidden layer 3
model.add(Dense(1))  # output layer (con attivazione lineare di default)


# Compilazione del modello: scelgo ottimizzatore, funzione di loss e metriche
model.compile(optimizer='sgd', loss='mse' 
              #,metrics=['mse']
              )
#loss è la funzione da minimizzare durante l'addestramento. Metrics invece 
#sono valori calcolati per valutare le prestazioni del modello, ma non influenzano direttamente l’addestramento.

#addestramento del modello
history = model.fit(x=x_train, y=y_train, 
          batch_size=32, epochs=5,
          shuffle=True, # mischio i dati di input ad ogni epoche
          validation_data=(x_valid, y_valid),
          verbose=0)#verbose = 0 non stampa a schero i dati dell'addestramento

weights, bias = model.layers[-1].get_weights() #prendo peso e bias dell'ultimo neurone
m_fit = weights[0][0]  # peso
b_fit = bias[0]        # bias/intercetta

#ordinate predette dalla rete neurale
y_pred = model.predict(x_valid)

# Valutazione finale del modello sui dati di validazione
score = model.evaluate(x_valid, y_valid, batch_size=32, verbose=0)

#stampa della loss
print('Test loss on validation data (MSE):', score)

#print('Test metric :', score[1])

#plotto dei dati di validazione, della funzione target e del predict dela rete neurale
plt.plot(x_valid, y_target, label='target')
plt.scatter(x_valid, y_valid, color='r', label='validation data')
plt.plot(x_valid, y_pred, label=f'Predicted fit (m={m_fit:.2f}, b={b_fit:.2f})')
plt.legend()
plt.grid(True)
plt.savefig("grafici/complex_funzione.png")

plt.clf()

#plotto la loss sui dati di validazione ed allenamento
plt.plot(history.history['loss'])# Loss sul training set per ogni epoca
plt.plot(history.history['val_loss'])   # Loss sul validation set per ogni epoca
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='best')
plt.grid(True)
plt.savefig("grafici/complex_loss.png")

#aumentando sigma ma tenendo fissi Ndata e Nepochs peggiora la loss essendo i dati più "sporcati" 
#aumentando Nepochs ma tenendo fissi Ndata e sigma la loss migliora e tende asintoticamente a un certo valore. teoricamente se aumento le epoche
#dovrei osservare un peggioramento della loss perchè rischio l'overfitting sui data train. In questo caso per osservare overfitting devo aumentare
#i Ndata train sensibilmente, avere sigma>0 e aumentare la complessità del modello con layer a più neuroni. un layer con 1 solo neurone è più difficili
#che overfitti per la poca complessità e anche perchè la funzione target è allineata con l'output del neurone stesso.
#se aumento Ndata ma tengo fissi Nepochs e sigma la loss migliora


