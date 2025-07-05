import tensorflow as tf
from tensorflow import keras
import numpy as np
import matplotlib.pyplot as plt
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from keras.layers import Flatten, Conv2D, MaxPooling2D
from PIL import Image, ImageEnhance, ImageOps #librerie per rendere le immagini più simili possibili a quelle di training


#imposto i seed per la replicabilità
seed=0
np.random.seed(seed) 
tf.random.set_seed(seed)

#dimensioni dell'immagine di input
img_rows, img_cols = 28, 28 # numero di pixels 
# output
num_classes = 10 # 10 digits

#dati di validazione e allenamento
(X_train, Y_train), (X_test, Y_test) = mnist.load_data()

#cambio la forma dei dati di input da (60000,28,28) a (60000,28,28,1), CNN vuole infatti il numero di canali, in questo caso 1(l'immagine è in scala di grigi)
if keras.backend.image_data_format() == 'channels_first':
    X_train = X_train.reshape(X_train.shape[0], 1, img_rows, img_cols)
    X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
    input_shape = (1, img_rows, img_cols)
else:
    X_train = X_train.reshape(X_train.shape[0], img_rows, img_cols, 1)
    X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)
    input_shape = (img_rows, img_cols, 1)

#trasforma i numeri interi delle etichette in vettori binari (one-hot)
Y_train = keras.utils.to_categorical(Y_train, num_classes)
Y_test = keras.utils.to_categorical(Y_test, num_classes)

def create_CNN():
    # inizializzo il modello
    model = Sequential()
    # primo convutional layer con 10 filtri 5X5
    model.add(Conv2D(20, kernel_size=(5, 5),
                     activation='relu',
                     input_shape=input_shape))
    #primo output 24X24X10
    model.add(MaxPooling2D(pool_size=(2, 2)))#riduce ogni dimensione spaziale di 2
    #secondo output 12X12X10
    model.add(Dropout(0.25))

    #secondo convutional layer con 20 filtri 5X5 
    model.add(Conv2D(35, kernel_size=(5, 5),
                     activation='relu'))
    #terzo output 8X8X20  
    model.add(MaxPooling2D(pool_size=(2, 2)))#riduce ogni dimensione spaziale di 2
    #quarto output 4X4X20
    model.add(Dropout(0.25))

    #aggiunta layer DNN
    model.add(Flatten())  # da (4,4,20) a vettore di dimensione 4*4*20=320

    #primo layer DNN
    model.add(Dense(4*4*20, activation='relu'))
    # apply dropout with rate 0.5
    model.add(Dropout(0.5))
    #secondo layer DNN
    model.add(Dense(160, activation='relu'))
    #output layer con softmax. Softmax trasforma i valori grezzi in probabilità
    model.add(Dense(num_classes, activation='softmax'))

    return model
 
    
#compilazione del modello
def compile_model():
    model=create_CNN()
    model.compile(loss=keras.losses.categorical_crossentropy,
                  optimizer='SGD',
                  metrics=['acc'])
    return model

# parametri di addestramento
batch_size = 64
epochs = 10

# creo la rete neurale 
model_CNN = compile_model()

# addestramento dei DNN
history = model_CNN.fit(X_train, Y_train,
              batch_size=batch_size,
              epochs=epochs,
              verbose=0,
              validation_data=(X_test, Y_test))

#valutazione della rete
score = model_CNN.evaluate(X_test, Y_test, verbose=0)
print(f"Test accuracy: {score[1]*100:.2f}%")

#carico le immagini e valuto la predizione della rete neurale su di esse
for i in range(0,9):
    digit_filename = f"pictures/{i}.png"
    
    # Carico immagine in scala di grigi
    digit_in = Image.open(digit_filename).convert('L')

    # Resize a 28x28
    digit_in = digit_in.resize((28, 28), Image.LANCZOS)

    # Converto in array e normalizza tra 0 e 1
    data = np.array(digit_in).astype(np.float32) / 255.0

    # Inverto colori solo se serve (se sfondo è chiaro)
    if np.mean(data) > 0.5:
        data = 1.0 - data

    # Reshape per input nel modello
    input_data = data.reshape(1, 28, 28, 1)

    # Predizione
    prediction = model_CNN.predict(input_data)
    predicted_class = np.argmax(prediction)
    confidence = np.max(prediction)
    
    # Visualizzo immagine con digit predetto sopra e confidenza
    plt.imshow(data, cmap='gray')
    plt.title(f'Predicted digit: {predicted_class}. Confidence: {confidence}', fontsize=16, color='red')
    plt.axis('off')
    
    # Salva immagine con predizione sovrapposta
    plt.savefig(f"grafici/{i}_pred.png")

    plt.close()  # Chiude la figura per liberare memoria