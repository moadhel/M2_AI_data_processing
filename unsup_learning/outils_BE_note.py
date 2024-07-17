import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import os
os.environ["KERAS_BACKEND"] = "torch"
from keras_core.datasets import mnist

# Fonctions utiles ------------------------------------------------------------------------------------------------

def showcentroids(M):
    """
    Affiche les centroides (lignes) dans la matrice M de taille k x 28**2, avec
    k un multiple de 10.
    """
    n = M.shape[0]
    fig, axarr=plt.subplots(nrows=int(n/10), ncols=10, figsize=(16, 10))
    axarr=axarr.flatten()
    
    for i in range(n):        
        axarr[i].imshow(np.reshape(M[i],(28,28)), cmap='gray')
        axarr[i].set_xticks([])
        axarr[i].set_yticks([])
        axarr[i].set_title("i=: {}".format(i))
    plt.show()

def showconf(labels, estimated_labels):
    """
    Affiche la matrice de confusion avec les erreurs de classification par classe.
    Entrées :
    - labels : vecteur contenant les vrais labels des données
    - estimated_labels : vecteur de même taille que labels, contenant les labels estimés 
    """    
    k0 = 10
    cm = confusion_matrix(labels, estimated_labels, labels=range(k0))
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=range(k0))
    disp.plot()


def loadmnist():
    """
    Charge les données d'entraînement et de test de MNIST
    """    
    (x_train, y_train), (x_test, y_test) = mnist.load_data()

    x_train = x_train.astype('float32') 
    x_test = x_test.astype('float32')

    # Normalisation
    x_train = x_train/255.0
    x_test = x_test/255.0

    # Vectorisation
    X_train = x_train.reshape(len(x_train),-1)
    X_test = x_test.reshape(len(x_test),-1)

    return X_train, y_train, X_test, y_test

def get_map_labels(estimated_labels, true_labels):
    """
    Renvoie un vecteur qui permet d'associer à chaque label produit par 
    l'algorithme de clustering l'un des vrais labels des données.
    Entrées :
      - estimated_labels : vecteur de labels estimés produits par l'algorithme
      - true_labels : vecteur de vrais labels pour les mêmes données que estimated_labels
    Sortie :
      - vecteur de taille k (nombre de labels produits par l'algorithme) avec les
        associations
    """        
    # Nombre de labels produits par l'algorithme
    k = int(estimated_labels.max()) + 1
    # Nombre de vrais labels 
    k0 = int(true_labels.max()) + 1
    map_labels = np.zeros(k)
    # On choisit un vrai label pour associer à chaque label estimé
    for i in range(k):
        index = np.where(estimated_labels == i,1,0)
        true_label = np.bincount(true_labels[index==1], None, k0).argmax()
        map_labels[i] = true_label
    return map_labels

