import numpy as np
from scipy.fft import fft, fftfreq


E = np.array([
    [1.0, 2.0, 3.0, 4.0],   
    [2.0, 3.0, 4.0, 5.0],   
    [3.0, 4.0, 5.0, 6.0],   
    [4.0, 5.0, 6.0, 7.0],  
    [5.0, 6.0, 7.0, 8.0]])  


F1 = np.fft.fftn(E ,norm='ortho')
F2 = fft(E, axis=1, norm='ortho')

print(abs(F1 - F2))