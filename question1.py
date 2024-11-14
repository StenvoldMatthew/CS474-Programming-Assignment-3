import numpy as np
import argparse
import matplotlib.pyplot as plt

def display(data_series, titles, xlabel="Index", ylabel="Amplitude"):
  num_plots = len(data_series)
  plt.figure(figsize=(12, 2.5 * num_plots))
    
  for i, series in enumerate(data_series):
    plt.subplot(num_plots, 1, i + 1)
    plt.stem(series) 
    plt.title(titles[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
  plt.tight_layout()
  plt.show()


def fft(data: list[float], nn: int, isign: int):
  #n, mmax, m, j, istep, i = None # uint32
  #wtemp, wr, wpr, wpi, wi, theta = None # double
  #tempr, tempi = None # float

  n = nn*2
  j = 1
  for i in (range(1, n, 2)):
    if j > i:
      data[j], data[i] = data[i], data[j]
      data[j+1], data[i+1] = data[i+1], data[j+1]
    m = nn
    while (m >= 2 and j > m):
      j -= m
      m //= 2
    j += m

  mmax = 2
  while (n > mmax):
    istep = mmax *2
    theta = isign * (2*np.pi / mmax)
    wtemp = np.sin(.5*theta)
    wpr = -2.0 * wtemp * wtemp
    wpi = np.sin(theta)
    wr = 1.0
    wi = 0.0
    for m in range(1, mmax, 2):
      for i in range(m, n+1, istep):
        j = i + mmax
        tempr = wr*data[j]-wi*data[j+1]
        tempi = wr*data[j+1]+wi*data[j]
        data[j] = data[i]-tempr
        data[j+1] = data[i+1]-tempi
        data[i] += tempr
        data[i+1] += tempi

      wtemp = wr
      wr = wtemp*wpr-wi*wpi+wr
      wi = wi*wpr+wtemp*wpi+wi
    mmax=istep

def fft2d(data, isign):
  # Apply 1D FFT on rows
  for i in range(data.shape[0]):
    fft(data[i, :].tolist(), data.shape[1] // 2, isign)

  # Apply 1D FFT on columns
  for j in range(data.shape[1]):
    fft(data[:, j].tolist(), data.shape[0] // 2, isign)

  return data

def convertToInput(arr, nn):
  data = [0] * (2 * nn + 1)
  for i in range(nn):
    data[2 * i + 1] = arr[i]  # Real 
    data[2 * i + 2] = 0       # Imaginary
  return data

def extractResults(data, nn):
  real = [data[2 * i + 1] for i in range(nn)]
  imag = [data[2 * i + 2] for i in range(nn)]
  magnitude = [np.sqrt(real[i]**2 + imag[i]**2) for i in range(nn)]
  return real, imag, magnitude


def question1a():
  data_series = []
  titles = ["Real Part of DFT", "Imaginary Part of DFT", "Magnitude of DFT", "Inverse of DFT"]

  # FFT
  input = [2, 3, 4, 4]
  nn = 4
  data = convertToInput(input, nn)
  fft(data, nn, -1)
  data_series.extend(extractResults(data, nn))
  
  # Inverse FFT
  for i in range(1, len(data)):
    data[i] /= nn
  fft(data, nn, 1)
  inverse = [data[2 * i + 1] for i in range(nn)]
  data_series.append(inverse)
  display(data_series, titles)

def question1b():
  dataSeries = []
  titles = ["Cosine", "Shifted"]

  # Creating the Cosine
  N = 128  # Number of samples
  u = 8    # Frequency in cycles
  x = np.arange(N)  # Sample indices
  f_x = np.cos(2 * np.pi * u * x / N)
  dataSeries.append(f_x)

  # FFT
  data = convertToInput(f_x, N)
  for i in range(N):
    data[2 * i + 1] *= (-1) ** i
  dataSeries.append(data.copy())
  
  fft(data, N, -1)
  real, imag, mag = extractResults(data, N)
  phase = [np.arctan2(imag[i], real[i]) for i in range(N)]
  dataSeries.extend([real, imag, mag, phase])
  titles.extend(["Real", "Imaginary", "Magnitude", "Phase"])
  
  display(dataSeries, titles)

def question1c():
  dataSeries = []
  titles = ["Rect Function", "Shifted"]

  N = 128
  data = np.loadtxt("Rect_128.dat")
  dataSeries.append(data.copy())
  data = convertToInput(data, N)
  for i in range(N):
    data[2 * i + 1] *= (-1) ** i
  dataSeries.append(data.copy())
  
  fft(data, N, -1)
  real, imag, mag = extractResults(data, N)
  phase = [np.arctan2(imag[i], real[i]) for i in range(N)]
  dataSeries.extend([real, imag, mag, phase])
  titles.extend(["Real", "Imaginary", "Magnitude", "Phase"])
  
  display(dataSeries, titles)
  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Image Converter')
    parser.add_argument('-q','--question', type=str, default = "1a", help='Which question you want run')
    args = parser.parse_args()
    match args.question:
      case "1a":
        question1a()
      case "1b":
        question1b()
      case "1c":
        question1c()
      case _:
        print("Question not recognized")