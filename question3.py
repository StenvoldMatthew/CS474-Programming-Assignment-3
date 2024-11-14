import numpy as np
import argparse
import matplotlib.pyplot as plt
from PIL import Image

def showImages(images, titles):
    # Determine the number of images
    num_images = len(images)
    
    # Calculate number of rows and columns
    if num_images > 4:
        num_rows = (num_images // 4) + (num_images % 4 > 0)  # Add an extra row if there are leftovers
        num_cols = 4
    else:
        num_rows = 1
        num_cols = num_images

    # Create subplots with the calculated rows and columns
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, num_rows * 5))
    axes = axes.flatten()  # Flatten the axes array for easy indexing

    # Plot each image
    for ax, img, title in zip(axes, images, titles):
        ax.imshow(img, cmap='gray', vmin=0, vmax=255)
        ax.set_title(title)
        ax.axis('off')

    # Hide any unused subplots
    for i in range(num_images, len(axes)):
        axes[i].axis('off')

    plt.tight_layout()  # Adjust layout to prevent overlap
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
  if isign == -1:
    data = convertToInput2d(data, data.shape[0], data.shape[1])
  else:
    data /= (data.shape[0]//2)*(data.shape[1]//2)
  
  # Each pixel is represented by the following matrix: | R  I |
  # The I values are set to the most recently          | I  0 |
  # updated value before computing

  # Apply 1D FFT on rows
  for i in range(1, data.shape[0], 2):
    data[i, 2::2] = data[i+1, 1::2] 
    fft(data[i, :], data.shape[1]//2, isign)

  # Apply 1D FFT on columns
  for j in range(1, data.shape[1], 2):
    data[2::2, j] = data[1::2, j+1]
    fft(data[:, j], data.shape[0]//2, isign)

  if isign == 1:
    return data[1::2, 1::2]
  return data

def convertToInput2d(arr, N, M):
  data = np.zeros((N*2+1, M*2+1))
  for i in range(N):
    for j in range(M):
      data[2 * i + 1, 2* j + 1] = arr[i, j]  # Real
  return data

def extractResults(data):
  real = data[1::2, 1::2]
  imag = data[2::2, 1::2]
  mag = np.zeros_like(real)
  for i in range(real.shape[0]):
    for j in range(real.shape[1]):   
      mag[i,j] = np.sqrt(real[i,j]**2 + imag[i,j]**2)
  return real, imag, mag

def runTest(filename):
  images = []
  titles = []
  image = np.array(Image.open(filename).convert('L'), dtype=np.float32)
  images.append(image)
  titles.append("Starting Image")

  fft = fft2d(image, -1)
  real, imag, mag = extractResults(fft)

  #3a 
  zeroPhase = np.zeros_like(fft)
  zeroPhase[1::2, 1::2] = mag

  inverseZP = fft2d(zeroPhase, 1)
  images.append(inverseZP)
  titles.append("Magnitude Only")

  #3b
  phaseOnly = np.zeros_like(fft)
  for i in range(real.shape[0]):
    for j in range(real.shape[1]):
      phaseOnly[2 * i + 1, 2* j + 1] = np.cos(np.arctan2(imag[i, j], real[i, j]))
      phaseOnly[2 * i + 2, 2* j + 2] = np.sin(np.arctan2(imag[i, j], real[i, j]))
  inversePO = fft2d(phaseOnly, 1)
  normalizedPO = 255 * (inversePO - inversePO.min()) / (inversePO.max() - inversePO.min())
  images.append(normalizedPO)
  titles.append("Phase Only")

  showImages(images, titles)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Image Converter')
    parser.add_argument('-f','--filename', type=str, default = "lenna.png", help='Which question you want run')
    args = parser.parse_args()
    runTest(args.filename)
