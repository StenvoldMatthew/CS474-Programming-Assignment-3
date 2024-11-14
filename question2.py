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
  

  # Apply 1D FFT on rows
  for i in range(data.shape[0]):
    fft(data[i, :], data.shape[1]//2, isign)

  # Apply 1D FFT on columns
  for j in range(data.shape[1]):
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
  imag = data[2::2, 2::2]
  mag = np.zeros_like(real)
  for i in range(real.shape[0]):
    for j in range(real.shape[1]):   
      mag[i,j] = np.sqrt(real[i,j]**2 + imag[i,j]**2)
  return real, imag, mag

def checkFFT2d():
  images = []
  titles = []
  image = np.array(Image.open("boat.png").convert('L'), dtype=np.float32)
  images.append(image)
  titles.append("Starting Image")

  # Run fft and ifft
  fft = fft2d(image, -1)
  ifft = fft2d(fft, 1)

  images.append(ifft)
  titles.append("After FFT and IFFT")
  showImages(images, titles)


def question2a():
  images = []
  titles = []
  imageSize = 512
  squareSize = 32
  image = np.zeros((imageSize, imageSize), int)
  start = imageSize//2 - squareSize//2
  for i in range(squareSize):
    for j in range(squareSize):
      image[start+i, start+j] = 256
  images.append(image)
  titles.append("Starting Image")
  
  # Do FFT
  fft = fft2d(image, -1)
  real, imag, mag = extractResults(fft)
  images.append(mag)
  titles.append("Magnitude Pre Adjust")

  #Shift Magnitude
  rows, cols = mag.shape
  top_left = mag[:rows//2, :cols//2]
  top_right = mag[:rows//2, cols//2:]
  bottom_left = mag[rows//2:, :cols//2]
  bottom_right = mag[rows//2:, cols//2:]

  adjMag = np.zeros_like(mag)

  adjMag[:rows//2, :cols//2] = bottom_right
  adjMag[:rows//2, cols//2:] = bottom_left
  adjMag[rows//2:, :cols//2] = top_right
  adjMag[rows//2:, cols//2:] = top_left
  print(adjMag[:10, :10])
  adjMag = np.log(1 + adjMag)
  print(adjMag[:10, :10])


  images.append(adjMag)
  titles.append("Magnitude Post Adjust")

  showImages(images, titles)


# Part 2b

def question2b():
  images = []
  titles = []
  imageSize = 512
  squareSize = 64
  image = np.zeros((imageSize, imageSize), int)
  start = imageSize//2 - squareSize//2
  for i in range(squareSize):
    for j in range(squareSize):
      image[start+i, start+j] = 256
  images.append(image)
  titles.append("Starting Image")
  
  # Do FFT
  fft = fft2d(image, -1)
  real, imag, mag = extractResults(fft)
  images.append(mag)
  titles.append("Magnitude Pre Adjust")

  #Shift Magnitude
  rows, cols = mag.shape
  top_left = mag[:rows//2, :cols//2]
  top_right = mag[:rows//2, cols//2:]
  bottom_left = mag[rows//2:, :cols//2]
  bottom_right = mag[rows//2:, cols//2:]

  adjMag = np.zeros_like(mag)

  adjMag[:rows//2, :cols//2] = bottom_right
  adjMag[:rows//2, cols//2:] = bottom_left
  adjMag[rows//2:, :cols//2] = top_right
  adjMag[rows//2:, cols//2:] = top_left
  print(adjMag[:10, :10])
  adjMag = np.log(1 + adjMag)
  print(adjMag[:10, :10])


  images.append(adjMag)
  titles.append("Magnitude Post Adjust")

  showImages(images, titles)

# Part 2c

def question2c():
  images = []
  titles = []
  imageSize = 512
  squareSize = 128
  image = np.zeros((imageSize, imageSize), int)
  start = imageSize//2 - squareSize//2
  for i in range(squareSize):
    for j in range(squareSize):
      image[start+i, start+j] = 256
  images.append(image)
  titles.append("Starting Image")
  
  # Do FFT
  fft = fft2d(image, -1)
  real, imag, mag = extractResults(fft)
  images.append(mag)
  titles.append("Magnitude Pre Adjust")

  #Shift Magnitude
  rows, cols = mag.shape
  top_left = mag[:rows//2, :cols//2]
  top_right = mag[:rows//2, cols//2:]
  bottom_left = mag[rows//2:, :cols//2]
  bottom_right = mag[rows//2:, cols//2:]

  adjMag = np.zeros_like(mag)

  adjMag[:rows//2, :cols//2] = bottom_right
  adjMag[:rows//2, cols//2:] = bottom_left
  adjMag[rows//2:, :cols//2] = top_right
  adjMag[rows//2:, cols//2:] = top_left
  #adjMag = np.log(1 + adjMag)


  images.append(adjMag)
  titles.append("Magnitude Post Adjust")

  showImages(images, titles)
  


  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Image Converter')
    parser.add_argument('-q','--question', type=str, default = "check2d", help='Which question you want run')
    args = parser.parse_args()
    match args.question:
      case "2a":
        question2a()
      case "2b":
        question2b()
      case "2c":
        question2c()
      case "check2d":
        checkFFT2d()
      case _:
        print("Question not recognized")