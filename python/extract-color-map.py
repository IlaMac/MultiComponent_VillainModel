import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

cmap = plt.get_cmap("hsv")
# rgb_cm = cmap.colors

arr = []
for i in range(0, 256):
  plt.plot(i, 0, ".", color=cmap(i))

  arr += [cmap(i)[0:3]]


print(arr)

# plt.show()