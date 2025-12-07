import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

df = pd.read_csv('wynik/fehlberg_lorenz.txt', sep=';')
# df = pd.read_csv('wynik/rk4_lorenz.txt', sep=';')
# df = pd.read_csv('wynik/fehlberg_lorenz.txt', sep=';')

df.columns = df.columns.str.strip()

y0 = df['y0']
y1 = df['y1']
y2 = df['y2']

fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')

ax.plot(y0, y1, y2, color='black', lw=1.5)
ax.set_xlabel('y0')
ax.set_ylabel('y1')
ax.set_zlabel('y2')

plt.tight_layout()
plt.show()
