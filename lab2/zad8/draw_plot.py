import matplotlib.pyplot as plt
import numpy as np

DATA_TYPE = "double"
DATA_FILE = "wyniki_" + DATA_TYPE + ".txt" 

def load_data(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            values = list(map(float, line.split()))
            data.append(values)
    return np.array(data)

data = load_data(DATA_FILE)

logh = data[:,0]
ef2 = data[:,1]
ef3 = data[:,2]
ec2 = data[:,3]
eb2 = data[:,4]
eb3 = data[:,5]

plt.figure(figsize=(10, 7))

plt.plot(logh, ef2, label="forward 2")
plt.plot(logh, ef3, label="forward 3")
plt.plot(logh, ec2, label="central 2")
plt.plot(logh, eb2, label="backward 2")
plt.plot(logh, eb3, label="backward 3")

plt.xlabel("log10(h)")
plt.ylabel("log10(|error|)")
plt.title(f"Błędy różnicowe dla typu {DATA_TYPE}")
plt.grid(True)
plt.legend()
plt.show()