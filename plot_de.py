import matplotlib.pyplot as plt
import pandas as pd
df = pd.read_csv("/nobackup/rmvn14/ham-chalk/output/timesteps.csv")
print(df)
plt.plot(df["time"],df["damage"]/df["damage"].max() ,label="Damage")
plt.plot(df["time"],df["energy"]/df["energy"].max() ,label="Energy")
plt.plot(df["time"],df["oobf"]/df["oobf"].max() ,label="OOBF")
plt.legend()
plt.show()
