from SimpleMD import SimpleMD


md = SimpleMD(100, 3, steps=150000, exeDir='thermo', force='harmonic', thermostat='white', temperature=1)
md.setup_positions(density=0.7)
md.setup_masses(1)
md.log_positions(period=5)
md.log_output()
md.execute()


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print(md.temperature)
plt.plot(md.out_times, md.temperature)
plt.axhline(y=np.mean(md.temperature))
plt.savefig('temperature.png')


md.clean_files()
