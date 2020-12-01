import matplotlib.pyplot as plt


dt = [0.5,0.4,0.3,0.2,0.1,0.05,0.01,0.001,0.0001]


# N2 data
N2vv = [None, None, 2184.3, 2233.9, 2162.1, 2145.8, 2139.3, 2139.5, 2139.6]
N2se = [None, None, 2184.3, 2195.6, 2148.5, 2142.6, 2139.3, 2139.5, 2139.6]
# relative error parameters
N2vvv0 = [None,None]+[abs(2350.4-i)/2350.4 for i in N2vv[2:]]
N2sev0 = [None,None]+[abs(2350.4-i)/2350.4 for i in N2se[2:]]


Evv = [None, 1.2982650789018275, 0.0038484992841378318, 0.0009446456263159236, 0.00019612604053244364,
       4.7112971985748675e-05, 1.8610224055486248e-06, 1.860064927066817e-08, 1.8601434126077032e-10]
Ese = [2.66601794217714, 1.7522962468170533, 0.02741005676044304, 0.010372556036301897, 0.0036747331687280874,
       0.0016052025511824784, 0.00029163059613603677, 2.8572656043010745e-05, 2.8514994600568654e-06]



# O2 data
O2vv = [None, None, 1499.0, 1456.2, 1440.2, 1432.3, 1433.1,1433.0,1433.0]
O2se = [None, None, 1456.2, 1440.2, 1432.3, 1432.3, 1432.3,1432.9,1432.9]
O2vvv0 = [None, None]+[abs(1433.0-i)/1433.0 for i in O2vv[2:]]
O2sev0 = [None, None]+[abs(1433.0-i)/1433.0 for i in O2se[2:]]

EvvO2 = [2.2654726180455977, 2.100498344515442, 0.050630247617975535, 0.01872956166570478, 0.004594774839335882,
         0.0011452031703401272, 4.576208032952731e-05, 4.576017935344402e-07, 4.576125199747562e-09]
EseO2 = [3.7580597861718332, 3.901277037285156, 0.28814762855545756, 0.1380286836786054, 0.05768889914821934,
         0.02685441075051831, 0.005104092675932562, 0.0005049265477007818, 5.043872664446501e-05]

#------------------------------------------------------------------------
# Frequency - time step

'''
# Plot Frequency-dt
plt.title('Nitrogen: Vibrational Frequency vs Time Step')
plt.xlabel('time step [Å(u/eV)^0.5]')
plt.ylabel('Frequency [cm^-1]')
plt.xscale('log')
plt.plot(dt, N2vv, 'x',color='#de5855', ls='-',label='Velocity Verlet')
plt.plot(dt, N2se, 'x',color='#546bd1', ls='--',label='Symplectic Euler')
plt.legend()
plt.show()
# '∆ν/ν0'
'''


# pure vs roto stored in special files **
t = []
pv = []
rv = []
with open('separation_O2_r.txt','r') as rv_file:
  lines = rv_file.readlines()
  for line in lines:
    time, v = line.split(' ')
    t.append(float(time))
    rv.append(float(v))
with open('separation_O2_p.txt','r') as pv_file:
  lines = pv_file.readlines()
  for line in lines:
    pv.append(float(line.split(' ')[1]))



# pure vibration vs roto_vibration
plt.title('Oxygen Separation vs Time')
plt.xlabel('time [Å(u/eV)^0.5]')
plt.ylabel('Separation [Å]')
plt.plot(t, pv, color='#de5855', label='Pure Vibration')
plt.plot(t, rv, color='#546bd1', label='Roto Vibration')
plt.legend()
plt.show()

#----------------------------------------------------------------
# Max energy error vs timestep
'''
# Plot max(∆E/E0)-dt
plt.title('Oxygen: Maximum Relative Energy Error vs Time Step')
plt.xlabel('time step [Å(u/eV)^0.5]')
plt.ylabel('max(∆E/E0)')
plt.xscale('log')
plt.plot(dt, EvvO2, 'x', color='#de5855', ls='-', label='Velocity Verlet')
plt.plot(dt, EseO2, 'x', color='#546bd1', ls='--', label='Symplectic Euler')
plt.legend()
plt.show()
'''

# centrifugal -> larger avg sep -> lower F -> larger T -> lower cm-1 
