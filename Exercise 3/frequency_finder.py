import numpy as np
from scipy.signal import find_peaks


def main():
  with open('separation.txt','r') as data_table:
    times,data = [],[]
    lines = data_table.readlines()
    for line in lines:
      split = line.split(' ')
      times.append(float(split[0]))
      data.append(float(split[1]))
    times = np.array(times)
    indices = find_peaks(data)[0]
    peak_times = times[indices]
    # period in units of amu*(Å/eV)^0.5
    period = (peak_times[-1]-peak_times[0])/(len(peak_times)-1)
  # Constants to calculate Period in seconds
  c = 299_792_458
  f = c/0.01
  T = 1/f
  L = 1.000_014_95*10**-10 
  m = 1.660_539_066_60*10**-27
  E = 1.602_176_634*10**-19
  # Calculate units that we have in seconds
  t = L*(m/E)**0.5

  # Chaning units of period (from ~10.18fs) into cm^-1
  period = T/(period*t)
  return period


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()


