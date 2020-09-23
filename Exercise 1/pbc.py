import math as m

def pbc(x,l):
  return [i%l for i in x]

def mic(x,l):
  x = pbc(x,l)
  corners = [[0,0,0],[0,0,l],[0,l,0],[l,0,0],[0,l,l],[l,0,l],[l,l,0],[l,l,l]]
  mic = l
  for corner in corners:
    mic = min(mic,m.sqrt(sum([(i-j)**2 for i,j in zip(x,corner)])))
  return mic

