import sys

t_ = [(1,5,11),(2,10, 15), (3,4,5)]
_t_ = [(124,15, 5), (199, 20, 4), (240, 45, 31)]
for t, _t in zip(t_, _t_):
  print(t, _t)
  p, d, m = t
  _p, _d, _m = _t
  print(p,d,m)
  print(_p,_d,_m)
  
  
  
  #for t, _t in zip(t_, _t_):