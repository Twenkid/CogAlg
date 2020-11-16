import time

def inc_rng(a,b, p_):
    print("inc_rng")    
    ip_ = p_  # to differentiate from new p_
    print("after aliasing")
    print("p_ = ", p_)
    print("ip_ = ", ip_)
    pri_s, I, D, V, rv, olp, p_, olp_ = 0, 0, 0, 0, 0, 0, [], [] 
    print("inc_rng")
    print("p_ = ", p_)
    print("ip_ = ", ip_)
    return p_, ip_

x = 5; y = 1;
p_ = [1,2,3,4,5] 
print(p_)
p1_ , ip1_= inc_rng(x, y , p_)
print("caller scope")
print(p1_)
print(ip1_)


def showtime(name, end, start, cnt):
	print(name + ":" + str(end-start) + "(" +  str((end-start)/cnt) + ")" + "[" + str(cnt) + " calls]")
    
start = time.time()   
kray = 10000000
dulgoime_alabalanica = 1.0
drugo_dulgo_imehahahahaha = 1.00001246
for nachalo in range (1, kray)
  dulgoime_alabalanica *=  drugo_dulgo_imehahahahaha
  
  
  