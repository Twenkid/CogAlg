# By Todor, 10/4/2019
# unfold_blob Python semantics test
# https://github.com/khanh93vn/CogAlg/blob/master/frame_2D_alg/intra_comp.py
# You shouldn't remove elements from a list while traversing it with for ... in? (Or in general don't modify it in place)
# It is changing the structure and it seems that elements are skipped. Notice the number of iterations in seg1.
# Maybe because the length of the structure gets shorter than
# as was expected during the initalization of the iteration and 
# some elements are skipped.
# I suggest using separate iterations as in seg3
# Simple test with simpler structures

'''
  Code under test:

        P_ = []                             # buffer for Ps at line y
        for seg in seg_:
            if y < seg[0] + seg[1][0]:      # y < y0 + Ly (y within segment):

                P_.append(seg[2][y - seg[0]])   # append P at line y of seg
            else:                           # y >= y0 + Ly (out of segment):
                seg_.remove(seg)
'''

def seg1(): #as done in unfold_blob ... https://github.com/khanh93vn/CogAlg/blob/master/frame_2D_alg/intra_comp.py
    P_ = []                            
    a = (1,2,3); b=(4,5,6); c = (7,8,9); d=(5,5,5); e = (9,9,9)
    seg_ = [a,b,c,d,e]
    print("==== seg1 ====")
    print("Append to P[]    if < 2")
    print("Remove from seg_ if >= 2")    
    print(seg_)
    i = 1
    for seg in seg_:
       print(i, seg_) 
       x,y,z = seg             
       if x < 2:
         print("O", seg)
         P_.append(seg)
       else:
         print("X", seg)
         seg_.remove(seg)
       i+=1

    print("seg_", seg_)
    print("P_", P_)
    
#as done in intra_comp https://github.com/khanh93vn/CogAlg/blob/master/frame_2D_alg/intra_comp.py
def seg2():
    P_ = []                            
    a = (1,2,3); b=(4,5,6); c = (7,8,9); d=(5,5,5)
    seg_ = [a,b,c,d]
    print("==== seg2 ====")
    print("Append to P[]    if > 2")
    print("Remove from seg_ if <= 2")    
    seg_ = [a,b,c,d]
    print(seg_)

    for seg in seg_:
       print(seg_)       
       x,y,z = seg   
       if x > 2:
         print("O", seg)
         P_.append(seg)
       else:
         print("X", seg)
         seg_.remove(seg)

    print("seg_", seg_)
    print("P_", P_)
    
def seg3(): #Separate traversals
    P_ = []                            
    a = (1,2,3); b=(4,5,6); c = (7,8,9); d=(5,5,5)
    seg_ = [a,b,c,d]
    print("==== seg3 SEPARATE TRAVERSALS ====")
    print("Append if > 2");
    print("Remove if <= 2");
    print("In separate traversals")
    seg_ = [a,b,c,d]
    print(seg_)

    for seg in seg_:       
       x,y,z = seg   
       if x > 2:
         print("O", seg)
         P_.append(seg)
         
    for seg in seg_:
       x,y,z = seg   
       if x <= 2:
         print("X", seg)
         seg_.remove(seg)
                     
    print("seg_", seg_)
    print("P_", P_)
    
def header():
  print("========================")
  print("Test: Don't remove elements from a container that is being traversed with for ... in?")    
  print("========================")
  print()
  

header()

seg1()

seg2()

seg3()
