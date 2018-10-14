#Author: Todor Arnaudov, 10-2018
#Tracking events and states from frame_dblobs (f4.py)

bNoForkContent = True #if False - prints the content of _fork_ - gets big; if True - prints only [] or not empty
global counter, PeqP

class Fork:
    def __init__(s, where, fork, id, roots, id_roots,  x, _x, y):
      s.where = where
      s.fork = fork #.copy()
      s.id = id #fork
      s.roots = roots
      s.id_roots = id_roots
      s.x = x
      s._x = _x
      s.y = y
      s.d = "\t\t"

    def __str__(s):
        '''
        ss = str(s.where) + " _fork_=";
        ss += str(s.fork) + ", id=" + str(s.id);
        ss+=  ",x=" + str(s.x);
        ss+=",_x="+str(s._x);
        ss+=", y=" + str(y);
        '''
        if (bNoForkContent):
          if (len(s.fork) > 0): forkContent = "#"+str(len(s.fork)) #"#>0#"
          else: forkContent = str(s.fork) #[]
        else: forkContent = str(s.fork)
        d = s.d
        #ss = str(s.where) + "fork_id=" + str(s.id) + ",  _fork_=" + forkContent + ", roots=" + str(s.roots) + ", id_roots=" + str(s.id_roots) + ",  x=" + str(s.x) + ",_x="+str(s._x) + ", y=" + str(s.y)
        #ss = str(s.where) + "\tfork_id=" + str(s.id) + "\t_fork_=" + forkContent + "\troots=" + str(s.roots) + "\tid_roots=" + str(s.id_roots) + "\tx=" + str(s.x) + "\tx=" + str(s._x) + "\ty=" + str(s.y)
        ss = str(s.where) + d + "fork_id=" + str(s.id) + d + "_fork_=" + forkContent + d + "roots=" + str(
          s.roots) + d + "id_roots=" + str(s.id_roots) + d + "x=" + str(s.x) + d + "x=" + str(s._x) + d + "y=" + str(s.y)


        #+ " _fork_=" + forkContent;

        return ss

class ForkLog:
    def __init__(s, max=1000):
        s.log = []
        s.count = 0
        s.max = max
        s.bPrint = True#False #print when reaching the limit

    def Add(s, fork):
      if (s.count > s.max):
        if (s.bPrint): s.Print("P:\\cog\\fork3.txt"); s.bPrint = False; return
      s.log.append(fork)

    def Print(s, path):
      print("ForkLog ", len(s.log))
      out =""
      #f = open("P:\\cog\\fork3.txt", "wt")
      f = open(path, "wt")
      import sys
      f.write(str(sys.exc_info()))
      f.write("\n")
      #f.write(forkLog.Print(), f)  # str(forkLog))

      for j,i in enumerate(s.log):
        #print(j,i)
        if (j%66==0): print(j)
        f.write(str(j) + ": "+ str(i) + "\n");

      f.close()
      return out


#forkLog = ForkLog(10000)
forkLog = ForkLog(1000)
