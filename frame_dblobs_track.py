# import cv2
# import argparse
from scipy import misc
from time import time
from collections import deque

'''
Todor, 11-10-2018: Logging code for checking ids and _fork_ ... classes Fork, ForkLog; forkLog, ...
Primitive for now, but can be checked simply via a text editor such as Notepad++: select a number of id and the editor would mark all other lines.
bNoForkContent - doesn't print the content, just if it's not empty
The exceptiion is at:

too many values to unpack (expected 7)
  File "f3.py", line 264, in scan_P_
    hP, frame = form_blob(hP, frame)  # blob (all connected blob segments) += blob segment at hP
  File "f3.py", line 309, in form_blob
    (s, Ls, Is, Ds, Dys, Vs, Vys), Py_, x, xd, root, fork_ = seg  # s, root are ignored
ForkLog 6768
```

The last logged line is that branch:

6767: if roots == 0 id=2109825995208,x=1024,_x=1019, y=410 _fork_=#>0#

That's out of bounds (1023)?


'''   
    frame() is my core algorithm of levels 1 + 2, modified for 2D: segmentation of image into blobs, then search within and between blobs.
    frame_blobs() is frame() limited to definition of initial blobs per each of 4 derivatives, vs. per 2 gradients in current frame().
    frame_dblobs() is updated version of frame_blobs with only one blob type: dblob, to ease debugging, currently in progress.

    Each performs several levels (Le) of encoding, incremental per scan line defined by vertical coordinate y, outlined below.
    value of y per Le line is shown relative to y of current input line, incremented by top-down scan of input image,
    prefix '_' denotes higher-line variable or pattern, vs. same-type lower-line variable or pattern,
    postfix '_' denotes array name, vs. same-name elements of that array:

    1Le, line y:    x_comp(p_): lateral pixel comparison -> tuple of derivatives ders ) array ders_
    2Le, line y- 1: y_comp(ders_): vertical pixel comp -> 2D tuple ders2 ) array ders2_ 
    3Le, line y- 1+ rng*2: form_P(ders2) -> 1D pattern P ) P_  
    4Le, line y- 2+ rng*2: scan_P_(P, _P) -> _P, roots: down-connections, fork_: up-connections between Ps 
    5Le, line y- 3+ rng*2: form_segment: merge vertically-connected _Ps in non-forking blob segments
    6Le, line y- 4+ rng*2+ segment depth: form_blob: merge connected segments into blobs

    These functions are tested through form_P, I am currently debugging scan_P_. 
    All 2D functions (y_comp, scan_P_, etc.) input two lines: higher and lower, convert elements of lower line 
    into elements of new higher line, and displace elements of old higher line into higher function.
    Higher-line elements include additional variables, derived while they were lower-line elements.

    Pixel comparison in 2D forms lateral and vertical derivatives: 2 matches and 2 differences per pixel. 
    They are formed on the same level because average lateral match ~ average vertical match.
    Each vertical and horizontal derivative forms separate blobs, suppressing overlapping orthogonal representations.
    They can also be summed to estimate diagonal or hypot derivatives, for blob orientation to maximize primary derivatives.
    Orientation increases primary dimension of blob to maximize match, and decreases secondary dimension to maximize difference.

    Subsequent union of lateral and vertical patterns is by strength only, orthogonal sign is not commeasurable?
    Initial input line could be 400 for debugging, that area in test image seems to be the most diverse  
'''


def lateral_comp(pixel_):  # comparison over x coordinate: between min_rng of consecutive pixels within each line

    ders_ = []  # tuples of complete derivatives: summation range = rng
    rng_ders_ = deque(maxlen=rng)  # array of tuples within rng of current pixel: summation range < rng
    max_index = rng - 1  # max index of rng_ders_
    pri_d, pri_m = 0, 0  # fuzzy derivatives in prior completed tuple

    for p in pixel_:  # pixel p is compared to rng of prior pixels within horizontal line, summing d and m per prior pixel
        for index, (pri_p, d, m) in enumerate(rng_ders_):

            d += p - pri_p  # fuzzy d: running sum of differences between pixel and all subsequent pixels within rng
            m += min(p, pri_p)  # fuzzy m: running sum of matches between pixel and all subsequent pixels within rng

            if index < max_index:
                rng_ders_[index] = (pri_p, d, m)
            else:
                ders_.append(
                    (pri_p, d + pri_d, m + pri_m))  # completed bilateral tuple is transferred from rng_ders_ to ders_
                pri_d = d;
                pri_m = m  # to complement derivatives of next rng_t_: derived from next rng of pixels

        rng_ders_.appendleft(
            (p, 0, 0))  # new tuple with initialized d and m, maxlen displaces completed tuple from rng_t_

    ders_ += reversed(rng_ders_)  # or tuples of last rng (incomplete, in reverse order) are discarded?
    return ders_


def vertical_comp(ders_, ders2__, _dP_, dframe):
    # comparison between rng vertically consecutive pixels, forming ders2: 2D tuple of derivatives per pixel

    dP = 0, 0, 0, 0, 0, 0, 0, []  # lateral difference pattern = pri_s, L, I, D, Dy, V, Vy, ders2_
    dP_ = deque()  # line y - 1+ rng*2
    dbuff_ = deque()  # line y- 2+ rng*2: _Ps buffered by previous run of scan_P_
    new_ders2__ = deque()  # 2D: line of ders2_s buffered for next-line comp

    x = 0  # lateral coordinate of current pixel
    max_index = rng - 1  # max ders2_ index
    min_coord = rng * 2 - 1  # min x and y for form_P input: ders2 from comp over rng*2 (bidirectional: before and after pixel p)
    dy, my = 0, 0  # for initial rng of lines, to reload _dy, _vy = 0, 0 in higher tuple

    for (p, d, m), (ders2_, _dy, _my) in zip(ders_,
                                             ders2__):  # pixel compared to rng _p s in ders2_, summing dy and my per higher pixel
        x += 1
        index = 0
        for (
        _p, _d, dy, _m, my) in ders2_:  # vertical derivatives are incomplete; prefix '_' denotes higher-line variable

            dy += p - _p  # fuzzy dy: running sum of differences between pixel and all lower pixels within rng
            my += min(p, _p)  # fuzzy my: running sum of matches between pixel and all lower pixels within rng

            if index < max_index:
                ders2_[index] = (_p, d, dy, m, my)

            elif x > min_coord and y > min_coord + 400:

                _v = _m - abs(
                    d) - ave  # _m - abs(d): projected m cancelled by negative d: d/2, + projected rdn value of overlapping dP: d/2
                vy = my + _my - abs(dy) - ave
                ders2 = _p, _d, dy + _dy, _v, vy
                dP, dP_, dbuff_, _dP_, dframe = form_P(ders2, x, dP, dP_, dbuff_, _dP_, dframe)

            index += 1

        ders2_.appendleft(
            (p, d, 0, m, 0))  # initial dy and my = 0, new ders2 replaces completed t2 in vertical ders2_ via maxlen
        new_ders2__.append(
            (ders2_, dy, my))  # vertically-incomplete 2D array of tuples, converted to ders2__, for next-line ycomp

    if y > min_coord + 400:  # not-terminated P at the end of each line is buffered or scanned:

        if y == rng * 2 + 400:  # _P_ initialization by first line of Ps, empty until vertical_comp returns P_
            dP_.append([dP, x, [], 0])  # empty _fork_ in the first line of hPs, x-1: delayed P displacement
        else:
            dP_, dbuff_, _dP_, dframe = scan_P_(x, dP, dP_, dbuff_, _dP_, dframe)  # scans higher-line Ps for contiguity

    return new_ders2__, dP_, dframe  # extended in scan_P_; net_s are packed into frames


def form_P(ders2, x, P, P_, buff_, _P_,
           frame):  # initializes, accumulates, and terminates 1D pattern: dP | vP | dyP | vyP

    p, d, dy, v, vy = ders2  # 2D tuple of derivatives per pixel, "y" denotes vertical derivatives:
    s = 1 if d > 0 else 0  # core = 0 is negative: no selection?

    if s == P[0] or x == rng * 2:  # s == pri_s or initialized: P is continued, else terminated:
        pri_s, L, I, D, Dy, V, Vy, ders2_ = P
    else:
        if y == rng * 2 + 400:  # _P_ initialization by first line of Ps, empty until vertical_comp returns P_
            P_.append([P, x - 1, [], 0])  # hP container to maintain fork refs, empty _fork_ in the first line of _Ps
        else:
            P_, buff_, _P_, frame = scan_P_(x - 1, P, P_, buff_, _P_, frame)  # scans higher-line Ps for contiguity
            # x-1: last P displacement
        L, I, D, Dy, V, Vy, ders2_ = 0, 0, 0, 0, 0, 0, []  # new P initialization

    L += 1  # length of pattern
    I += p  # summed input and derivatives are accumulated as P and alt_P parameters, continued or initialized:
    D += d  # lateral D
    Dy += dy  # vertical D
    V += v  # lateral V
    Vy += vy  # vertical V
    ders2_.append(ders2)  # ders2s are buffered for oriented rescan and incremental range | derivation comp

    P = s, L, I, D, Dy, V, Vy, ders2_  # P is disposable, only _P[0] id is conserved

    return P, P_, buff_, _P_, frame  # accumulated within line, P_ is a buffer for conversion to _P_


#@Todor, 11-10-2018

bNoForkContent = True #if False - prints the content of _fork_ - gets big; if True - prints only [] or not empty
class Fork:
    def __init__(s, where, fork, id,  x, _x, y):
      s.where = where
      s.fork = fork #.copy()
      s.id = id
      s.x = x
      s._x = _x
      s.y = y
    def __str__(s):
        '''
        ss = str(s.where) + " _fork_=";
        ss += str(s.fork) + ", id=" + str(s.id);
        ss+=  ",x=" + str(s.x);
        ss+=",_x="+str(s._x);
        ss+=", y=" + str(y);
        '''
        if (bNoForkContent):
          if (len(s.fork) > 0): forkContent = "#>0#"
          else: forkContent = str(s.fork) #[]
        else: forkContent = str(s.fork)

        ss = str(s.where) + "id=" + str(s.id) + ",x=" + str(s.x) + ",_x="+str(s._x) + ", y=" + str(y) + " _fork_=" + forkContent;
        return ss


class ForkLog:
    def __init__(s, max=1000):
        s.log = []
        s.count = 0
        s.max = max

    def Add(s, fork):
      if (s.count > s.max): return
      s.log.append(fork)

    def Print(s, path):
      print("ForkLog ", len(s.log))
      out =""
      #f = open("P:\\cog\\fork.txt", "wt")
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


forkLog = ForkLog(10000)

def scan_P_(x, P, P_, _buff_, _P_, frame):  # P scans shared-x-coordinate _Ps in _P_, forms overlaps

    buff_ = deque()  # new buffer for displaced hPs, for scan_P_(next P)
    fork_ = []  # hPs connected to input P
    ini_x = 0  # only to start while loop

    while ini_x <= x:  # while horizontal overlap between P and _P, after that: P -> P_
        if _buff_:
            hP = _buff_.popleft()  # load _P buffered in prior run of scan_P_, if any
            _P, _x, _fork_, roots = hP;  forkLog.Add(Fork("if _buff_", _fork_, id(_fork_), x, _x, y))
        elif _P_:
            hP = _P_.popleft()
            _P, _x, _fork_, roots = hP; forkLog.Add(Fork("elif _P_:", _fork_, id(_fork_), x, _x, y))  # higher-line P container, for return as fork ref
            # roots = 0: number of Ps connected to current _P(pri_s, L, I, D, Dy, V, Vy, ders2_)
        else:
            break
        _L = _P[1];
        ini_x = _x - _L;
        ave_x = _x - _L // 2  # initial and average x coordinates of _P

        if P[0] == _P[0]:  # if s == _s: core sign match, + selective inclusion if contiguity eval?
            roots += 1  # always > 1?
            hP[3] = roots  # nothing else is modified
            fork_.append(hP)  # P-connected hPs, appended with blob and converted to Py_ after P_ scan

        if _x > x - P[
            1]:  # x overlap between hP and next P: hP is buffered for next scan_P_, else _P is included in unique blob segment
            buff_.append(hP)
        else:
            if roots == 0:  # test after full scan over P_, never happens beyond y==5?
                id(roots)
            ini = 1
            forkLog.Add(Fork("before if y>rng*2:", _fork_, id(_fork_), x, _x, y))
            if y > rng * 2 + 1 + 400:  # beyond 1st line of _fork_ Ps, else: blob segment ini only
                if len(_fork_) == 1:
                    if len(_fork_[0]) == 6:  # len seg == 6, len hP == 4
                        if _fork_[0][4][
                            0] == 1:  # _fork roots, see if ini = 1, second [] is a _P container; never happens
                            seg = form_seg(_P, _fork_[0], ave_x)  # _P is added to blob segment at _fork_[0]
                            del (hP[:]);
                            hP += seg
                            ini = 0
                    else:
                        break
            if ini == 1:
                del (hP[:])  # blob segment [seg_vars, Py_, ave_x, Dx, root, _fork_] is initialized at not-included hP:
                hP += _P, [(_P, ave_x, 0)], ave_x, 0, [roots, [], (
                0, 0, 0, 0, 0, 0, 0, 0, 0)], _fork_  # also replaces fork_[len(fork_)]

            try:
                if roots == 0:  # only happens in initialized _Ps: always len(_fork_) == 0?
                    forkLog.Add(Fork("if roots == 0 ", _fork_, id(_fork_), x, _x, y))
                    if len(_fork_):  # blob ini per seg, above
                        if len(hP) == 6:
                            hP, frame = form_blob(hP, frame)  # blob (all connected blob segments) += blob segment at hP
                    else:
                        frame = form_frame(hP[4][1], hP[4][2], frame,
                                           ave_x)  # (root_, blob, frame, x): blob, root_ packed in frame
            except Exception as e:
                import sys, traceback
                print(e)
                a,b,trace = sys.exc_info()
                print(a)
                print(b)
                traceback.print_tb(trace)
                #print(t)
                forkLog.Print("P:\\cog\\fork2.txt")
                exit(0)
                #f = open("P:\\cog\\fork.txt", "wt")
                #f.write(str(forkLog))
                #f.close()

    if len(fork_) == 0:  # test after full scan over _P_, never happens beyond y==5?
        id(fork_)  # roots are reciprocals in _Ps, over multiple Ps

    buff_ += _buff_  # _buff_ is likely empty
    P_.append([P, x, fork_, 0])  # P with no overlap to next _P is buffered for next-line scan_P_, converted to _P

    return P_, buff_, _P_, frame  # _P_ and buff_ contain only _Ps with _x => next x


def form_seg(P, seg, x):  # continued or initialized blob segment is incremented by attached _P
    s, L, I, D, Dy, V, Vy, ders2_ = P  # s is identical, ders2_ is a replacement, x, fork_, roots are ignored
    (s, Ls, Is, Ds, Dys, Vs, Vys,
     ignore), Py_, _x, xD, root, fork_ = seg  # fork_ assigned at ini only, roots at form_blob?
    Ls += L
    Is += I
    Ds += D
    Dys += Dy
    Vs += V
    Vys += Vy
    xd = x - _x
    xD += xd  # for segment normalization and orientation eval, | += |xd| for curved max_L norm, orient?
    Py_.append(((s, L, I, D, Dy, V, Vy, ders2_), x, xd))
    seg = [(s, Ls, Is, Ds, Dys, Vs, Vys, ignore), Py_, x, xD, root, fork_]  # initial ders2_ is ignored
    return seg


def form_blob(seg, frame):  # continued or initialized blob is incremented by attached blob segment and its root_
    (s, Ls, Is, Ds, Dys, Vs, Vys), Py_, x, xd, root, fork_ = seg  # s, root are ignored
    for index, _seg in enumerate(fork_):  # _segment per fork = _P, Py_, ave_x, Dx, root, _fork_

        _roots, _root_, _blob = _seg[0][4]  # if len(seg[0]) == 6, or all displaced forks, if any?
        s, Lb, Ib, Db, Dyb, Vb, Vyb, xD, yD = _blob  # incomplete blob doesn't have ave_x
        Lb += Ls
        Ib += Is
        Db += Ds
        Dyb += Dys
        Vb += Vs
        Vyb += Vys
        xD += xd
        yD += len(Py_)  # only if max Py_?
        _blob = s, Lb, Ib, Db, Dyb, Vb, Vyb, xD, yD
        _roots -= 1  # root segment is terminated and appended:
        _root_.append(seg)

        if _roots == 0:
            if len(_seg[5]):  # _fork_
                _seg, frame = form_blob(_seg,
                                        frame)  # recursive higher-level segment -> blob inclusion and termination test
            else:
                frame = form_frame(_seg[0][4][1], _seg[0][4][2], frame,
                                   x)  # all connected roots and forks terminate, blob packed in frame

                # also roots and forks of co-forks per seg: right_count and left_1st (for transfer at roots+forks term)
                # right: cont roots and len(_fork_), each summed at current subb term, for blob term eval?
                # subb inclusion is unique, in leftmost continued fork?

        _seg[0][4] = [_roots, _root_, _blob]
        seg[5][index] = _seg  # return to fork
    return [seg, frame]  # top segment includes rep of partial blob


def form_frame(root_, blob, frame, x):
    s, Lb, Ib, Db, Dyb, Vb, Vyb, xD, yD = blob
    (Lf, If, Df, Dyf, Vf, Vyf, xDf,
     yDf), blob_ = frame  # to compute averages of dframe, redundant for same-scope alt_frames?
    Lf += Lb
    If += Ib
    Df += Db
    Dyf += Dyb
    Vf += Vb
    Vyf += Vyb
    xDf += xD  # for frame normalization, orient eval, += |xd| for curved max_L?
    yDf += yD
    blob_.append([(Lb, Ib, Db, Dyb, Vb, Vyb, x - xD // 2, xD, y, yD), root_])  # xD to normalize blob before comp_P
    frame = (Lf, If, Df, Dyf, Vf, Vyf, xDf, yDf), blob_
    return frame


def image_to_blobs(
        image):  # postfix '_' denotes array vs. element, prefix '_' denotes higher-line vs. lower-line variable

    _P_ = deque()  # higher-line same- d-, v-, dy-, vy- sign 1D patterns
    frame = (0, 0, 0, 0, 0, 0, 0, 0), []  # (Lf, If, Df, Dyf, Vf, Vyf, xDf, yDf), blob_
    global y
    y = 400  # initial input line, set at 400 as that area in test image seems to be the most diverse

    ders2_ = deque(maxlen=rng)  # vertical buffer of incomplete derivatives tuples, for fuzzy ycomp
    ders2__ = []  # vertical buffer + horizontal line: 2D array of 2D tuples, deque for speed?
    pixel_ = image[0, :]  # first line of pixels at y == 0
    ders_ = lateral_comp(pixel_)  # after part_comp (pop, no t_.append) while x < rng?

    for (p, d, m) in ders_:
        ders2 = p, d, 0, m, 0  # dy, my initialized at 0
        ders2_.append(ders2)  # only one tuple per first-line ders2_
        ders2__.append((ders2_, 0, 0))  # _dy, _my initialized at 0

    for y in range(401, Y):  # or Y-1: default term_blob in scan_P_ at y = Y?
        if (y%30 ==0): print(y + ": ", end="")
        pixel_ = image[y, :]  # vertical coordinate y is index of new line p_
        ders_ = lateral_comp(pixel_)  # lateral pixel comparison
        ders2__, _P_, frame = vertical_comp(ders_, ders2__, _P_, frame)  # vertical pixel comparison

    # frame ends, last vertical rng of incomplete ders2__ is discarded,
    # vertically incomplete P_ patterns are still inputted in scan_P_?
    return frame  # frame of 2D patterns to be outputted to level 2


# pattern filters: eventually updated by higher-level feedback, initialized here as constants:

rng = 2  # number of leftward and upward pixels compared to each input pixel
ave = 63 * rng * 2  # average match: value pattern filter
ave_rate = 0.25  # average match rate: ave_match_between_ds / ave_match_between_ps, init at 1/4: I / M (~2) * I / D (~2)

image = misc.face(gray=True)  # read image as 2d-array of pixels (gray scale):
image = image.astype(int)
Y, X = image.shape  # image height and width

# or:
# argument_parser = argparse.ArgumentParser()
# argument_parser.add_argument('-i', '--image', help='path to image file', default='./images/raccoon.jpg')
# arguments = vars(argument_parser.parse_args())
# image = cv2.imread(arguments['image'], 0).astype(int)

start_time = time()
frame_of_blobs = image_to_blobs(image)
end_time = time() - start_time
print(end_time)



 #forkLog.Add(Fork("if roots == 0 ", _fork_, id(_fork_), x, _x, y))
