# ************ MAIN FUNCTIONS *******************************************************************************************
# -image_to_blobs()
# -comp_pixel()
# -form_P_()
# -scan_P_()
# -form_seg_()
# -form_blob()
# ***********************************************************************************************************************
from scipy import misc
from time import time
from collections import deque, namedtuple
import numpy as np
import cv2  #Todor

# 13.6.2019
# Todor's modifications for visualisation of the Blobs in color
# Added:
# Collecting blob's on termination in a simple container with:
# global id, blob = [], ave = 25
# OpenCV
# Rendering
# ave = 25
# Best run from cmd line: python frame_blobs.py

def image_to_blobs(image):  # root function, postfix '_' denotes array vs element, prefix '_' denotes higher- vs lower- line variable

    frame = [[0, 0, 0, 0], [], image.shape]  # params, blob_, shape
    dert__ = comp_pixel(image)  # vertically and horizontally bilateral comparison of adjacent pixels
    seg_ = deque()  # buffer of running segments

    for y in range(1, height - 1):  # first and last row are discarded
        P_ = form_P_(dert__[y])  # horizontal clustering
        P_ = scan_P_(P_, seg_, frame)
        seg_ = form_seg_(y, P_, frame)

    while seg_:  form_blob(seg_.popleft(), frame)  # frame ends, last-line segs are merged into their blobs
    return frame  # frame of 2D patterns

    # ---------- image_to_blobs() end -----------------------------------------------------------------------------------


def comp_pixel(p__):  # bilateral comparison between vertically and horizontally consecutive pixels within image

    dert__ = np.empty(shape=(height, width, 4), dtype=int)  # initialize dert__

    dy__ = p__[2:, 1:-1] - p__[:-2, 1:-1]  # vertical comp between rows, first and last column are discarded
    dx__ = p__[1:-1, 2:] - p__[1:-1, :-2]  # lateral comp between columns, first and last row are discarded
    g__ = np.abs(dy__) + np.abs(dx__) - ave  # deviation of gradient, initially approximated as |dy| + |dx|

    dert__[:, :, 0] = p__
    dert__[1:-1, 1:-1, 1] = g__
    dert__[1:-1, 1:-1, 2] = dy__  # first row, last row, first column and last-column are discarded
    dert__[1:-1, 1:-1, 3] = dx__

    return dert__

    # ---------- comp_pixel() end ---------------------------------------------------------------------------------------


def form_P_(dert_):  # horizontally cluster and sum consecutive pixels and their derivatives into Ps

    P_ = deque()  # row of Ps
    i, g, dy, dx = dert_[1]  # first dert
    x0, I, G, Dy, Dx, L = 1, i, g, dy, dx, 1  # P params
    P_dert_ = [(i, g, dy, dx)]  # dert buffer
    _s = g > 0  # sign

    for x, (i, g, dy, dx) in enumerate(dert_[2:-1], start=2):
        s = g > 0
        if s != _s:  # P is terminated and new P is initialized
            P_.append([_s, x0, I, G, Dy, Dx, L, P_dert_])
            x0, I, G, Dy, Dx, L = x, 0, 0, 0, 0, 0
            P_dert_ = []
        # accumulate P params:
        I += i
        G += g
        Dy += dy
        Dx += dx
        L += 1
        P_dert_.append((i, g, dy, dx))
        _s = s  # prior sign

    P_.append([_s, x0, I, G, Dy, Dx, L, P_dert_])  # last P in row
    return P_

    # ---------- form_P_() end ------------------------------------------------------------------------------------------


def scan_P_(P_, seg_, frame):  # integrate x overlaps (forks) between same-sign Ps and _Ps into blob segments

    new_P_ = deque()

    if P_ and seg_:            # if both are not empty
        P = P_.popleft()       # input-line Ps
        seg = seg_.popleft()   # higher-line segments,
        _P = seg[2][-1]        # last element of each segment is higher-line P
        stop = False
        fork_ = []

        while not stop:
            x0 = P[1]          # first x in P
            xn = x0 + P[-2]    # first x in next P
            _x0 = _P[1]        # first x in _P
            _xn = _x0 +_P[-2]  # first x in next _P

            if P[0] == _P[0] and _x0 < xn and x0 < _xn:  # test for sign match and x overlap
                seg[3] += 1
                fork_.append(seg)  # P-connected segments are buffered into fork_

            if xn < _xn:  # _P overlaps next P in P_
                new_P_.append((P, fork_))
                fork_ = []
                if P_:
                    P = P_.popleft()  # load next P
                else:  # terminate loop
                    if seg[3] != 1:  # if roots != 1: terminate seg
                        form_blob(seg, frame)
                    stop = True
            else:  # no next-P overlap
                if seg[3] != 1:  # if roots != 1: terminate seg
                    form_blob(seg, frame)

                if seg_:  # load next _P
                    seg = seg_.popleft()
                    _P = seg[2][-1]
                else:  # if no seg left: terminate loop
                    new_P_.append((P, fork_))
                    stop = True

    while P_:  # terminate Ps and segs that continue at line's end
        new_P_.append((P_.popleft(), []))  # no fork
    while seg_:
        form_blob(seg_.popleft(), frame)  # roots always == 0

    return new_P_

    # ---------- scan_P_() end ------------------------------------------------------------------------------------------


def form_seg_(y, P_, frame):  # convert or merge every P into segment, merge blobs
    new_seg_ = deque()

    while P_:
        P, fork_ = P_.popleft()
        s, x0 = P[:2]
        params = P[2:-1]      # I, G, Dy, Dx, L, Ly
        xn = x0 + params[-1]  # next-P x0 = x0 + L
        params.append(1)      # add Ly

        if not fork_:  # new_seg is initialized with initialized blob
            blob = [s, [0] * (len(params)), [], 1, [y, x0, xn]]  # s, params, seg_, open_segments, box
            new_seg = [y, params, [P], 0, fork_, blob]  # y0, params, Py_, roots, fork_, blob
            blob[2].append(new_seg)
        else:
            if len(fork_) == 1 and fork_[0][3] == 1:  # P has one fork and that fork has one root
                new_seg = fork_[0]
                I, G, Dy, Dx, L, Ly = params
                Is, Gs, Dys, Dxs, Ls, Lys = new_seg[1]  # fork segment params, P is merged into segment:
                new_seg[1] = [Is + I, Gs + G, Dys + Dy, Dxs + Dx, Ls + L, Lys + Ly]
                new_seg[2].append(P)  # Py_: vertical buffer of Ps
                new_seg[3] = 0        # reset roots
                blob = new_seg[-1]

            else:  # if > 1 forks, or 1 fork that has > 1 roots:
                blob = fork_[0][5]
                new_seg = [y, params, [P], 0, fork_, blob]  # new_seg is initialized with fork blob
                blob[2].append(new_seg)   # segment is buffered into blob

                if len(fork_) > 1:        # merge blobs of all forks
                    if fork_[0][3] == 1:  # if roots == 1: fork hasn't been terminated
                        form_blob(fork_[0], frame)  # merge seg of 1st fork into its blob

                    for fork in fork_[1:len(fork_)]:  # merge blobs of other forks into blob of 1st fork
                        if fork[3] == 1:
                            form_blob(fork, frame)

                        if not fork[5] is blob:
                            params, seg_, open_segs, box = fork[5][1:]  # merged blob, omit sign
                            blob[1] = [par1 + par2 for par1, par2 in zip(blob[1], params)]  # sum merging blobs
                            blob[3] += open_segs
                            blob[4][0] = min(blob[4][0], box[0])  # extend box y0
                            blob[4][1] = min(blob[4][1], box[1])  # extend box x0
                            blob[4][2] = max(blob[4][2], box[2])  # extend box xn
                            for seg in seg_:
                                if not seg is fork:
                                    seg[5] = blob  # blobs in other forks are references to blob in the first fork
                                    blob[2].append(seg)  # buffer of merged root segments
                            fork[5] = blob
                            blob[2].append(fork)
                        blob[3] -= 1  # open_segments -= 1: shared with merged blob

        blob[4][1] = min(blob[4][1], x0)  # extend box x0
        blob[4][2] = max(blob[4][2], xn)  # extend box xn
        new_seg_.append(new_seg)

    return new_seg_

    # ---------- form_seg_() end --------------------------------------------------------------------------------------------

def form_blob(term_seg, frame):  # terminated segment is merged into continued or initialized blob (all connected segments)
    global id
    global blobs

    y0s, params, Py_, roots, fork_, blob = term_seg
    blob[1] = [par1 + par2 for par1, par2 in zip(params, blob[1])]
    blob[3] += roots - 1  # number of open segments

    if not blob[3]:  # if open_segments == 0: blob is terminated and packed in frame

        s, [I, G, Dy, Dx, L, Ly], seg_, open_segs, (y0, x0, xn) = blob
        yn = y0s + params[-1]  # yn = y0 + Ly
        map = np.zeros((yn - y0, xn - x0), dtype=bool)  # local map of blob
        for seg in seg_:
            seg.pop()  # remove references to blob
            for y, P in enumerate(seg[2], start=seg[0]):
                x0P = P[1]
                LP = P[-2]
                xnP = x0P + LP
                map[y - y0, x0P - x0:xnP - x0] = True

        frame[0][0] += I
        frame[0][3] += G
        frame[0][1] += Dy
        frame[0][2] += Dx
        frame[1].append(nt_blob(
                                Derts= [I, [[ (G, Dy, Dx, L, Ly, 1, 0, []) ]]],  # Derts[0] = I, Dert[1] = single blob,
                                # rng=1 for comp_range, also layer index = derts[-(rng-1|2)][fa]:
                                # fa=0: sub_layer index: 0 g | 1 ga, none for hypot_g
                                # sub_blob_= [], nested to depth = Derts[index]
                                sign=s,
                                box= (y0, yn, x0, xn),  # boundary box
                                map= map,  # blob boolean map, to compute overlap
                                root_blob=[blob],
                                seg_=seg_,
                                ))

        #Collection for display
        blobs.append([id, (y0, yn, x0, xn), seg_])
        print("terminate form_blob", id)
        id+=1

    # ---------- form_blob() end ----------------------------------------------------------------------------------------

# ************ PROGRAM BODY *********************************************************************************************

ave = 25
DEBUG = True
blobs = [] #collected terminated blobs
path = 'D:\\Cog\\sc.png' #Adjust path to local or whatever

# Load inputs --------------------------------------------------------------------
# image = misc.imread('./../images/raccoon_eye.jpg', flatten=True).astype(int)  # will not be supported by future versions of scipy
#from utils import imread

#image = imread('./../images/raccoon_eye.jpg').astype(int)

#image = cv2.imread('./../images/raccoon_eye.jpg', 0).astype(int)
#image = cv2.imread('D:/Cog/raccoon.jpg',0).astype(int)
#image = cv2.imread('D:\\Cog\\sc.png',0).astype(int)

image = cv2.imread('D:\\Cog\\sc.png',1).astype(int)
b,g,r = cv2.split(image)
channel = b


imageShow = cv2.imread(path, 1) #Images for display - int format is not appropriate
bs,gs,rs = cv2.split(imageShow)
result =  imageShow #np.empty_like(imageShow) #imageShow.copy()
cv2.imshow("source", imageShow) #color
cv2.imshow("b", bs)  #Take the blue channel

cv2.waitKey(2000)

height, width = channel.shape #image.shape

# Main ---------------------------------------------------------------------------
start_time = time()

nt_blob = namedtuple('blob', 'Derts sign box map root_blob seg_')
Dert = namedtuple('Dert', 'G, Dy, Dx, L, Ly, sub_blob_')
Pattern = namedtuple('Pattern', 'sign, x0, I, G, Dy, Dx, L, dert_')
Segment = namedtuple('Segment', 'y, I, G, Dy, Dx, L, Ly, Py_')
Blob = namedtuple('Blob', 'I, Derts, sign, alt, rng, dert__, box, map, root_blob, seg_')
Frame = namedtuple('Frame', 'Dert, dert__')

id = 0

#frame_of_blobs = image_to_blobs(image)
frame_of_blobs = image_to_blobs(channel)

r1 = 50; g1 = 0;  b1 = 0;
r0 = 0; g0 = 0; b0 = 50
cv2.line(result, (0, 0), (200, 200), (255, 255, 255), 1)
bLine = False

for b in blobs:
    i, (y0, yn, x0, xn), seg_ = b
    #   print(b)
    r1 = 50 + i%205
    b0 = 50 + i%205
    for seg in seg_:
        #print(seg)
        #y, [I, G, Dy, Dx, L, Ly], Py_, _, _ = seg
        y, _, Py_, _, _ = seg
        print(y)
        for P in Py_:
            #print(P)
            #s, x0, I, G, Dy, Dx, L, dert_ = P
            s, x0, I, G, Dy, Dx, L, _ = P
            if bLine:
                if s:
                    cv2.line(result, (x0, y), (x0 + L-2, y), (r0, g0, b0), 1)
                    #cv2.line(result, (x0,y0),(x0+L, y0), (r1,g1,b1), 1)
                else:
                    cv2.line(result, (x0, y), (x0 + L-2, y), (r1, g1, b1), 1)
                    #cv2.line(result, (x0, y0), (x0 + L, y0), (r0, g0, b0), 1)
            else:
                for x in range(x0,x0+L):
                    if s:
                        result[y][x] = (r0,g0,b0)
                    else:
                        result[y][x] = (r1, g1, b1)
            y+=1;

cv2.imshow("blobs", result)
#cv2.imwrite("E:\\blobs_cogalg.png", result)
cv2.imwrite("blobs_cogalg.png", result)
cv2.waitKey(0)
# from intra_blob_debug import intra_blob_hypot  # not yet functional, comment-out to run
# frame_of_blobs = intra_blob_hypot(frame_of_blobs)  # evaluate for deeper clustering inside each blob, recursively

# DEBUG --------------------------------------------------------------------------
'''
if DEBUG:
    from utils import *
    draw('./../debug/root_blobs', map_frame(frame_of_blobs))

    from intra_blob_test import intra_blob
    intra_blob(frame_of_blobs[1])
'''
# END DEBUG -----------------------------------------------------------------------

end_time = time() - start_time
print(end_time)
# ************ PROGRAM BODY END ******************************************************************************************
