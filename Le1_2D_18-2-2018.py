from scipy import misc
from collections import deque
import numpy as np

''' Level 1 with patterns defined by the sign of quadrant gradient: modified core algorithm of levels 1+2.
    Pixel comparison in 2D forms lateral and vertical derivatives: 2 matches and 2 differences per pixel. 
    They are formed on the same level because average lateral match ~ average vertical match.
    Pixels are discrete samples of continuous image, so rightward and downward derivatives per pixel are 
    equally representative samples of 0-90 degree quadrant gradient: minimal unique unit of 2D gradient. 
    Such gradient is computed as the average of these two orthogonally diverging derivatives.
    2D blobs are defined by same-sign quadrant gradient, of value for vP or difference for dP.
    Level 1 performs several steps of incremental encoding, per line defined by vertical coordinate y:
    y: comp(p_): lateral comp -> tuple t,
    y- 1: ycomp(t_): vertical comp -> quadrant t2,
    y- 1: form_P(P): lateral combination -> 1D pattern P,  
    y- 2: scan_P_(P, _P): computes vertical continuity between 1D Ps  
    y- >2: form_blob: displaces y-2 P into linked 2D-contiguous blob
    y- >3: term_blob: terminated blobs are evaluated for comp_P and form_PP, -> 2D patterns,
           resulting PPs are evaluated for blob re-orientation, re-scan, and PP re-consolidation 
    All 2D functions (ycomp, scan_P_, etc.) input two lines: relatively higher and lower.
    Higher-line patterns include additional variables, derived while they were lower-line patterns
    
    postfix '_' denotes array name, vs. same-name elements of that array 
    prefix '_' denotes higher-line variable or pattern '''


def comp(p_):  # comparison of consecutive pixels within line forms tuples: pixel, match, difference

    t_ = []  # complete fuzzy tuples: summation range = rng
    it_ = deque(maxlen=rng)  # incomplete fuzzy tuples: summation range < rng

    for p in p_:
        index = 0

        for it in it_:  # incomplete tuples, with summation range from 0 to rng
            pri_p, fd, fm = it

            d = p - pri_p  # difference between pixels
            m = min(p, pri_p)  # match between pixels

            fd += d  # fuzzy d: sum of ds between p and all prior ps within it_
            fm += m  # fuzzy m: sum of ms between p and all prior ps within it_

            it_[index] = pri_p, fd, fm
            index += 1

        if len(it_) == rng:  # or while x < rng: icomp(){ p = pop(p_).., no t_.append?
            t_.append((pri_p, fd, fm))  # completed tuple is transferred from it_ to t_

        it_.appendleft((p, 0, 0))  # new prior tuple, fd and fm are initialized at 0

    t_ += it_  # last number = rng of tuples remain incomplete
    return t_


def ycomp(t_, t2__, _vP_, _dP_):  # vertical comparison between pixels, forms 2D tuples t2

    vP_ = []; vP = [0,0,0,0,0,0,0,0,[]]  # pri_s, I, D, Dy, M, My, G, Olp, e_
    dP_ = []; dP = [0,0,0,0,0,0,0,0,[]]  # pri_s, I, D, Dy, M, My, G, Olp, e_

    vg_blob_, dg_blob_ = [],[]  # output line of blobs, vertical concat -> frame in frame()

    x = 0; new_t2__ = []   # t2_ buffer: 2D array
    olp, ovG, odG = 0,0,0  # len of overlap between vP and dP, gs summed over olp, all shared

    for t, t2_ in zip(t_, t2__):  # compares vertically consecutive pixels, forms quadrant gradients

        p, d, m = t; index = 0; x += 1

        for t2 in t2_:
            pri_p, _d, fdy, _m, fmy = t2

            dy = p - pri_p  # vertical difference between pixels
            my = min(p, pri_p)  # vertical match between pixels

            fdy += dy  # fuzzy dy: sum of dys between p and all prior ps within t2_
            fmy += my  # fuzzy my: sum of mys between p and all prior ps within t2_

            t2_[index] = pri_p, _d, fdy, _m, fmy
            index += 1

        if len(t2_) == rng:  # or while y < rng: i_ycomp(){ t2_ = pop(t2__), t = pop(t_)., no form_P?

            dg = _d + fdy  # d gradient
            vg = _m + fmy - ave  # v gradient
            t2 = pri_p, _d, fdy, _m, fmy  # complete 2D tuple moved from t2_ to form_P:

            # form 1D patterns vP and dP: horizontal spans of same-sign vg or dg, with associated vars:

            olp, ovG, odG, vP, dP, vP_, _vP_, vg_blob_ = \
            form_P(1, t2, vg, dg, olp, ovG, odG, vP, dP, vP_, _vP_, vg_blob_, x)

            olp, odG, ovG, dP, vP, dP_, _dP_, dg_blob_ = \
            form_P(0, t2, dg, vg, olp, odG, ovG, dP, vP, dP_, _dP_, dg_blob_, x)

        t2_.appendleft((p, d, 0, m, 0))  # initial fdy and fmy = 0, new t2 replaces completed t2 in t2_
        new_t2__.append(t2_)
        
    # line ends, vP and dP term, no init, inclusion with incomplete lateral fd and fm:

    if olp:  # if vP x dP overlap len > 0, incomplete vg - ave / (rng / X-x)?

        odG *= ave_k; odG = odG.astype(int)  # ave_k = V / I, to project V of odG

        if ovG > odG:  # comp of olp vG and olp dG, == goes to vP: secondary pattern?
            dP[7] += olp  # overlap of lesser-oG vP or dP, or P = P, Olp?
        else:
            vP[7] += olp  # to form rel_rdn = alt_rdn / len(e_)

    if y + 1 > rng:  # starting with the first line of complete t2s

        vP_, _vP_, vg_blob_ = scan_P_(0, vP, vP_, _vP_, vg_blob_, x)  # returns empty _vP_
        dP_, _dP_, dg_blob_ = scan_P_(1, dP, dP_, _dP_, dg_blob_, x)  # returns empty _dP_

    return new_t2__, _vP_, _dP_, vg_blob_, dg_blob_  # blob_s are extended in scan_P_

    # poss alt_: top P alt = Olp, oG, alt_oG: to remove if hLe demotion and alt_oG < oG?
    # P_ can be redefined as np.array ([P, alt_, roots, forks) to increment without init?


def form_P(typ, t2, g, alt_g, olp, oG, alt_oG, P, alt_P, P_, _P_, blob_, x):

    # forms 1D dP or vP, then scan_P_ adds forks in _P fork_s and accumulates blob_

    p, d, dy, m, my = t2  # 2D tuple of quadrant variables per pixel
    pri_s, I, D, Dy, M, My, G, Olp, e_ = P  # initial pri_ vars = 0, or skip form?

    s = 1 if g > 0 else 0  # g = 0 is negative?
    if s != pri_s and x > rng + 2:  # P is terminated

        if typ: alt_oG *= ave_k; alt_oG = alt_oG.astype(int)  # ave V / I, to project V of odG
        else: oG *= ave_k; oG = oG.astype(int)               # same for h_der and h_comp eval?

        if oG > alt_oG:  # comp between overlapping vG and dG
            Olp += olp  # olp is assigned to the weaker of P | alt_P, == -> P: local access
        else:
            alt_P[7] += olp

        P = pri_s, I, D, Dy, M, My, G, Olp, e_  # no ave * alt_rdn / e_: adj < cost
        P_, _P_, blob_ = scan_P_(typ, P, P_, _P_, blob_, x)  # scan over contiguous higher-level _Ps

        I, D, Dy, M, My, G, Olp, e_ = 0,0,0,0,0,0,0,[]  # P and olp initialization
        olp, oG, alt_oG = 0,0,0

    # continued or initialized vars are accumulated:

    olp += 1  # len of overlap to stronger alt-type P, accumulated until P or _P terminates
    oG += g; alt_oG += alt_g  # for eval to assign olp to alt_rdn of vP or dP

    I += p    # pixels summed within P
    D += d    # lateral D, for P comp and P2 orientation
    Dy += dy  # vertical D, for P2 normalization
    M += m    # lateral D, for P comp and P2 orientation
    My += my  # vertical M, for P2 orientation
    G += g    # d or v gradient summed to define P value, or V = M - 2a * W?

    if typ:
        e_.append((p, g, alt_g))  # g = v gradient, for selective incremental range comp
    else:
        e_.append(g)  # g = d gradient and pattern element, for selective incremental derivation

    P = [s, I, D, Dy, M, My, G, Olp, e_]
    return olp, oG, alt_oG, P, alt_P, P_, _P_, blob_  # accumulated in ycomp


def scan_P_(typ, P, P_, _P_, blob_, x):  # P scans overlapping _Ps in _P_, forms overlapping Gs

    buff_, root_, selmax_ = [],[],[]
    s, I, D, Dy, M, My, G, Olp, e_ = P  # unused Olp: 1D overlap by stronger alt Ps, no unpack?

    ix = x - len(e_)  # initial x of P
    _ix = 0  # initialized ix of _P displaced from _P_ by last scan_P_

    while x >= _ix:  # P to _P match eval, while horizontal overlap between P and _P_:

        ex = x  # ex is lateral coordinate of loaded P element
        _P = _P_.popleft()  # _P = _P in y-2, blob in y-3, fork_ in y-1

        if s == _P[0][0]:  # if s == _s: vg or dg sign match, fork_.append eval

            oG = 0  # fork gradient overlap: oG += g (distinct from alt_P oG)
            while ex > _P[0][1]: # ex > _ix
        
                for e in e_:  # accumulation of oG between P and _P:

                    if typ: oG += e[1]  # if vP: e = p, g, alt_g
                    else: oG += e  # if dP: e = g
                    ex += 1

            if oG > ave * 16: # !max _P: likely termination, new blob and fork cost?
                _P[2].append((oG, P))  # fork_.append(P)

            elif oG > ave * 4: # max _P: summation only, form_blob cost?
                selmax_.append((oG, P))  # for select and _fork_.append after full scan

        if _P[0][2] > ix:  # if _x > ix:
            buff_.append(_P)  # _P is buffered for scan_P_(next P)

        else:  # no overlap between _P and next P, term_blob eval, form_blob:

            if (_P[2] != 1 and y > rng) or y == Y - 1:  # if fork_==0 | >1: segment term

                blob = term_blob(typ, _P, blob_)  # blob = _P[1]: single segment per _P
                blob_.append(blob); init = 1  # blob_ is line y - >2

            else: init = 0
            form_blob(_P, init)  # default, typ is of linked _P, for potential comp_P?

    # no overlap between P and next _P

    if root_== 0 and selmax_:  # selection of root _P by max oG, for symmetric fork_

        root = max(selmax_) # same as root = max(selmax_, key=lambda selmax: selmax[0])?
        root[1].append((root[0], P))  # (oG, P) added to fork_ of max root _P

    P = s, ix, x, I, D, Dy, M, My, G, Olp, e_  # P becomes _P

    P_.append((P, [], []))  # initial blob & fork_, _P_ = P_ for next-line scan_P_()
    buff_ += _P_  # excluding displaced _Ps

    return P_, buff_, blob_  # _P_ = buff_ for scan_P_(next P)


''' sequential displacement and higher-line inclusion at line end:
    y- 1: P
    y- 2: _P, fork_
    y- >2: _P blob of variable depth
    y- >3: terminated blob segment, fork_; global term if root_==0: same blob OG eval?
    
    no rdn P rep and select:
    root_.sort(key=lambda fork: fork[0], reverse=True) # max-to-min oG,
    select_ = deque(); rdn = 1  # number of select forks per P, or no rdn?
    root = root_.pop()
    while root_ and (root[0] > ave * 12):  # oG = fork[0], no oG rdn, fixed ref and 11 var cost?
        select_.appendleft(root); rdn += 1  # inclusion if match, no neg forks?
        root = root_.pop()  # no: fork = rdn + alt_rdn / len(e_), _P: float cost?
    init = 0 if select_ == 1 else 1  # new blob if select_== 0, new segment if select_> 1?
    for root in select_:
        root = form_blob(P, root, rdn, init)  # P added to fork segment, in root-connected blob
        root_.appendleft(root)  # not-selected roots are out of root_
'''

def form_blob(_P, init):  # _P inclusion into blob, init if fork_ != 1

    oG, _P = _P  # or oG arg?
    _P, blob, fork_ = _P  # fork_ is unused
    s, ix, x, I, D, Dy, M, My, G, Olp, e_ = _P  # no rdn = e_ / Olp + blob_rdn

    if init:  # new blob or segment:

        L2 = len(e_)  # no separate e2_: Py_( P( e_? overlap / comp_P only?
        I2 = I
        D2 = D; Dy2 = Dy
        M2 = M; My2 = My
        G2 = G
        OG = oG  # vertical contiguity for comp_P eval
        Olp2 = Olp
        Py_ = [_P]  # + oG? vertical array of patterns within a blob

    else:  # blob extend:
        L2, I2, D2, Dy2, M2, My2, G2, OG, Olp2, rdn2, Py_ = blob

        L2 += len(e_)
        I2 += I
        D2 += D; Dy2 += Dy
        M2 += M; My2 += My
        G2 += G
        OG += oG
        Olp2 += Olp
        Py_.append(_P)  # + oG?

    blob = s, L2, I2, D2, Dy2, M2, My2, G2, OG, Olp2, Py_, fork_

    return blob # _P summed into blob


def term_blob(typ, _P, blob_):  # blob eval for comp_P, only if complete term: root_ and fork_ == 0?

    P, blob, fork_ = _P  # top-down fork_, no root_: no rdn, eval of _fork_ cost only?
    Py_ = blob[11]

    if blob[8] > ave*8 and Py_ > 2:  # blob value: OG += oG, cost = ave*8 or ave_OG?

        vPP, dPP = [], []  # 2D value P and difference P
        _P = Py_.popleft()  # initial comparand, top-down order

        while Py_: # comp_P starts from 2nd P

            P = Py_.popleft()
            _P, _vs, _ds = comp_P(typ, P, _P)  # per blob, before orient

            while Py_:  # form_PP starts from 3rd P

                P = Py_.popleft()
                P, vs, ds = comp_P(typ, P, _P)

                vPP = form_PP(typ, P, vs, _vs, vPP); blob_.append(vPP)
                dPP = form_PP(typ, P, ds, _ds, dPP); blob_.append(dPP)

                _P = P; _vs = vs; _ds = ds

    return blob_  # blob | PP_, comp_P may continue over fork_, after comp_segment?


def comp_P(typ, P, _P):  # forms vertical derivatives of P vars, also from conditional DIV comp

    s, ix, x, I, D, Dy, M, My, G, oG, Olp, e_ = P
    _s, _ix, _x, _I, _D, _Dy, _M, _My, _G, _oG, _Olp, _e_ = _P

    ddx = 0  # optional, 2Le norm / D? s_ddx and s_dL correlate, s_dx position and s_dL dimension don't?
    ix = x - len(e_)  # initial coordinate of P; S is generic for summed vars I, D, M:

    dx = x - len(e_)/2 - _x - len(_e_)/2  # Dx? comp(dx), ddx = Ddx / h? dS *= cos(ddx), mS /= cos(ddx)?
    mx = x - _ix; if ix > _ix: mx -= ix - _ix  # mx = x olp, - a_mx -> vxP, distant P mx = -(a_dx - dx)?

    dL = len(e_) - len(_e_); mL = min(len(e_), len(_e_))  # relative olp = mx / L? ext_miss: Ddx + DL?
    dI = I - _I; mI = min(I, _I)
    dD = D - _D; mD = min(D, _D)
    dM = M - _M; mM = min(M, _M)  # no G comp: y-derivatives are incomplete, no alt_ comp: rdn only?

    Pd = ddx + dL + dI + dD + dM  # defines dPP; var_P form if PP form, term if var_P or PP term;
    Pm = mx + mL + mI + mD + mM   # defines vPP; comb rep value = Pm * 2 + Pd? group by y_ders?

    if dI * dL > div_a: # DIV comp: cross-scale d, neg if cross-sign, nS = S * rL, ~ rS,rP: L defines P
                        # no ndx, yes nmx: summed?

        rL = len(e_) / len(_e_)  # L defines P, SUB comp of rL-normalized nS:
        nI = I * rL; ndI = nI - _I; nmI = min(nI, _I)  # vs. nI = dI * nrL?
        nD = D * rL; ndD = nD - _D; nmD = min(nD, _D)
        nM = M * rL; ndM = nM - _M; nmM = min(nM, _M)

        Pnm = mx + nmI + nmD + nmM  # normalized m defines norm_vPP, as long as rL is computed
        if Pm > Pnm: nvPP_rdn = 1; vPP_rdn = 0 # added to rdn, or diff alt, olp, div rdn?
        else: vPP_rdn = 1; nvPP_rdn = 0

        Pnd = ddx + ndI + ndD + ndM  # normalized d defines norm_dPP or ndPP
        if Pd > Pnd: ndPP_rdn = 1; dPP_rdn = 0  # value = D | nD
        else: dPP_rdn = 1; ndPP_rdn = 0

        div_f = 1
        nvars = Pnm, nmI, nmD, nmM, vPP_rdn, nvPP_rdn, \
                Pnd, ndI, ndD, ndM, dPP_rdn, ndPP_rdn

    else:
        div_f = 0  # DIV comp flag
        nvars = 0  # DIV + norm derivatives

    vs = 1 if Pm > ave * 5 > 0 else 0  # comp cost = ave * 5, or rep cost: n vars per P?
    ds = 1 if Pd > 0 else 0

    P_ders = Pm, Pd, mx, dx, mL, dL, mI, dI, mD, dD, mM, dM, div_f, nvars
    P = P, P_ders

    return P, vs, ds  # for inclusion in vPP_, dPP_ by form_PP:

''' 
    no DIV comp(L): match is insignificant and redundant to mS, mLPs and dLPs only?:
    if dL: nL = len(e_) // len(_e_)  # L match = min L mult
    else: nL = len(_e_) // len(e_)
    fL = len(e_) % len(_e_)  # miss = remainder 
    form_PP at fork_eval after full rdn: A = a * alt_rdn * fork_rdn * norm_rdn, 
    form_pP (parameter pattern) in +vPPs only, then cost of adjust for pP_rdn?
    comp_P is not fuzzy: x & y vars are already fuzzy?  
    eval per fork, PP, or yP, not per comp
    no comp aS: m_aS * rL cost, minor cpr / nL? no DIV S: weak nS = S // _S; fS = rS - nS  
    or aS if positive eV (not eD?) = mx + mL -ave:
    aI = I / len(e_); dI = aI - _aI; mI = min(aI, _aI)  
    aD = D / len(e_); dD = aD - _aD; mD = min(aD, _aD)  
    aM = M / len(e_); dM = aM - _aM; mM = min(aM, _aM)
    d_aS comp if cs D_aS, _aS is aS stored in _P, S preserved to form hP SS?
    iter dS - S -> (n, M, diff): var precision or modulo + remainder?  '''


def form_PP(typ, P, Ps, _Ps, PP):  # forms vPPs | dPPs, and pPs within each

    P, P_ders = P  # + fork_ per segment: possibly merged into PP?
    s, ix, x, I, D, Dy, M, My, G, Olp, e_ = P
    Pm, Pd, mx, dx, mL, dL, mI, dI, mD, dD, mM, dM, div_f, nvars = P_ders

    if Ps != _Ps:
        PP = term_PP(PP)  # then eval for reorient, rescan, recursion?

        dxP_ = []; dx2 = dx, dxP_ # or ddx_P: known x match?
        mxP_ = []; mx2 = mx, mxP_

        LP_ = []; L2 = len(e_), LP_  # no rL, fL?
        IP_ = []; I2 = I, IP_
        DP_ = []; D2 = D, DP_
        MP_ = []; M2 = M, MP_
        Dy2 = Dy
        My2 = My
        G2 = G
        Olp2 = Olp
        Py_ = deque(P)

    else:  # continued PP vars are accumulated:

        L2, I2, D2, Dy2, M2, My2, G2, Olp2, Py_ = PP
        L2 += len(e_)
        I2 += I
        D2 += D; Dy2 += Dy
        M2 += M; My2 += My
        G2 += G
        Olp2 += Olp
        Py_.appendleft(P)

        # mx, mL, mI, mD, mM; ddx, dL, dI, dD, dM;
        # also norm pPs?

    form_pP(dx, dxP_); form_pP(mx, mxP_)  # same form_pP?
    form_pP(len(e_), LP_)
    form_pP(I, IP_)
    form_pP(D, DP_)
    form_pP(M, MP_)

    PP = s, L2, I2, D2, Dy2, M2, My2, G2, Olp2, Py_  # pP_s are packed in parameters

    # fork_ per blob, comp_P(_P, P_)?
    # rdn: alt_P, alt_PP, fork, alt_pP?
    return PP

''' 
        if typ: alt_oG *= ave_k; alt_oG = alt_oG.astype(int)  # ave V / I, to project V of odG
        else: oG *= ave_k; oG = oG.astype(int)               # same for h_der and h_comp eval?
        if oG > alt_oG:  # comp between overlapping vG and dG
            Olp += olp  # olp is assigned to the weaker of P | alt_P, == -> P: local access
        else:
            alt_P[7] += olp 
        olp, oG, alt_oG = 0,0,0
        P_, _P_, blob_ = scan_P_(typ, P, P_, _P_, blob_, x)  # no scan over hLe _Ps 
        olp += 1  # len of overlap to stronger alt-type P, accumulated until P or _P terminates
        G += g; alt_oG += alt_g  # for eval to assign olp to alt_rdn of vP or dP 
        after orient:
        if typ: e_.append((P, Pm, Pd))  # selective incremental-range comp_P, if min Py_? 
        else: e_.append(P)  # selective incremental-derivation comp_P?
'''


def form_pP(par, pP_):  # forming parameter patterns within PP

    # a_mx = 2; a_mw = 2; a_mI = 256; a_mD = 128; a_mM = 128: feedback to define vpPs: parameter value patterns
    # a_PM = a_mx + a_mw + a_mI + a_mD + a_mM  or A * n_vars, rdn accum per pP, alt eval per vertical overlap?

    # LIDV per dx, L, I, D, M? select per term?
    # alt2_: fork_ alt_ concat, to re-compute redundancy per PP


def term_PP(PP):  # vPP | dPP eval for blob | PP rotation, re-scan, re-comp, recursion, accumulation

    P2_ = []
''' 
    conversion of root to term, sum into wider fork, also sum per frame?
    
    dimensionally reduced axis: vP PP or contour: dP PP; dxP is direction pattern
    PP = PP, root_, blob_, _vPP_, _dPP_?
    vPP and dPP included in selected forks, rdn assign and form_PP eval after fork_ term in form_blob?
    blob= 0,0,0,0,0,0,0,0,0,0,[],[]  # crit, rdn, W, I2, D2, Dy2, M2, My2, G2, rdn2, alt2_, Py_
    vPP = 0,0,0,0,0,0,0,0,0,0,[],[]
    dPP = 0,0,0,0,0,0,0,0,0,0,[],[]  # P2s are initialized at non-matching P transfer to _P_?
    np.array for direct accumulation, or simply iterator of initialization?
    P2_ = np.array([blob, vPP, dPP],
        dtype=[('crit', 'i4'), ('rdn', 'i4'), ('W', 'i4'), ('I2', 'i4'), ('D2', 'i4'), ('Dy2', 'i4'),
        ('M2', 'i4'), ('My2', 'i4'), ('G2', 'i4'), ('rdn2', 'i4'), ('alt2_', list), ('Py_', list)]) 
    
    mean_dx = 1  # fractional?
    dx = Dx / H
    if dx > a: comp(abs(dx))  # or if dxP Dx: fixed ddx cost?  comp of same-sign dx only
    vx = mean_dx - dx  # normalized compression of distance: min. cost decrease, not min. benefit?
    
    
    eval of d,m adjust | _var adjust | x,y adjust if projected dS-, mS+ for min.1D Ps over max.2D
        if dw sign == ddx sign and min(dw, ddx) > a: _S /= cos (ddx)  # to angle-normalize S vars for comp
    if dw > a: div_comp (w): rw = w / _w, to width-normalize S vars for comp: 
        if rw > a: pn = I/w; dn = D/w; vn = V/w; 
            comp (_n) # or default norm for redun assign, but comp (S) if low rw?
            if d_n > a: div_comp (_n) -> r_n # or if d_n * rw > a: combined div_comp eval: ext, int co-variance?
        comp Dy and My, /=cos at PP term?  default div and overlap eval per PP? not per CP: sparse coverage?
        
    rrdn = 1 + rdn_w / len(e_)  # redundancy rate / w, -> P Sum value, orthogonal but predictive
    
    S = 1 if abs(D) + V + a * len(e_) > rrdn * aS else 0  # rep M = a*w, bi v!V, rdn I?
'''


def frame(f):  # postfix '_' denotes array vs. element, prefix '_' denotes higher-line variable

    global ave; ave = 127  # filters, ultimately set by separate feedback, then ave *= rng
    global rng; rng = 1

    global div_a; div_a = 127  # not justified
    global ave_k; ave_k = 0.25  # average V / I

    global Y; global X; Y, X = f.shape  # Y: frame height, X: frame width
    global y; y = 0

    _vP_, _dP_, frame_ = [], [], []

    t2_ = deque(maxlen=rng)  # vertical buffer of incomplete pixel tuples, for fuzzy ycomp
    t2__ = []  # vertical buffer + horizontal line: 2D array of 2D tuples, deque for speed?
    p_ = f[0, :]  # first line of pixels
    t_ = comp(p_)  # after part_comp (pop, no t_.append) while x < rng?

    for t in t_:
        p, d, m = t
        t2 = p, d, 0, m, 0  # fdy and fmy initialized at 0
        t2_.append(t2)  # only one tuple per first-line t2_
        t2__.append(t2_)  # in same order as t_

    # part_ycomp (pop, no form_P) while y < rng?

    for y in range(1, Y):  # vertical coordinate y is index of new line p_

        p_ = f[y, :]
        t_ = comp(p_)  # lateral pixel comparison
        t2__, _vP_, _dP_, vg_blob_, dg_blob_ = ycomp(t_, t2__, _vP_, _dP_) # vertical pixel comp

        frame_.append((vg_blob_, dg_blob_))  # line of blobs is added to frame of blobs

    return frame_  # frame of 2D patterns is outputted to level 2

f = misc.face(gray=True)  # input frame of pixels
f = f.astype(int)
frame(f)
