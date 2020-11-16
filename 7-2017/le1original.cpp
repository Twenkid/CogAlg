from scipy import misc

/*'''
Level 1:
Cross-comparison between consecutive pixels within horizontal scan line (row).
Resulting difference patterns dPs: spans of pixels forming same-sign differences,
and relative match patterns vPs: spans of pixels forming same-sign predictive value,
overlap within each line of pixels.
'''
*/
/*
 C++ conversion by Todor Arnaudov
 
 Guidelines:
 
 Leading T_ means "type" or "reTurn Type" - for structs; return and param
 Trailing p_ in Boris convention means a list ("array"). For now STL vector is used (dynamic resizing, holds size). A simple representation (pure array as a pointer) may be also used, but then the maximum size should be allocated and the lenght should be send as an additional parameter, e.g.:
 
  inc_rng(..., p_, p_len)
  
  X = len(p_) ==> X = p_.size()
*/

//Initial tests with simple usage of float. Further optimisation may consider integer for most parts, until the divisions.

#include <misc.h>;


struct T_inc_rng{  //inc_rng
     vector<Item> vP_, dP_;
};

typedef T_inc_rng T_vP_dP_;

struct T_Le1{
 vector<> FP_;
};

T_inc_rng inc_rng(Item a,Item aV,Item aD,Item min_r,Item A,Item AV,Item AD,Item r, vector<Item> p_):


    if (r > min_r){  // A, AV, AD inc.to adjust for redundancy to patterns formed by prior comp:
        A += a;    // a: min m for inclusion into positive vP
        AV += aV;  // aV: min V for initial comp() recursion, AV: min V for higher recursions
    }

    if (r > (min_r-1)){  // default range is shorter for d_[w]: redundant ds are smaller than ps
        AD += aD;    // aV: min |D| for comp() recursion over d_[w], AD: min |D| for recursion
    }
    
    //keep a list of declarations, if not declared - declare ... 
    int X = p_.size(); //X = len(p_);
    
    //not declared, declare --> it's a list (p_) --> vector; pixels can be integer, but for now - float
    //copy constructor
    vector<Item> ip_(p_);
    //ip_ = p_  // to differentiate from new p_ // test C++!

    //declare vP_, dP_ 
    vector<Item> vP_, dP_;
    
    //vP_, dP_ = [],[]  // r was incremented in higher-scope p_
    
    Item pri_s=0, I=0, D=0, V=0, rv=0, olp=0;
    vector<Item> p_, olp_;
    
    #pri_s, I, D, V, rv, olp, p_, olp_ = 0, 0, 0, 0, 0, 0, [], []  // tuple vP=0
    #pri_sd, Id, Dd, Vd, rd, dolp, d_, dolp_ = 0, 0, 0, 0, 0, 0, [], []  // tuple dP=0
    
    Item pri_sd=0, Id=0, Dd=0, Vd=0, rd=0, dolp=0;
    vector<Item> d_, dolp_;
    

    //in Python range(a, b) is until < of b !
    for (x=r+1; x<X; x++) //range(r+1, X):
    {
        p, fd, fv = ip_[x]       // compared to a pixel at x-r-1:
        pp, pfd, pfv = ip_[x-r]  // previously compared p(ignored), its fd, fv to next p
        fv += pfv  // fv is summed over extended-comp range
        fd += pfd  // fd is summed over extended-comp range

        pri_p, pri_fd, pri_fv = ip_[x-r-1]  // for comp(p, pri_p), pri_fd and pri_fv ignored

        pri_s, I, D, V, rv, p_, olp, olp_, pri_sd, Id, Dd, Vd, rd, d_, dolp, dolp_, vP, dP_ = \
        comp(p, pri_p, fd, fv, x, X,
             pri_s, I, D, V, rv, p_, olp, olp_,
             pri_sd, Id, Dd, Vd, rd, d_, dolp, dolp_,
             a, aV, aD, min_r, A, AV, AD, r, vP_, dP_)
    }

    T_inc_rng vP_dP_(vP_, dP_);
    return vP_dP_;
    
    #return vP_, dP_  // local vPs and dPs to replace p_, A, AV, AD accumulated per comp recursion


def inc_der(a, aV, aD, min_r, A, AV, AD, r, d_):

    if (r > min_r){ A += a; AV += aV; }
    if (r > min_r-1) AD += aD;

    X = d_.size(); //len(d_)
    #ip_ = d_  // to differentiate from new d_
    vector<Item> ip_(d_);

/*     
    fd, fv, r, vP_, dP_ = 0, 0, 0, [], []  // r is initialized for each d_
    pri_s, I, D, V, rv, olp, p_, olp_ = 0, 0, 0, 0, 0, 0, [], []  // tuple vP=0,
    pri_sd, Id, Dd, Vd, rd, dolp, d_, dolp_ = 0, 0, 0, 0, 0, 0, [], []  // tuple dP=0
 */
 
    Item fd =0, fv =0, r = 0;
    vector<Item> vP_, dP; 
    
    Item pri_s=0, I=0, D=0, V=0, rv=0, olp;
    vector<Item> p_, olp_;
    
    Item pri_sd=0, Id=0, Dd=0, Vd=0, rd=0, dolp=0;
    vector<Item> d_, dolp_ ;
    
 
    pri_p = ip_[0]

    #for x in range(1, X):
    for(x=1; x<X; x++){

        p = ip_[x]  // better than pop()?

        pri_s, I, D, V, rv, p_, olp, olp_, pri_sd, Id, Dd, Vd, rd, d_, dolp, dolp_, vP, dP_ = \
        comp(p, pri_p, fd, fv, x, X,
             pri_s, I, D, V, rv, p_, olp, olp_,
             pri_sd, Id, Dd, Vd, rd, d_, dolp, dolp_,
             a, aV, aD, min_r, A, AV, AD, r, vP_, dP_)

        pri_p = p
    }

    return vP_, dP_  // local vPs and dPs to replace d_


def comp(Item p, Item pri_p, Item fd, Item fv, Item x, Item X,  // input variables

         Item pri_s, Item I, Item D, Item V, Item rv, vector<Item> p_, Item olp, vector<Item> olp_,  // variables of vP
		 
         Item pri_sd, Item Id, Item Dd, Item Vd, Item rd, vector<Item> d_, Item dolp, vector<Item> dolp_,  // variables of dP
		 
         Item a, Item aV, Item aD, Item min_r, Item A, Item AV, Item AD, Item r, vector<Item> vP_, vector<Item> dP_)  // filter variables and output patterns
{

    Item d = p - pri_p      // difference between consecutive pixels
    Item m = min(p, pri_p)  // match between consecutive pixels
    Item  v = m - A          // relative match (predictive value) between consecutive pixels

    Item fd += d  // fd includes all shorter + current- range ds between comparands
    Item fv += v  // fv includes all shorter + current- range vs between comparands


    // formation of value pattern vP: span of pixels forming same-sign v s:

    Item s;
    (v > 0) ? s=1 : s=0;  // s: positive sign of v
	
    //if x > r+2 and (s != pri_s or x == X-1):  // if derived pri_s miss, vP is terminated
	
	if ((x>r+2) && ((s!=pri_s) || (x==X-1)))
    {
        if len(p_) > r+3 and pri_s == 1 and V > AV:  // min 3 comp over extended distance within p_:

            r += 1  // r: incremental range-of-comp counter
            rv = 1  // rv: incremental range flag:
            p_.append(inc_rng(a, aV, aD, min_r, A, AV, AD, r, p_))

        p = I / len(p_); d = D / len(p_); v = V / len(p_)  // default to eval overlap, poss. div.comp?
        vP = pri_s, p, I, d, D, v, V, rv, p_, olp_
        vP_.append(vP)  // output of vP, related to dP_ by overlap only, no discont comp till Le3?

        o = len(vP_), olp  // len(P_) is index of current vP
        dolp_.append(o)  // indexes of overlapping vPs and olp are buffered at current dP

        I, D, V, rv, olp, dolp, p_, olp_ = 0, 0, 0, 0, 0, 0, [], []  // initialization of new vP and olp
	}

    pri_s = s   // vP (span of pixels forming same-sign v) is incremented:
    olp += 1    // overlap to current dP
    I += pri_p  // ps summed within vP
    D += fd     // fds summed within vP into fuzzy D
    V += fv     // fvs summed within vP into fuzzy V
    pri = pri_p, fd, fv
    p_.append(pri)  // buffered within vP for selective extended comp


    // formation of difference pattern dP: span of pixels forming same-sign d s:

    sd = 1 if d > 0 else 0  // sd: positive sign of d;
    if x > r+2 and (sd != pri_sd or x == X-1):  // if derived pri_sd miss, dP is terminated

        if len(d_) > 3 and abs(Dd) > AD:  // min 3 comp within d_:

            rd = 1  // rd: incremental derivation flag:
            d_.append(inc_der(a, aV, aD, min_r, A, AV, AD, r, d_))

        pd = Id / len(d_); dd = Dd / len(d_); vd = Vd / len(d_)  // so all olp Ps can be directly evaluated
        dP = pri_sd, pd, Id, dd, Dd, vd, Vd, rd, d_, dolp_
        dP_.append(dP)  // output of dP

        o = len(dP_), dolp  // len(P_) is index of current dP
        olp_.append(o)  // indexes of overlapping dPs and dolps are buffered at current vP

        Id, Dd, Vd, rd, olp, dolp, d_, dolp_ = 0, 0, 0, 0, 0, 0, [], []  // initialization of new dP and olp

    pri_sd = sd  // dP (span of pixels forming same-sign d) is incremented:
    dolp += 1    // overlap to current vP
    Id += pri_p  // ps summed within dP
    Dd += fd     // fds summed within dP
    Vd += fv     // fvs summed within dP
    d_.append(fd)  // prior fds are buffered within dP, all of the same sign

    return pri_s, I, D, V, rv, p_, olp, olp_, pri_sd, Id, Dd, Vd, rd, d_, dolp, dolp_, vP_, dP_
    // for next p comparison, vP and dP increment, and output
}


#def Le1(Fp_): // last '_' distinguishes array name from element name
Ret_Le1 Le1(Fp_): // last '_' distinguishes array name from element name
    
    FP_ = vector<LP>; //lines of patterns
    
    #FP_ = []  // output frame of vPs: relative match patterns, and dPs: difference patterns
    
    unsigned w,h;
    #Y, X = Fp_.shape  // Y: frame height, X: frame width
    Fp_.shape(&w, &h);

    int a = 127  // minimal filter for vP inclusion
    int aV = 63  // minimal filter for incremental-range comp
    int aD = 63  // minimal filter for incremental-derivation comp
    int min_r=0  // default range of fuzzy comparison, initially 0

    #for y in range(Y):
    int y;
    for(y=0; y<Y; y++){

        ip_ = Fp_[y, :]  // y is index of new line ip_
        //take a line

        if min_r == 0: A = a; AV = aV  // actual filters, incremented per comp recursion
        else: A = 0; AV = 0  // if r > min_r

        if min_r <= 1: AD = aD
        else: AD = 0

        fd, fv, r, x, vP_, dP_ = 0, 0, 0, 0, [], []  // i/o tuple
        pri_s, I, D, V, rv, olp, p_, olp_ = 0, 0, 0, 0, 0, 0, [], []  // vP tuple
        pri_sd, Id, Dd, Vd, rd, dolp, d_, dolp_ = 0, 0, 0, 0, 0, 0, [], []  // dP tuple

        pri_p = ip_[0]

        for x in range(1, X):  // cross-compares consecutive pixels

            p = ip_[x]  // new pixel for comp to prior pixel, could use pop()?

            pri_s, I, D, V, rv, p_, olp, olp_, pri_sd, Id, Dd, Vd, rd, d_, dolp, dolp_, vP_, dP_ = \
            comp(p, pri_p, fd, fv, x, X,
                 pri_s, I, D, V, rv, p_, olp, olp_,
                 pri_sd, Id, Dd, Vd, rd, d_, dolp, dolp_,
                 a, aV, aD, min_r, A, AV, AD, r, vP_, dP_)

            pri_p = p  // prior pixel, pri_ values are always derived before use

        LP_ = vP_, dP_
        FP_.append(LP_)  // line of patterns is added to frame of patterns, y = len(FP_)
    }

    return FP_  // output to level 2

int main(){

  Misc misc;
  bool gray = true;
  #vector<Item> 
  Frame frame = misc.faceFloat(gray);
  T_Le1 ret_le1 = Le1(f);
  
  //f = misc.face(gray=True)  // input frame of pixels
  //f = f.astype(int)
  //Le1(f)
}

// at vP term: print ('type', 0, 'pri_s', pri_s, 'I', I, 'D', D, 'V', V, 'rv', rv, 'p_', p_)
// at dP term: print ('type', 1, 'pri_sd', pri_sd, 'Id', Id, 'Dd', Dd, 'Vd', Vd, 'rd', rd, 'd_', d_)