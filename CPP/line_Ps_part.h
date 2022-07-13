#pragma once

#include <vector>
#include "core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "core\core.hpp"
#include "../ConsoleApplication1/Utils/VariousUtils.h"  //timing

//#include <algorithm>
//#include <array>
//#include <math.h>

/* 12.July.2022, 13.7.2022
 C++ Port by Twenkid of Alex Pitertsev's port to Julia of Boris Kazachenko et al. Python level 1D algorithm from CogAlg:
 https://github.com/alex-pitertsev/CogAlg/blob/4ce0256d6a20a571d0e7243cdb0042d7037eee5d/line_1D_alg/julia/line_Ps_part.jl
 
 This is just a first quick draft done for 2 hours and a bit which is to show to the team how it looks like etc.
 PmPd are empty, why? - fixed, I've missed to add the creation of the new dert_;
 A few additional checks had to be added due to the lack of "nothing" in C++ in order to avoid attempt to acceess an empty P_ or something
 
 Tests with several tweaks with references and sending only the row of the Mat instead of slicing, but it seems the difference may be rather from different background loading.
 Best time of the main cycle traversing 768 lines of the raccoon.jpg is about 0.292 s, it varies to up to 0.3xx. 

*/ 

namespace Line1D {
    SYSTEMTIME mStartTime, mEndTime;
    CTiming timing; //Utils/VariousUtils.h    
   

    typedef cv::Point_<uint8_t> Pixel;

    /*
    CogAlg comments:

    line_Ps is a principal version of 1st - level 1D algorithm
    Operations :
    -
    -Cross - compare consecutive pixels within each row of image, forming dert_ : queue of derts, each a tuple of derivatives per pixel.
    dert_ is then segmented into patterns Pmsand Pds : contiguous sequences of pixels forming same - sign match or difference.
    Initial match is inverse deviation of variation : m = ave_ | d | -| d | ,
    rather than a minimum for directly defined match : albedo of an object doesn't correlate with its predictive value.
    -
    -Match patterns Pms are spans of inputs forming same - sign match.Positive Pms contain high - match pixels, which are likely
    to match more distant pixels.Thus, positive Pms are evaluated for cross - comp of pixels over incremented range.
    -
    -Difference patterns Pds are spans of inputs forming same - sign ds.d sign match is a precondition for d match, so only
    same - sign spans(Pds) are evaluated for cross - comp of constituent differences, which forms higher derivatives.
    (d match = min: rng + comp value : predictive value of difference is proportional to its magnitude, although inversely so)
    -
    Both extended cross - comp forks are recursive : resulting sub - patterns are evaluated for deeper cross - comp, same as top patterns.
    These forks here are exclusive per P to avoid redundancy, but they do overlap in line_patterns_olp.
    """
    */

    //Conversion from Alex code to C++, Todor, #12-7-2022, bug fixes and timing: 13-7-2022

    //using Images, ImageView, CSV

    //# instead of the class in Python version in Julia values are stored in the struct

    struct Cdert_ {
        int i;
        int p;
        int d;
        int m;
        bool mrdn;
        Cdert_(int i, int  p, int  d, int  m, bool  mrdn) : i{ i }, p{ p }, d{ d }, m{ m }, mrdn{ mrdn } {
        }
    };

    //dert_ = Cdert_[] # line - wide i_, p_, d_, m_, mrdn_  --> local

    //std::array<Cdert_, 3> dert_;
   // std::vector<Cdert_> dert_;

    struct CP {
        int L, I, D, M, Rdn, x0;
        //  M::Int32  # summed ave - abs(d), different from D
        // Rdn::Int32  # 1 + binary dert.mrdn cnt / len(dert_)       
        std::vector<Cdert_> dert_; //dert_::Vector{ Cdert_ }  # contains(i, p, d, m, mrdn)
        CP(int L, int I, int D, int M, int Rdn, int x0, std::vector<Cdert_> dert_) : L{ L }, I{ I }, D{ D }, M{ M }, Rdn{ Rdn }, x0{ x0 }, dert_{ dert_ } {
        }
        /*
         # subset = list  # 1st sublayer' rdn, rng, xsub_pmdertt_, _xsub_pddertt_, sub_Ppm_, sub_Ppd_
    # # for layer-parallel access and comp, ~ frequency domain, composition: 1st: dert_, 2nd: sub_P_[ dert_], 3rd: sublayers[ sub_P_[ dert_]]:
          # sublayers = list  # multiple layers of sub_P_s from d segmentation or extended comp, nested to depth = sub_[n]
            # subDertt_ = list  # m, d' [L,I,D,M] per sublayer, conditionally summed in line_PPs
            # derDertt_ = list  # for subDertt_s compared in line_PPs
            */
    };

    bool verbose = false;

    //# pattern filters or hyper - parameters: eventually from higher - level feedback, initialized here as constants :
    int ave = 15; // #  | difference | between pixels that coincides with average value of Pm
    int ave_min = 2;  // for m defined as min | d | : smaller ?
    int ave_M = 20;  // min M for initial incremental - range comparison(t_), higher cost than der_comp ?
    int ave_D = 5;  // min | D | for initial incremental - derivation comparison(d_)
    int ave_nP = 5;  // average number of sub_Ps in P, to estimate intra - costs ? ave_rdn_inc = 1 + 1 / ave_nP // 1.2
    float ave_rdm = 0.5;  // obsolete: average dm / m, to project bi_m = m * 1.5
    int ave_splice = 50;  // to merge a kernel of 3 adjacent Ps
    int init_y = 501;  // starting row, set 0 for the whole frame, mostly not needed
    int halt_y = 501;  // ending row, set 999999999 for arbitrary image.
    //#! these values will be one more (+1) in the Julia version because of the numbering specifics
    /* """
    Conventions:
    postfix 't' denotes tuple, multiple ts is a nested tuple
    postfix '_' denotes array name, vs.same - name elements
    prefix '_'  denotes prior of two same - name variables
    prefix 'f'  denotes flag
    1 - 3 letter names are normally scalars, except for Pand similar classes,
    capitalized variables are normally summed small - case variables,
    longer names are normally classes
    """
    */
    bool logging = true;

    //std::vector<CP> form_P_(void* rootP, std::vector<Cdert_> dert_, int rdn, int rng, int fPd)
    std::vector<CP> form_P_(void* rootP, std::vector<Cdert_> dert_, int rdn, int rng, int fPd)
    {  // rdn, rng, fPd) {
        //function form_P_(rootP, dert_; rdn, rng, fPd)  # after semicolon in Julia keyword args are placed
        //# initialization:
          //P_ = CP[]  # structure to store all the form_P(layer1) output
        
        auto P_ = std::vector<CP>(); // CP[]  # structure to store all the form_P(layer1) output
        int x = 0;
        bool _sign=false, sign=!_sign; // = nothing  # to initialize 1st P, (nothing != True) and (nothing != False) are both True
        bool b_first_run = true;

        for (auto const dert : dert_) {
            //# segment by sign
            if (fPd) sign = dert.d > 0;
            else
                sign = dert.m > 0;

            if (sign != _sign || b_first_run) { //first are set to be diff
                //# sign change, initializeand append P
                auto L = 1;
                auto I = dert.p;
                auto D = dert.d;
                auto M = dert.m;
                auto Rdn = dert.mrdn + 1;
                auto x0 = x;//# sublayers = [], # Rdn starts from 1
                    //push!(P_, CP(L, I, D, M, Rdn, x0, [dert])) // # save data in the struct
                auto new_dert_ = std::vector<Cdert_>({ dert });
                auto cp = CP(L, I, D, M, Rdn, x0, new_dert_); // std::vector<Cdert_>({ dert }
                P_.push_back(cp);
                b_first_run = false; 
            }
            else
            {//# accumulate params :  
                //if empty - error when trying to access P_.back()! it shoudn't be empty?
                if (P_.size() == 0) {
                    auto end = P_.front(); // .back();
                    end.L += 1;
                    end.I += dert.p;
                    end.D += dert.d;
                    end.M += dert.m;
                    end.Rdn += dert.mrdn;
                    //push!(P_[end].dert_, dert)
                    end.dert_.push_back(dert);
                }
                else {
                    auto end = P_.back();
                    end.L += 1;
                    end.I += dert.p;
                    end.D += dert.d;
                    end.M += dert.m;
                    end.Rdn += dert.mrdn;
                    //push!(P_[end].dert_, dert)
                    end.dert_.push_back(dert);
                }
            }
            x += 1;
            _sign = sign;
        } //for ...

        if (logging == 1)
            if (!fPd) {
                //CSV.write("./layer1_Pm_log_jl.csv", P_)
                //printf("log layer1_Pm_log....");
            }
            else
            {
                //CSV.write("./layer1_Pd_log_jl.csv", P_)
                //printf("Pd_log... {_");
            }
        return P_; // # used only if not rootP, else packed in rootP.sublayers
    }

   //returning a tuple (pair) instead of a vector, #13-7-2022
    //std::pair <std::vector<CP>, std::vector<CP>> line_Ps_root_pair(cv::Mat pixel_) {
    std::pair <std::vector<CP>, std::vector<CP>> line_Ps_root_pair(cv::Mat pixel_) { //only one line is sent
        //printf("line_Ps_root");
        std::vector<Cdert_> dert_; // //local  # line - wide i_, p_, d_, m_, mrdn_
        //auto_i = frame.at<Pixel>(r, c) = pixel;  //pixel_[0];
        auto _i = pixel_.at<Pixel>(0, 0).x; // = pixel;  //pixel_[0];    
        //auto shift = pixel_.colRange(1, pixel_.cols);

        for (auto c = 1; c < pixel_.cols; c++) {
            int i = pixel_.at<Pixel>(0, c).x; // = pixel_.at<Pixel>(0, c);  //pixel;
            int d = i - _i;
            int p = i + _i;
            int m = ave - abs(d);
            int mrdn = m < 0;
            //push!(dert_, Cdert_(i, p, d, m, mrdn))  # save data in the struct
            //_i = i      
            //auto new_dert_ = std::vector<Cdert_>({ dert });
            //auto cp = CP(L, I, D, M, Rdn, x0, new_dert_);
            dert_.push_back(Cdert_(i, p, d, m, mrdn));
        }
        //# form patterns, evaluate them for rng + and der + sub - recursion of cross_comp :
        auto Pm_ = form_P_(NULL, dert_, 1, 1, false); // rdn = 1, rng = 1, fPd = false)  # rootP = None, eval intra_P_(calls form_P_)
        auto Pd_ = form_P_(NULL, dert_, 1, 1, true); // rdn = 1, rng = 1, fPd = true); ==> to the global dert_

        if (logging == 1) {
            //printf("CSV.write(. / layer0_log_jl.csv, dert_)");
            //CSV.write("./layer0_log_jl.csv", dert_)
        }
        //std::vector<std::vector<CP>> PmPd = std::vector<std::vector<CP>>({ Pm_, Pd_ });
        std::pair <std::vector<CP>, std::vector<CP>> PmPd(Pm_, Pd_);
        return PmPd;
    }


    //returning a tuple (pair) instead of a vector, #13-7-2022
     //std::pair <std::vector<CP>, std::vector<CP>> line_Ps_root_pair(cv::Mat pixel_) {
    std::pair <std::vector<CP>, std::vector<CP>> line_Ps_root_pair_rowonly(cv::Mat& pixel_, int y) {
        //printf("line_Ps_root");
        std::vector<Cdert_> dert_; // //local  # line - wide i_, p_, d_, m_, mrdn_
        //auto_i = frame.at<Pixel>(r, c) = pixel;  //pixel_[0];
        auto _i = pixel_.at<Pixel>(y, 0).x; // = pixel;  //pixel_[0];    
        //auto shift = pixel_.colRange(1, pixel_.cols);

        for (auto c = 1; c < pixel_.cols; c++) {
            int i = pixel_.at<Pixel>(y, c).x; // = pixel_.at<Pixel>(0, c);  //pixel;
            int d = i - _i;
            int p = i + _i;
            int m = ave - abs(d);
            int mrdn = m < 0;
            //push!(dert_, Cdert_(i, p, d, m, mrdn))  # save data in the struct
            //_i = i      
            //auto new_dert_ = std::vector<Cdert_>({ dert });
            //auto cp = CP(L, I, D, M, Rdn, x0, new_dert_);
            dert_.push_back(Cdert_(i, p, d, m, mrdn));
        }
        //# form patterns, evaluate them for rng + and der + sub - recursion of cross_comp :
        auto Pm_ = form_P_(NULL, dert_, 1, 1, false); // rdn = 1, rng = 1, fPd = false)  # rootP = None, eval intra_P_(calls form_P_)
        auto Pd_ = form_P_(NULL, dert_, 1, 1, true); // rdn = 1, rng = 1, fPd = true); ==> to the global dert_

        if (logging == 1) {
            //printf("CSV.write(. / layer0_log_jl.csv, dert_)");
            //CSV.write("./layer0_log_jl.csv", dert_)
        }
        //std::vector<std::vector<CP>> PmPd = std::vector<std::vector<CP>>({ Pm_, Pd_ });
        std::pair <std::vector<CP>, std::vector<CP>> PmPd(Pm_, Pd_);
        return PmPd;
    }


    //std::vector<std::vector<CP>> line_Ps_root(cv::Mat pixel_) {
    std::vector<std::vector<CP>> line_Ps_root_vector(cv::Mat pixel_) {
        //printf("line_Ps_root");
        std::vector<Cdert_> dert_; // //local  # line - wide i_, p_, d_, m_, mrdn_
        //auto_i = frame.at<Pixel>(r, c) = pixel;  //pixel_[0];
        auto _i = pixel_.at<Pixel>(0, 0).x; // = pixel;  //pixel_[0];    
        //auto shift = pixel_.colRange(1, pixel_.cols);

        for (auto c = 1; c < pixel_.cols; c++) {
            int i = pixel_.at<Pixel>(0, c).x; // = pixel_.at<Pixel>(0, c);  //pixel;
            int d = i - _i;
            int p = i + _i;
            int m = ave - abs(d);
            int mrdn = m < 0;
            //push!(dert_, Cdert_(i, p, d, m, mrdn))  # save data in the struct
            //_i = i      
            //auto new_dert_ = std::vector<Cdert_>({ dert });
            //auto cp = CP(L, I, D, M, Rdn, x0, new_dert_);
            dert_.push_back(Cdert_(i, p, d, m, mrdn));
        }
        //# form patterns, evaluate them for rng + and der + sub - recursion of cross_comp :
        auto Pm_ = form_P_(NULL, dert_, 1, 1, false); // rdn = 1, rng = 1, fPd = false)  # rootP = None, eval intra_P_(calls form_P_)
        auto Pd_ = form_P_(NULL, dert_, 1, 1, true); // rdn = 1, rng = 1, fPd = true); ==> to the global dert_

        if (logging == 1) {
            //printf("CSV.write(. / layer0_log_jl.csv, dert_)");
            //CSV.write("./layer0_log_jl.csv", dert_)
        }
        std::vector<std::vector<CP>> PmPd = std::vector<std::vector<CP>>({ Pm_, Pd_ });        
        return PmPd;
        //return[Pm_, Pd_] // # input to 2nd level
    }

    /*
    function line_Ps_root(pixel_)  # Ps: patterns, converts frame_of_pixels to frame_of_patterns, each pattern may be nested
    local _i = pixel_[1]  #!differs from the python version # cross_comparison :
        for i in pixel_[2:end]  # pixel i is compared to prior pixel _i in a row :
    d = i - _i  # accum in rng
    p = i + _i  # accum in rng
    m = ave - abs(d)  # for consistency with deriv_comp output, else redundant
    mrdn = m < 0  # 1 if abs(d) is stronger than m, redundant here
        push!(dert_, Cdert_(i, p, d, m, mrdn))  # save data in the struct
        _i = i
        end

        # form patterns, evaluate them for rng + and der + sub - recursion of cross_comp :
    Pm_ = form_P_(nothing, dert_, rdn = 1, rng = 1, fPd = false)  # rootP = None, eval intra_P_(calls form_P_)
    Pd_ = form_P_(nothing, dert_, rdn = 1, rng = 1, fPd = true)

    if logging == 1
    CSV.write("./layer0_log_jl.csv", dert_)
     end

    return[Pm_, Pd_]  # input to 2nd level
    end
    */

    void Line1D(std::string path) {
        int render = 0;
        int fline_PPs = 0;
        int frecursive = 0;
        int logging = 1; // # logging of local functions variables
        //# logging = 0  # logging of local functions variables

        //# image_path = "./line_1D_alg/raccoon.jpg";
        std::string image_path = path; // "Z:\\r.jpg";  // "/home/alex/Python/CogAlg/line_1D_alg/raccoon.jpg";
        cv::Mat image; // // nothing

        //# check if image exist
            //if isfile(image_path)
       //image = cv::imread(image_path,1 ); //color
        image = cv::imread(image_path, 0); //grayscale -- no need to convert
        //image = load(image_path)  # read as N0f8(normed 0...1) type array of cells
        //else
        //println("ERROR: Image not found!")
        //end
        //gray_image = Gray.(image)  # convert rgb N0f8 to gray N0f8 array of cells
        //img_channel_view = channelview(gray_image)  # transform to the array of N0f8 numbers
        //gray_image_int = convert.(Int16, trunc.(img_channel_view .* 255))  # finally get array of 0...255 numbers

        cv::imshow("image", image);
        cv::waitKey();
        //image.shape
        unsigned X = image.cols;
        unsigned Y = image.rows;

        printf("\nX, Y = %d, %d", X, Y);
        //# Main
         //Y, X = size(gray_image) # Y: frame height, X : frame width

          //std::vector< > 
            //frame = []

         // # y is index of new row pixel_, we only need one row, use init_y = 1, halt_y = Y for full frame
        //init_y = 500;
        init_y = 0;
        halt_y = 1000;

        //unsigned y_end = halt_y < Y ? halt_y : Y; //why sets it to 1000, the bigger value?
        unsigned y_end = halt_y < Y ? halt_y : Y;

        cout << "halt_y=" << y_end;
        printf("\nTIMING for %d, %d", init_y, y_end); // halt_y);
        mStartTime = timing.now();        
        //for (int y = init_y; y < min(halt_y, Y); y++) {
        for (unsigned y = init_y; y < y_end; y++) {
            //printf("\ny=%u, y_end=%u", y, y_end);
            //line_Ps_root(gray_image_int[y, :])  # line = [Pm_, Pd_]
            //line_Ps_root(gray_image_int[y, :])  # line = [Pm_, Pd_]        
            //NO, col: x, row: y line_Ps_root(image.colRange(y, y+1).rowRange(0, X)); // gray_image_int[y, :])  # line = [Pm_, Pd_]
            
            //line_Ps_root_vector(image.rowRange(y, y+1).colRange(0, X)); // gray_image_int[y, :])  # line = [Pm_, Pd_]          
             line_Ps_root_pair(image.rowRange(y, y+1).colRange(0, X)); // gray_image_int[y, :])  # line = [Pm_, Pd_]
            //line_Ps_root_pair_rowonly(image, y); // .rowRange(y, y + 1).colRange(0, X)); //0.331, 0.330, 0.328. 0:0:0.335
            //Similar time, 0.309, 0.327, 0.360, 0.344, 0.362; line 0:0:0.332, 0:0:0.326, 0:0:0.365, 0:0:0.373
            //turning off some programs, pair: 0:0:0.293, 0:0:0.303, 0:0:0.307, 0:0:0.302s
            //y,y or y,y+1?
            //System Idle process is loading my computer, timing seems to vary more depending on some background load on my PC, rather than the tweaks with references,
            //Sending the original image Mat and the row only etc.
            //Min 0.292, Average about 0.3x s, depending on the system processes load on i5-6500 Win10            
            // Channel_16SC1(Range(0, height), Range(1, width)).copyTo(maskX(Range(0, height), Range(0, width - 1))); //X
        }
        mEndTime = timing.now();
        cout << endl << "line_Ps_root_pair" << endl;
        timing.timeDiff(mStartTime, mEndTime, 1); //print time       

        mStartTime = timing.now();
        //for (int y = init_y; y < min(halt_y, Y); y++) {
        for (unsigned y = init_y; y < y_end; y++) {
            
            line_Ps_root_vector(image.rowRange(y, y+1).colRange(0, X)); // gray_image_int[y, :])  # line = [Pm_, Pd_]                      
        }
        mEndTime = timing.now();
        cout << endl << "line_Ps_root_vector" << endl;
        timing.timeDiff(mStartTime, mEndTime, 1); //print time       

        mStartTime = timing.now();
        //for (int y = init_y; y < min(halt_y, Y); y++) {
        for (unsigned y = init_y; y < y_end; y++) {            
           line_Ps_root_pair_rowonly(image, y); // .rowRange(y, y + 1).colRange(0, X)); //0.331, 0.330, 0.328. 0:0:0.335
        }
        cout << endl << "line_Ps_root_pair_rowonly" << endl;
        mEndTime = timing.now();
        timing.timeDiff(mStartTime, mEndTime, 1); //print time               
    }

}
