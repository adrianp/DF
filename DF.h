/* 
 * This is an implementation of the DF1 and DF2 dynamic optimization problem 
 * generators presented in the paper: A Test Problem Generator for 
 * Non-Stationary Environments, Morrison, R., W., De Jong, K., A, 1999 IEEE.
 * 
 * This is implemented as an experiment for the AEON framework.
 * 
 * Author: Adrian Tudor Panescu, adrian.panescu@epfl.ch
*/

// TODO: use conventions for ifndef ...
#ifndef DF_H_
#define DF_H_

// Lib includes
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <vector>

// AEON includes
#include "experiment/IParametricExperiment.h"
#include "experiment/result/Result.h"
#include "solution/ParametricSolution.h"

class DF : public virtual aeon::IParametricExperiment {

    private:
        std::vector<double> H; // height
        std::vector<double> R; // slope
        std::vector<double> X; // peaks' abscissae
        std::vector<double> Y; // peaks' ordinates
        unsigned int N; // total number of peaks
        unsigned int local_optimas; // number of local optimas
        double HBase; // minimum value for peak heights
        double HLimit; // maximum value for peak heights
        double RBase; // minimum value for peak slopes
        double RLimit; // maximum value for peak heights
        bool cones; // true if the landscape consists of cones, false for hills

        double XBase; // minimum value for peak abscissa
        double XLimit; // // maximum value for peak abscissa
        double YBase; // minimum value for peak ordinate
        double YLimit; // maximum value for peak ordinate
        
        double bottom; // the minimum of the function

        unsigned int dynamic_peaks; // number of dynamic (moving) peaks

        std::vector<int> peaks_to_change; // the peaks that are moving
        
        // the following vectors are used for hybrid dynamics (some slow peaks
        // and some fast ones)
        std::vector<int> slow_peaks_to_change;
        std::vector<int> fast_peaks_to_change;
        
        // The A constants used by the logistics function
        double Ah_slow;
        double Ah_fast;
        double Ar_slow;
        double Ar_fast;
        double Ax_slow;
        double Ax_fast;
        double Ay_slow;
        double Ay_fast;
        
        // The scaling factors used when changing the landscape
        double scaleH_slow;
        double scaleH_fast;
        double scaleR_slow;
        double scaleR_fast;
        double scaleX_slow;
        double scaleX_fast;
        double scaleY_slow;
        double scaleY_fast;
        
        // the old values of Y in the logistics function
        std::vector<double> Yheight;
        std::vector<double> Yslope;
        std::vector<double> Yx;
        std::vector<double> Yy;
        
        // the sign of the changes
        std::vector<int> Hflag;
        std::vector<int> Rflag;
        std::vector<int> Xflag;
        std::vector<int> Yflag;

        std::string dynamic_type; // the type of dynamic transformation
        bool debug; // true for outputting the landscape
        
        
        int frequency; // the frequency of the change

        // variables used by the AEON framework
        double result_;
        unsigned int nParams_;
        
        // used to keep track of generation number inside the experiment
        int generation; 
        
        double evaluateFunction(const std::vector<double>& x);
        
        // types of dynamic changes
        void changeHeights(double A, double scale);
        void changeSlopes(double A, double scale);
        void changeX(double A, double scale);
        void changeY(double A, double scale);
        void slow_rotation();
        void fast_rotation();
        void hybrid_rotation();
        void dynamic_scaling();
        void dynamic_local_optima();
        void dynamic_elevate_optima();
        void dynamic_global_optima();
        void dynamic_downgrade_optima();

        // the next functions are used only for debugging purposes
        //void getFirstChanged();
        void printOptimum();
        void reportLandscape(int generation);
        void printPosition();
        
        // elements used for reporting statistics
        void report_statistics();
        std::vector<bool> covered_peaks;
        double solutions_on_peaks;
        int pop_size;
        std::ofstream stats_file;
        double best; // unused
        
    public:
        DF();
        virtual ~DF();
        
        void evaluate(std::vector<boost::shared_ptr
                                  <const aeon::ParametricSolution> >& solution);
        
        void getParameters(std::vector<std::string>& parameters);
        
        bool init(unsigned int seed, 
                  const std::string& argv, 
                  const std::string& mainFolder);
        
        bool getResults(std::vector<boost::shared_ptr<const aeon::Result> >& 
                                                                       results);
        void setup(int evolutionaryTime, 
                   aeon::EvolutionaryTimeType::Value evolutionaryTimeType, 
                   int currentEvaluation);
        
        bool isDynamic();
};
#endif // DF_H_
