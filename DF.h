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

        unsigned int dynamic_peaks; // number of dynamic (moving) peaks


        std::vector<int> peaks_to_change; // the peaks that are moving
        // The A constants used by the logistics function
        double Ah;
        double Ar;
        double Ax;
        double Ay;
        // The scaling factors used when changing the landscape
        double scaleH;
        double scaleR;
        double scaleX;
        double scaleY;
        
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
        
        // the frequency of the change
        int frequency;

        // variables used by the AEON framework
        double result_;
        unsigned int nParams_;
        
        // used to keep track of generation number inside the experiment
        int generation; 

        double evaluateFunction(const std::vector<double>& x);
        
        // types of dynamic changes
        void changeHeights();
        void changeSlopes();
        void changeX();
        void changeY();
        void dynamic_rotation();
        void dynamic_scaling();
        void dynamic_local_optima();
        void dynamic_elevate_optima();
        void dynamic_global_optima();
        void dynamic_downgrade_optima();

        // the next functions are used only for debugging purposes
        double getFirstChanged();
        void printOptimum();
        void reportLandscape(int generation);

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