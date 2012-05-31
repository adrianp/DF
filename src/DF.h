/* 
 * This is an implementation of the DF1 and DF2 dynamic optimization problem 
 * generators presented in the paper: A Test Problem Generator for 
 * Non-Stationary Environments, Morrison, R., W., De Jong, K., A, 1999 IEEE.
 * 
 * This is implemented as an experiment for the AEON framework.
 * 
 * Author: Adrian Tudor Panescu, adrian.panescu@epfl.ch
 */

#ifndef DF_H_
#define DF_H_

// Lib includes
#include <fstream>
#include <vector>

// Boost includes
#include <boost/shared_ptr.hpp>

// AEON includes
#include "experiment/EvolutionaryTimeType.h"
#include "experiment/IParametricExperiment.h"
#include "experiment/result/Result.h"
#include "solution/ParametricSolution.h"

class DF: public virtual aeon::IParametricExperiment {

private:
	/**
	 * The heights of the peaks
	 */
	std::vector<double> H;

	/**
	 * The slopes of the peaks
	 */
	std::vector<double> R;

	/**
	 * The abscissae of the peaks
	 */
	std::vector<double> X;

	/**
	 * The ordinates of the peaks
	 */
	std::vector<double> Y;

	/**
	 * The total number of peaks
	 */
	unsigned int N;

	/**
	 * The number of local optima
	 */
	unsigned int local_optima;

	/**
	 * The lowest peak height
	 */
	double HBase;

	/*
	 * The highest peak height
	 */
	double HLimit;

	/**
	 * The lowest peak slope
	 */
	double RBase;

	/**
	 * The highest peak slope
	 */
	double RLimit;

	/*
	 * If set to true, the DF1 generator is used (cones); else, the DF2
	 * generator is used (hills)
	 */
	bool cones;

	/**
	 * The lowest value of the abscissae
	 */
	double XBase;

	/**
	 * The highest value of the abscissae
	 */
	double XLimit;

	/**
	 * The lowest value for the ordinates
	 */
	double YBase;

	/**
	 * The highest value for the ordinates
	 */
	double YLimit;

	/**
	 * The minimum value of the optimzation function
	 */
	double bottom;

	/**
	 * The number of dynamic peaks
	 */
	unsigned int dynamic_peaks;

	/**
	 * The IDs of the dynamic peaks
	 */
	std::vector<int> peaks_to_change;

	/**
	 * The IDs of the peaks that are changing slowly
	 */
	std::vector<int> slow_peaks_to_change;

	/**
	 * The IDs of the peaks that are changing faster
	 */
	std::vector<int> fast_peaks_to_change;

	/**
	 * The A value for changing the heights of slow peaks
	 */
	double Ah_slow;

	/**
	 * The A value for changing the heights of fast peaks
	 */
	double Ah_fast;

	/**
	 * The A value for changing the slopes of slow peaks
	 */
	double Ar_slow;

	/**
	 * The A value for changing the slopes of fast peaks
	 */
	double Ar_fast;

	/**
	 * The A value for changing the abscissae of slow peaks
	 */
	double Ax_slow;

	/**
	 * The A value for changing the abscissae of fast peaks
	 */
	double Ax_fast;

	/**
	 * The A value for changing the ordinates of slow peaks
	 */
	double Ay_slow;

	/**
	 * The A value for changing the ordinates of fast peaks
	 */
	double Ay_fast;


	/**
	 * The scaling factor value for changing the heights of slow peaks
	 */
	double scaleH_slow;

	/**
	 * The scaling factor value for changing the heights of fast peaks
	 */
	double scaleH_fast;

	/**
	 * The scaling factor value for changing the slopes of slow peaks
	 */
	double scaleR_slow;

	/**
	 * The scaling factor value for changing the slopes of fast peaks
	 */
	double scaleR_fast;

	/**
	 * The scaling factor value for changing the abscissae of slow peaks
	 */
	double scaleX_slow;

	/**
	 * The scaling factor value for changing the abscissae of fast peaks
	 */
	double scaleX_fast;

	/**
	 * The scaling factor value for changing the ordinates of slow peaks
	 */
	double scaleY_slow;

	/**
	 * The scaling factor value for changing the ordinates of fast peaks
	 */
	double scaleY_fast;


	/**
	 * The old values of Y (logistics function) for the height changes
	 */
	std::vector<double> Yheight;

	/**
	 * The old values of Y (logistics function) for the slope changes
	 */
	std::vector<double> Yslope;

	/**
	 * The old value of Y (logistics function) for the abscissae changes
	 */
	std::vector<double> Yx;

	/**
	 * The old value of Y (logistics function) for the ordinate changes
	 */
	std::vector<double> Yy;


	/**
	 * The current direction of height changes
	 */
	std::vector<int> Hflag;

	/**
	 * The current direction of slope changes
	 */
	std::vector<int> Rflag;

	/**
	 * The current direction of abscissae changes
	 */
	std::vector<int> Xflag;

	/**
	 * The current direction of ordinate changes
	 */
	std::vector<int> Yflag;

	/**
	 * The type of dynamic change
	 */
	std::string dynamic_type;

	/**
	 * Flag for outputting the landscape (linear interpolation)
	 */
	bool debug;

	/**
	 * The frequency of the change
	 */
	int frequency;

	/**
	 * The result of the evaluation [AEON]
	 */
	double result_;

	/**
	 * The number of parameters (dimensions) [AEON]
	 */
	unsigned int nParams_;

	/**
	 * The current evolutionary time
	 */
	int generation;

	/**
	 * Marks the covered peaks
	 */
	std::vector<bool> covered_peaks;

	/**
	 * The number of current solutions on peaks
	 */
	double solutions_on_peaks;

	/**
	 * The threshold used for considering a solution as beeing on a peak
	 */
	double peak_threshold;

	/**
	 * The current population size (as seen by the experiment)
	 */
	int pop_size;

	/**
	 * The file used for writing the statistics
	 */
	std::ofstream stats_file;

	/**
	 * Calculate the value of the optimization function
	 */
	double evaluateFunction(const std::vector<double>& x);

	/**
	 * Dynamic change: heights
	 */
	void changeHeights(double A, double scale);

	/**
	 * Dynamic change: slopes
	 */
	void changeSlopes(double A, double scale);

	/**
	 * Dynamic change: abscissae
	 */
	void changeX(double A, double scale);

	/**
	 * Dynamic change: ordinates
	 */
	void changeY(double A, double scale);

	/**
	 * Dynamic change: slow rotation
	 */
	void slow_rotation();

	/**
	 * Dynamic change: fast rotation
	 */
	void fast_rotation();

	/**
	 * Dynamic change: hybrid rotation
	 */
	void hybrid_rotation();

	/**
	 * Scale ALL the peaks (?)
	 */
	void dynamic_scaling();

	/**
	 * Changes to local optima, globals remain the same (?)
	 */
	void dynamic_local_optima();

	/**
	 * Local optima become global optima (?)
	 */
	void dynamic_elevate_optima();

	/**
	 * Changes to global optima, globals remain the same (?)
	 */
	void dynamic_global_optima();

	/*
	 * Global optima become local optima (?)
	 */
	void dynamic_downgrade_optima();

	/**
	 * Output the landscape using linear interpolation
	 */
	void reportLandscape(int generation);

	/**
	 * Print the parameters of the first changing peak
	 */
	void printFirstPeak();

	/**
	 * Report the statistics calculated by the benchmark
	 */
	void report_statistics();





public:

	/**
	 * Constructor
	 */
	DF();

	/**
	 * Destructor
	 */
	virtual ~DF();

	/**
	 * Function called by AEON for evaluating individuals [AEON]
	 */
	void evaluate(
			std::vector<boost::shared_ptr<const aeon::ParametricSolution> >& solution);

	/**
	 * Get the parameters used by the optimization function [AEON]
	 */
	void getParameters(std::vector<std::string>& parameters);

	/**
	 * Initialize the experiment [AEON]
	 */
	bool init(unsigned int seed, const std::string& argv,
			const std::string& mainFolder);

	/**
	 * Get the evaluation results [AEON]
	 */
	bool getResults(
			std::vector<boost::shared_ptr<const aeon::Result> >& results);

	/**
	 * Called by AEON before evaluating an individual [AEON]
	 */
	void setup(int evolutionaryTime,
			aeon::EvolutionaryTimeType::Value evolutionaryTimeType,
			int currentEvaluation, int idIndividual);

	/**
	 * Is this a dynamic experiment? [AEON]
	 */
	bool isDynamic();

};
#endif // DF_H_
