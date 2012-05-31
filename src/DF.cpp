/* 
 * This is an implementation of the DF1 and DF2 dynamic optimization problem 
 * generators presented in the paper: A Test Problem Generator for 
 * Non-Stationary Environments, Morrison, R., W., De Jong, K., A, 1999 IEEE.
 * 
 * This is implemented as an experiment for the AEON framework.
 * 
 * Author: Adrian Tudor Panescu, adrian.panescu@epfl.ch
 */

// Lib includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// Boost includes
#include <boost/assign/std/vector.hpp> // used for vector += element1, element2;
#include <boost/shared_ptr.hpp>

// AEON includes
#include "solution/ParametricSolution.h"

// AEON util includes
#include "json/json.h"
#include "utils/ConfFileParser.h"
#include "utils/Filesystem.h"

// Dynamic experiment includes
#include "DF.h"
#include "dynamic_helpers.h"
#include "random_helpers.h"

using namespace aeon;
using namespace boost::assign;
// used for vector += element1, element2;
using namespace std;

DF::DF() {
	this->generation = 0;
}

DF::~DF() {
	stats_file.close();
}

void DF::evaluate(
		std::vector<boost::shared_ptr<const aeon::ParametricSolution> >& solution) {
	this->result_ = this->evaluateFunction(solution[0]->getParameters());
}

bool DF::getResults(std::vector<boost::shared_ptr<const Result> >& results) {
	std::vector<double> v;
	v.push_back(this->result_); //hack for ResultData compliance
	results.push_back(
			boost::shared_ptr<const aeon::Result>(
					new aeon::Result(
							boost::shared_ptr<const aeon::ResultData>(
									new aeon::ResultData(this->result_, v)))));
	return true;
}

void DF::setup(int evolutionaryTime,
		EvolutionaryTimeType::Value /*evolutionaryTimeType*/,
		int /*currentEvaluation*/, int /*idIndividual*/) {

	if (evolutionaryTime > this->generation) {
		if (this->dynamic_peaks > 0
				&& evolutionaryTime % this->frequency == 0) {

			if (this->dynamic_type == "slow_rotation") {
				this->slow_rotation();
			} else if (this->dynamic_type == "fast_rotation") {
				this->fast_rotation();
			} else if (this->dynamic_type == "hybrid_rotation") {
				this->hybrid_rotation();
			} else if (this->dynamic_type == "scaling") {
				this->dynamic_scaling();
			} else if (this->dynamic_type == "local_optima") {
				this->dynamic_local_optima();
			} else if (this->dynamic_type == "elevate_optima") {
				this->dynamic_elevate_optima();
			} else if (this->dynamic_type == "global_optima") {
				this->dynamic_global_optima();
			} else if (this->dynamic_type == "downgrade_optima") {
				this->dynamic_downgrade_optima();
			}

			if (this->debug) {
				this->reportLandscape(evolutionaryTime); // this will make things slow
			}

		}

		this->report_statistics();
		this->solutions_on_peaks = 0;
		this->pop_size = 0;
		this->generation = evolutionaryTime;

		for (unsigned int i = 0; i < this->N - this->local_optima; ++i) {
			this->covered_peaks[i] = false;
		}

	}
}

void DF::report_statistics() {
	int aux = 0;
	for (unsigned int i = 0; i < this->N - this->local_optima; ++i) {
		if (this->covered_peaks[i]) {
			aux += 1;
		}
	}
	stats_file << this->generation << "\t" << this->HLimit << "\t"
			<< this->bottom << "\t" << this->frequency << "\t"
			<< this->solutions_on_peaks << "\t" << this->pop_size << "\t" << aux
			<< "\t" << (this->N - this->local_optima) << std::endl;
}

void DF::getParameters(std::vector<std::string>& parameters) {
	for (unsigned int i = 0; i < this->nParams_; ++i) {
		parameters.push_back("x" + i);
	}
}

bool DF::isDynamic() {
	if (this->dynamic_peaks > 0) {
		return true;
	} else {
		return false;
	}
}

bool DF::init(unsigned int, const std::string&, const std::string&) {
	unsigned int i;
	Json::Value rootNode;
	std::string value;

	// we need to initialize the ConfFileParser
	std::map<std::string, std::string> confFileVariables;
	utils::ConfFileParser::init(confFileVariables);

	if (!utils::ConfFileParser::parse("./user_conf/dynamic_config.json",
			rootNode, true)) {
		std::cout << "The dynamic experiment configuration file could not be "
				"parsed." << std::endl;
		return false;
	} else {

		if (!rootNode["N"]) {
			std::cout << "Value for N could not be parsed." << std::endl;
			return false;
		} else {
			this->N = rootNode["N"].asInt();
		}

		if (!rootNode["local_optima"]) {
			std::cout << "Value for local_optima could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->local_optima = rootNode["local_optima"].asInt();
		}

		if (!rootNode["HBase"]) {
			std::cout << "Value for HBase could not be parsed." << std::endl;
			return false;
		} else {
			this->HBase = rootNode["HBase"].asDouble();
		}
		if (!rootNode["HLimit"]) {
			std::cout << "Value for HLimit could not be parsed." << std::endl;
			return false;
		} else {
			this->HLimit = rootNode["HLimit"].asDouble();
		}

		if (!rootNode["bottom"]) {
			std::cout << "Value for bottom could not be parsed." << std::endl;
			return false;
		} else {
			this->bottom = rootNode["bottom"].asDouble();
		}

		if (!rootNode["RBase"]) {
			std::cout << "Value for RBase could not be parsed." << std::endl;
			return false;
		} else {
			this->RBase = rootNode["RBase"].asDouble();
		}
		if (!rootNode["RLimit"]) {
			std::cout << "Value for RLimit could not be parsed." << std::endl;
			return false;
		} else {
			this->RLimit = rootNode["RLimit"].asDouble();
		}

		if (!rootNode["XBase"]) {
			std::cout << "Value for XBase could not be parsed." << std::endl;
			return false;
		} else {
			this->XBase = rootNode["XBase"].asDouble();
		}
		if (!rootNode["XLimit"]) {
			std::cout << "Value for XLimit could not be parsed." << std::endl;
			return false;
		} else {
			this->XLimit = rootNode["XLimit"].asDouble();
		}

		if (!rootNode["YBase"]) {
			std::cout << "Value for YBase could not be parsed." << std::endl;
			return false;
		} else {
			this->YBase = rootNode["YBase"].asDouble();
		}
		if (!rootNode["YLimit"]) {
			std::cout << "Value for YLimit could not be parsed." << std::endl;
			return false;
		} else {
			this->YLimit = rootNode["YLimit"].asDouble();
		}

		if (!rootNode["dynamic_peaks"]) {
			std::cout << "Value for dynamic_peaks could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->dynamic_peaks = rootNode["dynamic_peaks"].asInt();
		}

		if (!rootNode["cones"]) {
			std::cout << "Value for cones could not be parsed." << std::endl;
			return false;
		} else {
			this->cones = rootNode["cones"].asBool();
		}

		if (!rootNode["Ah_slow"]) {
			std::cout << "Value for Ah_slow could not be parsed." << std::endl;
			return false;
		} else {
			this->Ah_slow = rootNode["Ah_slow"].asDouble();
		}
		if (!rootNode["Ah_fast"]) {
			std::cout << "Value for Ah_fast could not be parsed." << std::endl;
			return false;
		} else {
			this->Ah_fast = rootNode["Ah_fast"].asDouble();
		}

		if (!rootNode["Ar_slow"]) {
			std::cout << "Value for Ar_slow could not be parsed." << std::endl;
			return false;
		} else {
			this->Ar_slow = rootNode["Ar_slow"].asDouble();
		}
		if (!rootNode["Ar_fast"]) {
			std::cout << "Value for Ar_fast could not be parsed." << std::endl;
			return false;
		} else {
			this->Ar_fast = rootNode["Ar_fast"].asDouble();
		}

		if (!rootNode["Ax_slow"]) {
			std::cout << "Value for Ax_slow could not be parsed." << std::endl;
			return false;
		} else {
			this->Ax_slow = rootNode["Ax_slow"].asDouble();
		}
		if (!rootNode["Ax_fast"]) {
			std::cout << "Value for Ax_fast could not be parsed." << std::endl;
			return false;
		} else {
			this->Ax_fast = rootNode["Ax_fast"].asDouble();
		}

		if (!rootNode["Ay_slow"]) {
			std::cout << "Value for Ay_slow could not be parsed." << std::endl;
			return false;
		} else {
			this->Ay_slow = rootNode["Ay_slow"].asDouble();
		}
		if (!rootNode["Ay_fast"]) {
			std::cout << "Value for Ay_fast could not be parsed." << std::endl;
			return false;
		} else {
			this->Ay_fast = rootNode["Ay_fast"].asDouble();
		}

		if (!rootNode["scaleH_slow"]) {
			std::cout << "Value for scaleH_slow could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleH_slow = rootNode["scaleH_slow"].asDouble();
		}
		if (!rootNode["scaleH_fast"]) {
			std::cout << "Value for scaleH_fast could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleH_fast = rootNode["scaleH_fast"].asDouble();
		}

		if (!rootNode["scaleR_slow"]) {
			std::cout << "Value for scaleR_slow could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleR_slow = rootNode["scaleR_slow"].asDouble();
		}
		if (!rootNode["scaleR_fast"]) {
			std::cout << "Value for scaleR_fast could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleR_fast = rootNode["scaleR_fast"].asDouble();
		}

		if (!rootNode["scaleX_slow"]) {
			std::cout << "Value for scaleX_slow could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleX_slow = rootNode["scaleX_slow"].asDouble();
		}
		if (!rootNode["scaleX_fast"]) {
			std::cout << "Value for scaleX_fast could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleX_fast = rootNode["scaleX_fast"].asDouble();
		}

		if (!rootNode["scaleY_slow"]) {
			std::cout << "Value for scaleY_slow could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleY_slow = rootNode["scaleY_slow"].asDouble();
		}
		if (!rootNode["scaleY_fast"]) {
			std::cout << "Value for scaleY_fast could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->scaleY_fast = rootNode["scaleY_fast"].asDouble();
		}

		if (!rootNode["debug"]) {
			std::cout << "Value for debug could not be parsed." << std::endl;
			return false;
		} else {
			this->debug = rootNode["debug"].asBool();
		}

		if (!rootNode["dynamic_type"]) {
			std::cout << "Value for dynamic_type could not be parsed."
					<< std::endl;
			return false;
		} else {
			std::vector<std::string> dynamic_types;
			dynamic_types += "slow_rotation", "fast_rotation", "hybrid_rotation", "scaling", "local_optima", "elevate_optima", "global_optima", "downgrade_optima";
			if (std::find(dynamic_types.begin(), dynamic_types.end(),
					rootNode["dynamic_type"].asString())
					!= dynamic_types.end()) {
				this->dynamic_type = rootNode["dynamic_type"].asString();
			} else {
				std::cout << "Specified dynamic change type could not be found."
						<< std::endl;
				return false;
			}
		}

		if (!rootNode["frequency"]) {
			std::cout << "Value for frequency could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->frequency = rootNode["frequency"].asInt();
		}

		if (!rootNode["seed"]) {
			std::cout << "Value for seed could not be parsed." << std::endl;
			return false;
		} else {
			srand(rootNode["seed"].asUInt());
		}

		if (!rootNode["peak_threshold"]) {
			std::cout << "Value for peak_threshold could not be parsed."
					<< std::endl;
			return false;
		} else {
			this->peak_threshold = rootNode["peak_threshold"].asDouble();
		}
	}

	// TODO: make the experiment work in more than 2 dimensions
	this->nParams_ = 2;

	if (this->local_optima > this->N || this->dynamic_peaks > this->N) {
		std::cout
				<< "BENCHMARK ERROR: N must agree with local_optima and dynamic_peaks"
				<< std::endl;
		return false;
	}

	// the global optima value will be equal to HLimit
	// TODO: add parameter for global optima value
	double globalOpt = this->HLimit;

	for (i = 0; i < this->N; ++i) {
		if (i < this->N - this->local_optima) {
			this->H.push_back(globalOpt);
			this->covered_peaks.push_back(false);
		} else {
			this->H.push_back(random_double_exclusive(this->HBase, globalOpt));
		}
		this->R.push_back(random_double(this->RBase, this->RLimit));
		this->X.push_back(random_double(this->XBase, this->XLimit));
		this->Y.push_back(random_double(this->YBase, this->YLimit));
	}

	// in hybrid rotation half of the peaks will move slow and half fast
	if (this->dynamic_type == "hybrid_rotation") {
		for (i = 0; i < this->N; ++i) {
			if (i % 2 == 0) {
				this->fast_peaks_to_change.push_back(i);
			} else {
				this->slow_peaks_to_change.push_back(i);
			}

			this->Yheight.push_back(random_double(0, 1));
			this->Yslope.push_back(random_double(0, 1));
			this->Yx.push_back(random_double(0, 1));
			this->Yy.push_back(random_double(0, 1));

			// normalization: flag must be either 1 or -1
			this->Hflag.push_back(rand() % 2);
			if (this->Hflag[i] == 0) {
				this->Hflag[i] = -1;
			}

			this->Rflag.push_back(rand() % 2);
			if (this->Rflag[i] == 0) {
				this->Rflag[i] = -1;
			}

			this->Xflag.push_back(rand() % 2);
			if (this->Xflag[i] == 0) {
				this->Xflag[i] = -1;
			}

			this->Yflag.push_back(rand() % 2);
			if (this->Yflag[i] == 0) {
				this->Yflag[i] = -1;
			}
		}
	} else { // not hybrid rotation
		if (this->dynamic_peaks > 0) {
			for (i = 0; i < this->dynamic_peaks; ++i) {
				this->peaks_to_change.push_back(i);
				this->Yheight.push_back(random_double(0, 1));
				this->Yslope.push_back(random_double(0, 1));
				this->Yx.push_back(random_double(0, 1));
				this->Yy.push_back(random_double(0, 1));

				// normalization: flag must be either 1 or -1
				this->Hflag.push_back(rand() % 2);
				if (this->Hflag[i] == 0) {
					this->Hflag[i] = -1;
				}

				this->Rflag.push_back(rand() % 2);
				if (this->Rflag[i] == 0) {
					this->Rflag[i] = -1;
				}

				this->Xflag.push_back(rand() % 2);
				if (this->Xflag[i] == 0) {
					this->Xflag[i] = -1;
				}

				this->Yflag.push_back(rand() % 2);
				if (this->Yflag[i] == 0) {
					this->Yflag[i] = -1;
				}
			}
		}
	}

	this->pop_size = 0;
	this->solutions_on_peaks = 0;
	stats_file.open("./statistics/benchmark.txt");
	stats_file << "GENERATION" << "\t" << "MAX" << "\t" << "MIN" << "\t"
			<< "FREQUENCY" << "\t" << "SOLS_ON_PEAKS" << "\t" << "SOLS" << "\t"
			<< "COVERED_PEAKS" << "\t" << "PEAKS" << std::endl;

	if (this->debug) {
		this->reportLandscape(0);
	}

	return true; // initialization successful
}

// TODO: may not have the same number of global and local optima
void DF::changeHeights(double A, double scale) {
	for (unsigned int i = 0; i < this->dynamic_peaks; ++i) {
		// first we calculate the change to be applied
		double s = waveStep(this->H[this->peaks_to_change[i]], this->HBase,
				this->HLimit, this->Yheight[i], scale);

		// we verify that the peak won't go out of bounds
		if (H[peaks_to_change[i]] + s * this->Hflag[i] > this->HLimit
				|| H[peaks_to_change[i]] + s * this->Hflag[i] < this->HBase) {
			this->Hflag[i] *= -1;
		}

		//The flag will remain the same until the next time the peak goes out of bounds

		// we apply the change
		H[peaks_to_change[i]] += this->Hflag[i] * s;

		// we calculate the new value of Y using the logistics function
		this->Yheight[i] = logisticsFunction(A, this->Yheight[i]);
	}
}

// TODO: heights might change when modifying the slopes
void DF::changeSlopes(double A, double scale) {
	for (unsigned int i = 0; i < this->dynamic_peaks; ++i) {
		double s = waveStep(this->R[this->peaks_to_change[i]], this->RBase,
				this->RLimit, this->Yslope[i], scale);

		if (R[peaks_to_change[i]] + s * this->Rflag[i] > this->RLimit
				|| R[peaks_to_change[i]] + s * this->Rflag[i] < this->RBase) {
			this->Rflag[i] *= -1;
		}
		R[peaks_to_change[i]] += this->Rflag[i] * s;
		this->Yslope[i] = logisticsFunction(A, this->Yslope[i]);
	}
}

// TODO: optima might change if peaks go "out-of-bounds"
void DF::changeX(double A, double scale) {
	for (unsigned int i = 0; i < this->dynamic_peaks; ++i) {
		double s = waveStep(this->X[this->peaks_to_change[i]], this->XBase,
				this->XLimit, this->Yx[i], scale);
		if ((X[peaks_to_change[i]] + s * this->Xflag[i] > this->XLimit)
				|| (X[peaks_to_change[i]] + s * this->Xflag[i] < this->XBase)) {
			this->Xflag[i] = this->Xflag[i] * -1;
		}

		X[peaks_to_change[i]] += this->Xflag[i] * s;
		this->Yx[i] = logisticsFunction(A, this->Yx[i]);
	}
}

//TODO: optima might change if peaks go "out-of-bounds"
void DF::changeY(double A, double scale) {
	for (unsigned int i = 0; i < this->dynamic_peaks; ++i) {
		double s = waveStep(this->Y[this->peaks_to_change[i]], this->YBase,
				this->YLimit, this->Yy[i], scale);
		if (Y[peaks_to_change[i]] + s * this->Yflag[i] > this->YLimit
				|| Y[peaks_to_change[i]] + s * this->Yflag[i] < this->YBase) {
			this->Yflag[i] = this->Yflag[i] * -1;
		}
		Y[peaks_to_change[i]] += this->Yflag[i] * s;
		this->Yy[i] = logisticsFunction(A, this->Yy[i]);
	}
}

void DF::fast_rotation() {
	this->changeX(this->Ax_fast, this->scaleX_fast);
	this->changeY(this->Ay_fast, this->scaleY_fast);
}

void DF::slow_rotation() {
	this->changeX(this->Ax_slow, this->scaleX_slow);
	this->changeY(this->Ay_slow, this->scaleY_slow);
}

void DF::hybrid_rotation() {
	double s;
	unsigned int i;
	for (i = 0; i < fast_peaks_to_change.size(); ++i) {
		// X
		s = waveStep(this->X[this->fast_peaks_to_change[i]], this->XBase,
				this->XLimit, this->Yx[i], this->scaleX_fast);
		if (X[fast_peaks_to_change[i]] + s * this->Xflag[i] > this->XLimit
				|| X[fast_peaks_to_change[i]] + s * this->Xflag[i] < this->XBase) {
			this->Xflag[i] *= -1;
		}
		X[fast_peaks_to_change[i]] += this->Xflag[i] * s;
		this->Yx[i] = logisticsFunction(Ax_fast, this->Yx[i]);

		// Y
		s = waveStep(this->Y[this->fast_peaks_to_change[i]], this->YBase,
				this->YLimit, this->Yy[i], this->scaleY_fast);
		if (Y[fast_peaks_to_change[i]] + s * this->Yflag[i] > this->YLimit
				|| Y[fast_peaks_to_change[i]] + s * this->Yflag[i] < this->YBase) {
			this->Yflag[i] *= -1;
		}
		Y[fast_peaks_to_change[i]] += this->Yflag[i] * s;
		this->Yy[i] = logisticsFunction(Ay_fast, this->Yy[i]);
	}
	for (i = 0; i < slow_peaks_to_change.size(); ++i) {
		// X
		s = waveStep(this->X[this->slow_peaks_to_change[i]], this->XBase,
				this->XLimit, this->Yx[i], this->scaleX_slow);
		if (X[slow_peaks_to_change[i]] + s * this->Xflag[i] > this->XLimit
				|| X[slow_peaks_to_change[i]] + s * this->Yflag[i] < this->XBase) {
			this->Xflag[i] *= -1;
		}
		X[slow_peaks_to_change[i]] += this->Xflag[i] * s;
		this->Yx[i] = logisticsFunction(Ax_slow, this->Yx[i]);
		// Y
		s = waveStep(this->Y[this->slow_peaks_to_change[i]], this->YBase,
				this->YLimit, this->Yy[i], this->scaleY_slow);
		if (Y[slow_peaks_to_change[i]] + s * this->Yflag[i] > this->YLimit
				|| Y[slow_peaks_to_change[i]] + s * this->Yflag[i] < this->YBase) {
			this->Yflag[i] *= -1;
		}
		Y[slow_peaks_to_change[i]] += this->Yflag[i] * s;
		this->Yy[i] = logisticsFunction(Ay_slow, this->Yy[i]);
	}
}

void DF::dynamic_scaling() {
}

void DF::dynamic_local_optima() {
}

void DF::dynamic_elevate_optima() {
}

void DF::dynamic_global_optima() {
}

void DF::dynamic_downgrade_optima() {
}

double DF::evaluateFunction(const std::vector<double>& v) {
	// The AEON framework tries to maximizes the fitness function.
	double res = this->bottom; // this will force the landscape to have a flat "bottom"
	unsigned int peak = this->N + 1;
	for (unsigned int i = 0; i < this->N; ++i) {
		double form = ((v[0] - this->X[i]) * (v[0] - this->X[i])
				+ (v[1] - this->Y[i]) * (v[1] - this->Y[i]));
		if (this->cones) {
			form = sqrt(form);
		}
		double aux = H[i] - R[i] * form;
		if (aux > res) {
			res = aux;
			peak = i;
		}
	}
	// checking if this solution is on a peak (global optima)
	if ((peak < this->N - this->local_optima)
			&& (res >= this->peak_threshold * this->HLimit)) {
		this->solutions_on_peaks += 1;
		this->covered_peaks[peak] = true;
	}
	this->pop_size += 1;
	return res;
}

// The next functions are used only for debugging purposes

void DF::printFirstPeak() {
	std::ofstream file("severity", ios::app);
	file << this->H[this->peaks_to_change[0]] << " "
			<< this->R[this->peaks_to_change[0]] << " "
			<< this->X[this->peaks_to_change[0]] << " "
			<< this->Y[this->peaks_to_change[0]] << std::endl;
	file.close();
}

void DF::reportLandscape(int time) {
	std::stringstream s;
	double accuracy = 0.01;
	s << "./landscape/evaluation" << time;
	std::ofstream file;
	file.open(s.str().c_str());
	for (double i = this->XBase; i <= this->XLimit; i += accuracy) {
		for (double j = this->YBase; j <= this->YLimit; j += accuracy) {
			std::vector<double> v;
			v.push_back(i);
			v.push_back(j);

			file << i << "\t" << j << "\t" << this->evaluateFunction(v)
					<< std::endl;
		}
	}
	file.close();
}
