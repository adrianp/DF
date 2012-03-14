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
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/assign/std/vector.hpp> // used for vector += element1, element2;
#include <boost/shared_ptr.hpp>

// AEON includes
#include "solution/ParametricSolution.h"
#include "json-cpp/json/json.h"
#include "utils/utils/ConfFileParser.h"

// Dynamic experiment includes
#include "DF.h"
#include "dynamic_helpers.h"
#include "random_helper.h" // TODO: use AEON random

using namespace aeon;
using namespace boost::assign; // used for vector += element1, element2;

DF::DF() {
    this->XBase = -1;
    this->XLimit = 1;
    this->YBase = -1;
    this->YLimit = 1;
    this->generation = 0;
}

DF::~DF() {
    // TODO: do we need to do something here?
}

void DF::evaluate(std::vector<boost::shared_ptr
        <const aeon::ParametricSolution> >& solution) {
    this->result_ = this->evaluateFunction(solution[0]->getParameters());
}

bool DF::getResults(std::vector<boost::shared_ptr<const Result> >& results) {
    // TODO: hack for ResultData compliance
    std::vector<double> v;
    v += this->result_;
    results.push_back(boost::shared_ptr<const aeon::Result>(
           new aeon::Result(boost::shared_ptr<const aeon::ResultData>(
                   new aeon::ResultData(this->result_, v)))));
   return true;
}

void DF::setup(int evolutionaryTime, aeon::EvolutionaryTimeType::Value, int) {
    if(evolutionaryTime > this->generation && 
       this->dynamic_peaks > 0 && 
       evolutionaryTime % this->frequency == 0) {
        if(this->dynamic_type == "rotation") {
            this->dynamic_rotation();
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
        } else if (this->dynamic_type == "heights") {
            this->changeHeights();
        } else if (this->dynamic_type == "slopes") {
            this->changeSlopes();
        } else if (this->dynamic_type == "x") {
            this->changeX();
        } else if (this->dynamic_type == "y") {
            this->changeY();
        }
        if(this->debug) {
            // this will make things slow
            this->reportLandscape(evolutionaryTime); 
        }
        this->generation = evolutionaryTime;
    }
}

void DF::getParameters(std::vector<std::string>& parameters) {
    for (unsigned int i = 0; i < this->nParams_; ++i) {
       parameters.push_back("x" + i);
    }
}

bool DF::isDynamic() {
    if(this->dynamic_peaks > 0) {
        return true;
    } else {
        return false;
    }
}

bool DF::init(unsigned int seed, const std::string&, const std::string&) {
    std::srand(seed); // TODO: use AEON random
    unsigned int i;
    Json::Value rootNode;
    std::string value;
    if(!utils::ConfFileParser::parse("user_conf/dynamic_config.json", 
            rootNode, true)) {
        std::cout << "The dynamic experiment configuration file could not be "
                        "parsed." << std::endl;
        return false;
    } else {
        if(!rootNode["N"]) {
            std::cout << "Value for N could not be parsed." << std::endl;
            return false;
        } else {
            this->N = rootNode["N"].asInt();
        }
        if(!rootNode["local_optimas"]) {
            std::cout << "Value for local_optimas could not be parsed." 
                        << std::endl;
            return false;
        } else {
            this->local_optimas = rootNode["local_optimas"].asInt();
        }
        if(!rootNode["HBase"]) {
            std::cout << "Value for HBase could not be parsed." << std::endl;
            return false;
        } else {
            this->HBase = rootNode["HBase"].asDouble();
        }
        if(!rootNode["HLimit"]) {
            std::cout << "Value for HLimit could not be parsed." << std::endl;
            return false;
        } else {
            this->HLimit = rootNode["HLimit"].asDouble();
        }
        if(!rootNode["RBase"]) {
            std::cout << "Value for RBase could not be parsed." << std::endl;
            return false;
        } else {
            this->RBase = rootNode["RBase"].asDouble();
        }
        if(!rootNode["RLimit"]) {
            std::cout << "Value for RLimit could not be parsed." << std::endl;
            return false;
        } else {
            this->RLimit = rootNode["RLimit"].asDouble();
        }
        if(!rootNode["dynamic_peaks"]) {
            std::cout << "Value for dynamic_peaks could not be parsed." 
                        << std::endl;
            return false;
        } else {
            this->dynamic_peaks = rootNode["dynamic_peaks"].asInt();
        }
        if(!rootNode["cones"]) {
            std::cout << "Value for cones could not be parsed." << std::endl;
            return false;
        } else {
            this->cones = rootNode["cones"].asBool();
        }
        if(!rootNode["Ah"]) {
            std::cout << "Value for Ah could not be parsed." << std::endl;
            return false;
        } else {
            this->Ah = rootNode["Ah"].asDouble();
        }
        if(!rootNode["Ar"]) {
            std::cout << "Value for Ar could not be parsed." << std::endl;
            return false;
        } else {
            this->Ar = rootNode["Ar"].asDouble();
        }
        if(!rootNode["Ax"]) {
            std::cout << "Value for Ax could not be parsed." << std::endl;
            return false;
        } else {
            this->Ax = rootNode["Ax"].asDouble();
        }
        if(!rootNode["Ay"]) {
            std::cout << "Value for Ay could not be parsed." << std::endl;
            return false;
        } else {
            this->Ay = rootNode["Ay"].asDouble();
        }
        if(!rootNode["scaleH"]) {
            std::cout << "Value for scaleH could not be parsed." << std::endl;
            return false;
        } else {
            this->scaleH = rootNode["scaleH"].asDouble();
        }
        if(!rootNode["scaleR"]) {
            std::cout << "Value for scaleR could not be parsed." << std::endl;
            return false;
        } else {
            this->scaleR = rootNode["scaleR"].asDouble();
        }
        if(!rootNode["scaleX"]) {
            std::cout << "Value for scaleX could not be parsed." << std::endl;
            return false;
        } else {
            this->scaleX = rootNode["scaleX"].asDouble();
        }
        if(!rootNode["scaleY"]) {
            std::cout << "Value for scaleY could not be parsed." << std::endl;
            return false;
        } else {
            this->scaleY = rootNode["scaleY"].asDouble();
        }
        if(!rootNode["debug"]) {
            std::cout << "Value for debug could not be parsed." << std::endl;
            return false;
        } else {
            this->debug = rootNode["debug"].asBool();
        }
        if(!rootNode["dynamic_type"]) {
            std::cout << "Value for dynamic_type could not be parsed." 
                        << std::endl;
            return false;
        } else {
            std::vector<std::string> dynamic_types;
            dynamic_types += "rotation", "scaling", "local_optima", 
                             "elevate_optima", "global_optima", 
                             "downgrade_optima", "heights", "slopes", "x", "y";
            if (std::find(dynamic_types.begin(), dynamic_types.end(), 
                    rootNode["dynamic_type"].asString())!=dynamic_types.end()) {
                this->dynamic_type = rootNode["dynamic_type"].asString();
            } else {
                std::cout << "Specified dynamic change type could not be found." 
                                << std::endl;
                return false;
            }
        }
        if(!rootNode["frequency"]) {
            std::cout << "Value for frequency could not be parsed." << std::endl;
            return false;
        } else {
            this->frequency = rootNode["frequency"].asInt();
        }
    }

    // TODO: make the experiment work in more than 2 dimensions
    this->nParams_ = 2;
        
    if(this->local_optimas > this->N || this->dynamic_peaks > this->N) {
        std::cout << "N, global_optimas, dynamic_peaks must agree" << std::endl;
        return false;
    }
    
    if(this->dynamic_peaks > 0) {
        for(i=0; i < this->dynamic_peaks; ++i) {
            this->peaks_to_change.push_back(i);
            this->Yheight.push_back(0.1);
            this->Yslope.push_back(0.1);
            this->Yx.push_back(0.1);
            this->Yy.push_back(0.1);
            this->Hflag.push_back(1);
            this->Rflag.push_back(1);
            this->Xflag.push_back(1);
            this->Yflag.push_back(1);
        }
    }

  //double globalOpt = random_double(this->HBase, this->HLimit);
  double globalOpt = this->HLimit;

  for(i = 0; i < this->N; ++i) {
      if(i < this->N - this->local_optimas)
      {
          this->H.push_back(globalOpt);
      }
      else
      {
          this->H.push_back(random_double_exclusive(this->HBase, globalOpt));
      }
      this->R.push_back(random_double(this->RBase, this->RLimit));
      this->X.push_back(random_double(this->XBase,this->XLimit));
      this->Y.push_back(random_double(this->YBase,this->YLimit));
  }

  return true; // initialization successful 
}

void DF::changeHeights() {
    for(unsigned int i = 0; i < this->dynamic_peaks; ++i) {
        double s = waveStep(this->H[this->peaks_to_change[i]],
                            this->HBase,
                            this->HLimit,
                            this->Yheight[i],
                            this->scaleH);

        if(H[peaks_to_change[i]] + s*this->Hflag[i] > this->HLimit ||
           H[peaks_to_change[i]] + s*this->Hflag[i] < this->HBase) {
            this->Hflag[i] *= -1;
        }
        H[peaks_to_change[i]] += this->Hflag[i]*s;
        this->Yheight[i] = logisticsFunction(this->Ah, this->Yheight[i]);
    }
}

// TODO: investigate height changes
void DF::changeSlopes() {
    for(unsigned int i = 0; i < this->dynamic_peaks; ++i)
    {
        double s = waveStep(this->R[this->peaks_to_change[i]],
                            this->RBase,
                            this->RLimit,
                            this->Yslope[i],
                            this->scaleR);

        if(R[peaks_to_change[i]] + s*this->Rflag[i] > this->RLimit ||
           R[peaks_to_change[i]] + s*this->Rflag[i] < this->RBase)
        {
            this->Rflag[i] *= -1;
        }
        R[peaks_to_change[i]] += this->Rflag[i]*s;
        this->Yslope[i] = logisticsFunction(this->Ar, this->Yslope[i]);
    }
}

// TODO: heights might change if peaks go "out-of-bounds"
void DF::changeX() {
    for(unsigned int i = 0; i < this->dynamic_peaks; ++i) {
        double s = waveStep(this->X[this->peaks_to_change[i]],
                            this->XBase,
                            this->XLimit,
                            this->Yx[i],
                            this->scaleX);

        if(X[peaks_to_change[i]] + s*this->Xflag[i] > this->XLimit ||
           X[peaks_to_change[i]] + s*this->Xflag[i] < this->XBase) {
            this->Xflag[i] *= -1;
        }
        X[peaks_to_change[i]] += this->Xflag[i]*s;
        this->Yx[i] = logisticsFunction(this->Ax, this->Yx[i]);
    }
}

// TODO: heights might change if peaks go "out-of-bounds"
void DF::changeY() {
    for(unsigned int i = 0; i < this->dynamic_peaks; ++i) {
        double s = waveStep(this->Y[this->peaks_to_change[i]],
                            this->YBase,
                            this->YLimit,
                            this->Yy[i],
                            this->scaleY);

        if(Y[peaks_to_change[i]] + s*this->Yflag[i] > this->YLimit ||
           Y[peaks_to_change[i]] + s*this->Yflag[i] < this->YBase) {
            this->Yflag[i] *= -1;
        }
        Y[peaks_to_change[i]] += this->Yflag[i]*s;
        this->Yy[i] = logisticsFunction(this->Ay, this->Yy[i]);
    }
}

void DF::dynamic_rotation() {
    this->changeX();
    this->changeY();
}

void DF::dynamic_scaling() {
    // TODO: scale ALL the peaks
}

void DF::dynamic_local_optima() {
    // TODO: changes to local optimas, globals remain the same
}

void DF::dynamic_elevate_optima() {
    // TODO: local optimas become global optimas
}
void DF::dynamic_global_optima() {
    // TODO: changes to global optimas, globals remain the same
}
void DF::dynamic_downgrade_optima() {
    // TODO: global optimas become local optimas
}

// TODO: make the experiment work in more than 2 dimensions
double DF::evaluateFunction(const std::vector<double>& v) {
    //double res = -std::numeric_limits<double>::max();
    double res = 0; // this will force the landscape to have a flat "bottom"
    for (unsigned int i = 0; i < this->N; ++i) {
        double form = ((v[0]-this->X[i])*(v[0]-this->X[i])+(v[1]-this->Y[i])*(v[1]-this->Y[i]));
        if(this->cones) {
            form = sqrt(form);
        }
        double aux = H[i]-R[i]*form;
        if(aux > res) {
            res = aux;
        }
    }
    return res;
}

// The next functions are used only for debugging purposes

// this method can be used to print the location of ONE optimum
void DF::printOptimum() {
    double res = -std::numeric_limits<double>::max();
    unsigned int index = 0;
    for(unsigned int i=0; i<this->N; ++i) {
        if(this->H[i] > res) {
            res = H[i];
            index = i;
        }
    }
    std::cout << "Optimum:" << " "
              << this->H[index] << " "
              << this->X[index] << " "
              << this->Y[index] << std::endl;
}

// this method can be used to plot the severity of the change
double DF::getFirstChanged() {
    return this->H[peaks_to_change[0]];
}

void DF::reportLandscape(int time) {
    std::stringstream s;
    // TODO: easy config of output format
    //s << "landscape/evaluation" << time << ".csv";
    s << "landscape/evaluation" << time;
    std::ofstream file;
    file.open(s.str().c_str());
    for(double i=-1; i<=1; i+=0.01) {
        for(double j=-1; j<=1; j+=0.01) {
            std::vector<double> v;
            v.push_back(i);
            v.push_back(j);
            //file << i << "," << j << "," << this->evaluateFunction(v) << std::endl;
            file << i << " " << j << " " << this->evaluateFunction(v) << std::endl;
        }
    }
    file.close();
}