#include "experiment/IExperiment.h"
#include "DF.h"

extern "C" {

aeon::IExperiment *createExperiment() {
	DF* peaks = new DF();
	return peaks;
}

void freeExperiment(aeon::IExperiment *experiment) {
	delete experiment;
}

}
