
#include <Rcpp.h>
#include "globals.h"
#include "Data.h"

Data::Data(Rcpp::IntegerMatrix& data, Rcpp::NumericVector& response) :
outcome(response), data(data) {

	this->n_obs = data.nrow();
	this->n_feat = data.ncol();

}

Data::~Data() {
}

size_t Data::getNumObservations() {
	return this->n_obs;
}

size_t Data::getNumFeatures() {
	return this->n_feat;
}

double Data::getOutcome(size_t sample) {
	return outcome(sample);
}

int Data::getFeature(size_t sample, size_t feature) {
	return data(sample, feature);
}
