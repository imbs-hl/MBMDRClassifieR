
#include <Rcpp.h>
#include "globals.h"
#include "Data.h"

Data::Data(Rcpp::IntegerMatrix& data, Rcpp::NumericVector& response) :
  outcome(response), data(data) {

  this->n_obs = data.nrow();
  this->n_feat = data.ncol();
  this->feature_names = std::unordered_map<std::string, size_t>();

  try {
    Rcpp::List dimnames = data.attr("dimnames");
    std::vector<std::string> names = dimnames[1];

    for(size_t j = 0; j < n_feat; ++j) {
      this->feature_names.emplace(std::make_pair(names[j], j));
    }
  } catch(std::exception& e) {
    Rcpp::Rcerr << "Failed to load colnames" << std::endl;
  }

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

int Data::getFeature(size_t sample, std::string feature) {
  return data(sample, feature_names[feature]);
}

std::string Data::getFeatureName(size_t feature) {
  Rcpp::List dimnames = data.attr("dimnames");
  std::vector<std::string> names = dimnames[1];
  return names[feature];
}

std::unordered_map<std::string, size_t> Data::getFeatureNames() {
  return feature_names;
}
