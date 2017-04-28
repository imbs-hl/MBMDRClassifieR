//
#include <Rcpp.h>
#include "globals.h"
#include "Mbmdr.h"
#include "Data.h"

// [[Rcpp::export]]
Rcpp::List mbmdrCpp(size_t pred_type,
                    Rcpp::IntegerMatrix& input_data, Rcpp::NumericVector& response,
                    size_t order, double alpha,
                    size_t max_results, size_t top_results,
                    size_t num_threads,
                    size_t verbose,
                    Rcpp::List saved_mbmdr) {

  Rcpp::List result;
  Mbmdr* mbmdr = 0;
  Data* data = 0;

  // Verbose streams
  std::vector<std::ostream*> v_levels;

  switch (verbose) {
  case 0:
    v_levels = {new std::stringstream, new std::stringstream, new std::stringstream};
    break;
  case 1:
    v_levels = {&std::cout, new std::stringstream, new std::stringstream};
    break;
  case 2:
    v_levels = {&std::cout, &std::cout, new std::stringstream};
    break;
  case 3:
    v_levels = {&std::cout, &std::cout, &std::cout};
    num_threads = 1;
    break;
  default:
    v_levels = {&std::cout, new std::stringstream, new std::stringstream};
  break;
  }

  try {

    // Initialize data
    *v_levels[1] << "Initializing data..." << std::endl;
    data = new Data(input_data, response);

    // Initialize MB-MDR
    *v_levels[1] << "Initializing MB-MDR..." << std::endl;
    if(saved_mbmdr.size() > 0) {
      // Load saved MB-MDR object
      *v_levels[2] << "Loading saved MB-MDR..." << std::endl;
      mbmdr = new Mbmdr(data, saved_mbmdr,
                        num_threads,
                        v_levels);

      // Prediction
      *v_levels[0] << "Predicting..." << std::endl;
      result.push_back(mbmdr->predict(), "predictions");
    } else {
      // Create new MB-MDR object
      *v_levels[2] << "Creating new MB-MDR..." << std::endl;
      mbmdr = new Mbmdr(data,
                        order,
                        alpha,
                        max_results,
                        pred_type,
                        num_threads,
                        v_levels);

      // Fit MB-MDR models
      *v_levels[0] << "Fitting models..." << std::endl;
      mbmdr->fit();

      // Export MB-MDR object
      *v_levels[1] << "Exporting top models..." << std::endl;
      Rcpp::List mbmdr_object;
      mbmdr_object.push_back(mbmdr->exportModels(), "models");
      mbmdr_object.push_back(mbmdr->getAlpha(), "alpha");
      mbmdr_object.push_back(mbmdr->getMaxModels(), "max_models");
      mbmdr_object.push_back(mbmdr->getMode(), "mode");
      mbmdr_object.push_back(mbmdr->getOrder(), "order");
      mbmdr_object.push_back(mbmdr->getModelFeatureNames(), "feature_names");
      mbmdr_object.attr("class") = "mbmdr";
      result.push_back(mbmdr_object, "mbmdr");
    }

  } catch(std::exception& e) {
    if(strcmp(e.what(), "User interrupt.") != 0) {
      Rcpp::Rcerr << e.what() << " MBMDRClassifieR will EXIT now." << std::endl;
    }
  }

  delete mbmdr;
  delete data;
  for (auto stream : v_levels) {
    if(stream != &std::cout) {
      delete stream;
    }
  }

  return result;

}
