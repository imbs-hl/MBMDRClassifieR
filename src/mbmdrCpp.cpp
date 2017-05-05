
#include <Rcpp.h>
#include "globals.h"
#include "Mbmdr.h"
#include "Data.h"
#include "Logger.h"

// [[Rcpp::export]]
Rcpp::List mbmdrCpp(size_t model_type,
                    Rcpp::IntegerMatrix& input_data, Rcpp::NumericVector& response,
                    size_t order, size_t min_cell_size, double alpha,
                    size_t max_results,
                    size_t num_threads,
                    size_t verbose,
                    std::string log_file,
                    Rcpp::List saved_mbmdr) {

  Rcpp::List result;
  Mbmdr* mbmdr = 0;
  Data* data = 0;

  // Logging
  Logger* logger = new Logger(log_file, verbose);

  try {

    // Initialize data
    logger->log(Config, "Initializing data", 1);
    data = new Data(input_data, response);

    // Initialize MB-MDR
    logger->log(Config, "Initializing MB-MDR", 1);
    if(saved_mbmdr.size() > 0) {
      // Load saved MB-MDR object
      logger->log(Info, "Loading saved MB-MDR", 1);
      mbmdr = new Mbmdr(data, saved_mbmdr,
                        num_threads,
                        logger);

      // Prediction
      logger->log(Info, "Predicting", 1);
      result.push_back(mbmdr->predict(), "predictions");
    } else {
      // Create new MB-MDR object
      logger->log(Info, "Creating new MB-MDR", 1);
      mbmdr = new Mbmdr(data,
                        order,
                        min_cell_size,
                        alpha,
                        max_results,
                        model_type,
                        num_threads,
                        logger);

      // Fit MB-MDR models
      logger->log(Info, "Fitting models", 2);
      mbmdr->fit();

      // Export MB-MDR object
      logger->log(Info, "Exporting", 3);
      Rcpp::List mbmdr_object;
      mbmdr_object.push_back(mbmdr->exportModels(), "models");
      mbmdr_object.push_back(mbmdr->getMinCellSize(), "min_cell_size");
      mbmdr_object.push_back(mbmdr->getAlpha(), "alpha");
      mbmdr_object.push_back(mbmdr->getMaxModels(), "max_models");
      mbmdr_object.push_back(mbmdr->getNumModels(), "num_models");
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
  delete logger;

  return result;

}
