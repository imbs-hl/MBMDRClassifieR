
#include <thread>
#include "Mbmdr.h"
#include "ModelClassification.h"

Mbmdr::Mbmdr() :
  data(0),
  order(0),
  n(0),
  alpha(0),
  max_models(0),
  num_models(0),
  feature_combination(0),
  predictions(0),
  logger(0) {
  this->mode = 1;
  this->models = std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers>();

  this->threads = 1;
}

Mbmdr::Mbmdr(Data* data,
             size_t order,
             double alpha,
             size_t max_models,
             size_t mode,
             size_t num_threads,
             Logger* logger) {
  this->mode = mode;
  this->data = data;
  this->order = order;
  this->alpha = alpha;
  this->max_models = max_models;
  this->num_models = 0;
  this->n = data->getNumObservations();
  this->models = std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers>();
  this->model_feature_names = std::unordered_set<std::string>();
  this->feature_combination = std::vector<size_t>(order);
  std::fill(this->feature_combination.begin(), this->feature_combination.end(), 0);
  this->predictions = std::vector<double>(0);

  if(num_threads == DEFAULT_NUM_THREADS) {
    // Detect cores
    this->threads = std::thread::hardware_concurrency();
  } else {
    this->threads = num_threads;
  }

  this->logger = logger;
}

Mbmdr::Mbmdr(Data* data, Rcpp::List saved_mbmdr,
             size_t num_threads,
             Logger* logger) {
  // Check if input is of correct type
  try {
    if(!saved_mbmdr.inherits("mbmdr")) {
      throw std::runtime_error("");
    }
  } catch (...) {
    throw std::runtime_error("Object must be a MB-MDR object.");
  }

  try {
    this->data = data;
    this->n = data->getNumObservations();
    std::fill(this->predictions.begin(), this->predictions.end(), 0);
    this->models = std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers>(CompareModelPointers(true));
    this->logger = logger;
  } catch (...) {
    throw std::runtime_error("Initialization of MB-MDR object failed.");
  }

  // Construct new models from saved MB-MDR object
  try {
    this->mode = saved_mbmdr["mode"];
    this->order = saved_mbmdr["order"];
    this->alpha = saved_mbmdr["alpha"];
    this->max_models = saved_mbmdr["max_models"];

    Rcpp::List models = saved_mbmdr["models"];
    this->predictions = std::vector<double>(n*models.size());

    for(int i = 0; i < models.size(); ++i) {
      Rcpp::List rcpp_model = models[i];

      double order = rcpp_model["order"];
      std::vector<size_t> features = rcpp_model["features"];
      std::vector<std::string> feature_names = rcpp_model["feature_names"];
      double alpha = rcpp_model["alpha"];
      std::vector<uint> in_cell = rcpp_model["in_cell"];
      std::vector<uint> out_cell = rcpp_model["out_cell"];
      std::vector<double> cell_predictions = rcpp_model["cell_predictions"];
      std::vector<double> cell_statistics = rcpp_model["cell_statistics"];
      std::vector<double> cell_pvalues = rcpp_model["cell_pvalues"];
      std::vector<int> cell_labels = rcpp_model["cell_labels"];
      double statistic = rcpp_model["statistic"];
      double pvalue = rcpp_model["pvalue"];

      Model* model;

      if(mode == 1) {
        // Create classification model
        logger->log(Info, "Creating classification model", 3);
        model = new ModelClassification(data,
                                        order,
                                        i,
                                        features,
                                        feature_names,
                                        alpha,
                                        logger);
      } else if(mode == 2) {
        // Create regression model
        logger->log(Info, "Creating regression model", 3);
      }

      model->loadModel(in_cell, out_cell, cell_predictions, cell_statistics, cell_pvalues, cell_labels, statistic, pvalue);

      this->models.push(model);

    }

    this->num_models = this->models.size();

  } catch(...) {
    throw std::runtime_error("Failed loading MB-MDR object.");
  }

  this->feature_combination = std::vector<size_t>(order);
  std::fill(this->feature_combination.begin(), this->feature_combination.end(), 0);

  if(num_threads == DEFAULT_NUM_THREADS) {
    // Detect cores
    this->threads = std::thread::hardware_concurrency();
  } else {
    this->threads = num_threads;
  }

}

Mbmdr::~Mbmdr() {
  while(!this->models.empty()) {
    Model* model = models.top();
    models.pop();
    delete model;
  }
}

void Mbmdr::possiblyAdd(Model* new_model) {

  logger->log(Info, "Comparing model statistics", 3);

  // Get model statistic of new model
  double new_model_statistic = new_model->getModelStatistic();

  // Insert in queue if not full yet
  logger->log(Info, "Statistic of new model: " + std::to_string(new_model_statistic), 3);
  if(models.size() < max_models) {
    logger->log(Info, "Queue not full yet. Adding new model.", 3);
    models.push(new_model);
  } else {
    // Get model statistic of model with lowest statistic in queue
    double top_model_statistic = models.top()->getModelStatistic();
    logger->log(Info, "Statistic of model at top position: " + std::to_string(top_model_statistic), 3);
    if(new_model_statistic > top_model_statistic) {
      // get rid of the root, i.e. feature model with lowest statistic
      logger->log(Info, "Adding new model", 3);
      Model* model = models.top();
      models.pop();
      delete model;

      // Add new model
      models.push(new_model);
    } else {
      // new model has lower statistic than all other models
      logger->log(Info, "Discarding new model", 3);
      delete new_model;
      return;
    }
  }

  // Adding feature names to set
  std::vector<size_t> features = new_model->getFeatures();
  for(size_t o = 0; o < order; ++o) {
    model_feature_names.insert(data->getFeatureName(features[o]));
  }
}

// [[Rcpp::plugins(cpp11)]]
bool Mbmdr::getNextFeatureCombination(size_t j) {

  size_t k = data->getNumFeatures();

  // Check if last feature combination is reached
  int last_feature_combination_sum = 0;
  int feature_combination_sum = 0;
  for(size_t i = 0; i < order; ++i) {
    last_feature_combination_sum += (k-i-1);
  }
  for(auto& n : feature_combination) {
    feature_combination_sum += n;
  }
  if (feature_combination_sum == last_feature_combination_sum) {
    return(false);
  }

  if (feature_combination[j] < (k-1) - ((order-1) - j)) {
    // Current index to increase is less than maximum at this index
    feature_combination[j] = feature_combination[j] + 1;
  } else {
    // Current index is at its maximum
    getNextFeatureCombination(j - 1);

    // Restart iteration in current index
    feature_combination[j] = feature_combination[j - 1] + 1;
  }

  return(true);

}

// [[Rcpp::plugins(cpp11)]]
void Mbmdr::fit() {

  // Create thread pool
  logger->log(Config, "Creating thread pool", 2);
  std::vector<std::thread> thread_pool;
  thread_pool.reserve(threads);

  // Add threads to thread pool
  for(size_t t=0; t<threads; ++t) {
    logger->log(Config, "Adding thread", 2);
    thread_pool.push_back(std::thread(&Mbmdr::fitModelInThread, this));
  }

  // Wait for completion
  logger->log(Info, "Waiting for threads to complete", 1);
  for(auto &thread : thread_pool) {
    thread.join();
  }


}

void Mbmdr::fitModelInThread() {

  std::vector<size_t> feature_combination;

  while(true) {

    // Feed thread with next feature combination
    std::unique_lock<std::mutex> lock(mutex);
    size_t model_index;
    if(getNextFeatureCombination(order-1)) {
      feature_combination = this->feature_combination;

      std::string message = "Calculating feature combination ";
      for(auto &col : feature_combination) {
        message += " " + std::to_string(col);
      }
      logger->log(Info, message, 3);

      // Increase model counter
      ++num_models;
      model_index = num_models;
    } else {
      logger->log(Info, "No further feature combinations available", 3);
      break;
    }
    lock.unlock();

    Model* model;

    if(mode == 1) {
      // Create classification model
      logger->log(Info, "Creating classification model", 3);
      model = new ModelClassification(data,
                                      order,
                                      model_index,
                                      feature_combination,
                                      alpha,
                                      logger);
    } else if(mode == 2) {
      // Create regression model
      logger->log(Info, "Creating regression model", 3);
    }

    // Fit the model
    logger->log(Info, "Model fit in progress...", 3);
    model->fit();

    // (Possibly) save model
    lock.lock();
    logger->log(Info, "Saving model", 3);
    possiblyAdd(model);
    lock.unlock();
  }

}

Rcpp::List Mbmdr::exportModels() {

  Rcpp::List export_object;

  while(!models.empty()) {
    // Get next model
    Model* model = models.top();

    // Initialize model export object
    Rcpp::List export_model_object;

    // Fill model export object
    export_model_object.push_back(model->getFeatures(), "features");
    export_model_object.push_back(model->getFeatureNames(), "feature_names");
    export_model_object.push_back(model->getObservationsInCell(), "in_cell");
    export_model_object.push_back(model->getObservationsOutCell(), "out_cell");
    export_model_object.push_back(model->getCellPredictions(), "cell_predictions");
    export_model_object.push_back(model->getCellStatistics(), "cell_statistics");
    export_model_object.push_back(model->getCellPValues(), "cell_pvalues");
    export_model_object.push_back(model->getCellLabels(), "cell_labels");
    export_model_object.push_back(model->getModelStatistic(), "statistic");
    export_model_object.push_back(model->getModelPValue(), "pvalue");
    export_model_object.push_back(model->getOrder(), "order");
    export_model_object.push_back(model->getAlpha(), "alpha");

    // Append model export object to list of models
    export_object.push_front(export_model_object);

    // Release model from queue
    models.pop();

    // Delete model
    delete model;
  }

  return export_object;

}


// [[Rcpp::plugins(cpp11)]]
std::vector<double> Mbmdr::predict() {

  // Create thread pool
  logger->log(Config, "Creating thread pool", 2);
  std::vector<std::thread> thread_pool;
  thread_pool.reserve(threads);

  // Add threads to thread pool
  for(size_t t=0; t<threads; ++t) {
    logger->log(Config, "Adding thread", 2);
    thread_pool.push_back(std::thread(&Mbmdr::predictInThread, this));
  }

  // Wait for completion
  logger->log(Info, "Waiting for threads to complete", 1);
  for(auto &thread : thread_pool) {
    thread.join();
  }

  return predictions;

}

void Mbmdr::predictInThread() {

  while(true) {

    // Feed thread with next model
    std::unique_lock<std::mutex> lock(mutex);
    Model* model;
    if(models.empty()) {
      break;
    } else {
      model = models.top();
      models.pop();
    }
    lock.unlock();

    // Predict sample outcomes in model
    logger->log(Info, "Predicting...", 3);
    std::vector<double> model_predictions = model->predict();

    // Add predictions
    lock.lock();
    logger->log(Info, "Saving predictions", 3);
    for(size_t i = 0; i < n; ++i) {
      predictions[i+n*model->getModelIndex()] += model_predictions[i];
    }
    lock.unlock();
    delete model;
  }

}

// Getters
uint Mbmdr::getMode() {
  return this->mode;
}
size_t Mbmdr::getOrder() {
  return this->order;
}
size_t Mbmdr::getN() {
  return this->n;
}
double Mbmdr::getAlpha() {
  return this->alpha;
}
size_t Mbmdr::getMaxModels() {
  return this->max_models;
}
size_t Mbmdr::getNumModels() {
  return this->num_models;
}
std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers> Mbmdr::getModels() {
  return this->models;
}
std::unordered_set<std::string> Mbmdr::getModelFeatureNames() {
  return this->model_feature_names;
}
