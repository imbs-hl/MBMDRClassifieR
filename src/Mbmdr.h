
#ifndef MBMDR_H_
#define MBMDR_H_

#include <Rcpp.h>
#include <queue>
#include <mutex>
#include <functional>
#include "globals.h"
#include "Model.h"
#include "Data.h"
#include "Logger.h"

class Mbmdr {

public:

  // Constructors
  Mbmdr();
  Mbmdr(Data* data,
        size_t order,
        size_t min_cell_size,
        double alpha,
        size_t max_models,
        size_t mode,
        size_t num_threads,
        Logger* logger);
  Mbmdr(Data* data,
        Rcpp::List saved_mbmdr,
        size_t num_threads,
        Logger* logger);

  // Destructor
  virtual ~Mbmdr();

  // Function to fit MB-MDR models
  void fit();

  // Function to predict MB-MDR models
  std::vector<double> predict();

  Rcpp::List exportModels();

  // Getters
  uint getMode() const;
  size_t getOrder() const;
  size_t getN() const;
  size_t getMinCellSize() const;
  double getAlpha() const;
  size_t getMaxModels() const;
  size_t getNumModels() const;
  std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers> getModels() const;
  std::unordered_set<std::string> getModelFeatureNames() const;

protected:

  // Pointer to original data
  Data* data;

  // Prediction mode
  uint mode;

  // Interaction order
  size_t order;

  // Total observations
  size_t n;

  // Minimum samples in cell to be relevant
  size_t min_cell_size;

  // Significance threshold for cell labels
  double alpha;

  // For multithreading
  std::mutex mutex;
  size_t threads;

  // Vector of feature models
  size_t max_models;
  size_t num_models;
  std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers> models;
  std::unordered_set<std::string> model_feature_names;

  // Current feature combination
  std::vector<size_t> feature_combination;

  // Predictions
  std::vector<double> predictions;

  // Add feature model to to priority queue, eventually replace model with lowest statistic
  void possiblyAdd(Model* new_model);

  // Iterate through all possible feature combinations
  bool getNextFeatureCombination(size_t j);

  // Logging
  Logger* logger;

private:
  void fitModelInThread();
  void predictInThread();

  DISALLOW_COPY_AND_ASSIGN(Mbmdr);

};

#endif /* MBMDR_H_ */
