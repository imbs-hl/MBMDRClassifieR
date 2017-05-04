
#include "ModelClassification.h"
#include <math.h>

ModelClassification::ModelClassification() :
  Model(),
  cases(0), controls(0), mu(0) {
}

ModelClassification::ModelClassification(Data* data,
                                         size_t order,
                                         size_t model_index,
                                         std::vector<size_t> features,
                                         double alpha,
                                         Logger* logger) :
  Model(data, order, model_index, features, alpha, logger),
  mu(0) {
  size_t idxs = pow(3, order);
  for(size_t i = 0; i < idxs; ++i) {
    this->cases.push_back(0);
    this->controls.push_back(0);
    this->case_prob_in_cell.push_back(0);
    this->case_prob_out_cell.push_back(0);
  }
}

ModelClassification::ModelClassification(Data* data,
                                         size_t order,
                                         size_t model_index,
                                         std::vector<size_t> features,
                                         std::vector<std::string> feature_names,
                                         double alpha,
                                         Logger* logger) :
  Model(data, order, model_index, features, feature_names, alpha, logger),
  mu(0) {
  size_t idxs = pow(3, order);
  for(size_t i = 0; i < idxs; ++i) {
    this->cases.push_back(0);
    this->controls.push_back(0);
    this->case_prob_in_cell.push_back(0);
    this->case_prob_out_cell.push_back(0);
  }
}

ModelClassification::~ModelClassification() {
  // Empty on purpose
}

void ModelClassification::fit() {
  logger->log(Info, "Calculating cell counts", 3);
  getCounts();
  logger->log(Info, "Classifying cells", 3);
  classifyCells();
  logger->log(Info, "Calculating model statistic", 3);
  calculateModelTestStatistic();
}

void ModelClassification::getCounts() {

  // Get total number of observations
  size_t n_obs = data->getNumObservations();

  // Initialize base vector
  std::vector<uint> bases(order);
  for(size_t i = 0; i < order; ++i) {
    bases[i] = pow(3, i);
  }

  // Iterate through all samples in dataset
  logger->log(Info, "Iterating through samples...", 3);
  for(size_t i = 0; i < n_obs; ++i) {

    // Get genotype combination
    std::vector<uint> genotype(order);
    for(size_t j = 0; j < order; ++j) {
      genotype[j] = data->getFeature(i, features[j]);
    }
    // Convert genotype combination to index
    int idx = std::inner_product(genotype.begin(), genotype.end(), bases.begin(), 0);

    // TODO: add NA handling
    ++n;
    if(data->getOutcome(i) == 1) { // Factor level 1 is assumed to code controls
      // Sample is a control
      controls[idx] += 1;
    } else {
      // Sample is a case
      cases[idx] += 1;
    }
  }

  // Count cases
  uint sum_of_cases = 0;
  for (auto& cell_cases : cases)
    sum_of_cases += cell_cases;

  // Overall mean case rate
  mu = (double)sum_of_cases / (double)n;

}

// Classify cells
void ModelClassification::classifyCells() {

  // Number of cells
  int cells = pow(3, order);

  // Count cases
  uint sum_of_cases = 0;
  for (auto& cell_cases : cases)
    sum_of_cases += cell_cases;

  // Iterate through cells
  for(int i = 0; i < cells; ++i) {
    // Number of observations in cell
    in_cell[i] = cases[i] + controls[i];

    // Number of observations not in cell
    out_cell[i] = n - in_cell[i];

    // Proportion of cases in cell
    case_prob_in_cell[i] = (double)cases[i] / (double)in_cell[i];
    cell_predictions[i] = case_prob_in_cell[i];

    // Proportion of cases not in cell
    case_prob_out_cell[i] = (double)(sum_of_cases - cases[i]) / (double)out_cell[i];

    // Calculate cell statistic
    cell_statistics[i] = 1 / (mu*(1-mu)) *
      (in_cell[i] * pow(case_prob_in_cell[i] - mu, 2) +
      out_cell[i] * pow(case_prob_out_cell[i] - mu, 2));

    // Calculate cell P value
    cell_pvalues[i] = 1 - std::erf(std::sqrt(cell_statistics[i] / 2)); // For one degree of freedom the CDF of the chi square distribution is just the error function

    // Determine cell label
    if(cell_pvalues[i] < alpha) {
      if(case_prob_in_cell[i] > case_prob_out_cell[i]) {
        // Classify as H
        cell_labels[i] = 1;
      } else {
        // Classify as L
        cell_labels[i] = -1;
      }
    } else {
      // Classify as O
      cell_labels[i] = 0;
    }
  }
}

// Calculate model test statistic
void ModelClassification::calculateModelTestStatistic() {

  int cells = pow(3, order);

  int in_H = 0;
  int cases_in_H = 0;
  int controls_in_H = 0;
  int in_L = 0;
  int cases_in_L = 0;
  int controls_in_L = 0;
  int in_O = 0;
  int cases_in_O = 0;
  int controls_in_O = 0;

  for(int i = 0; i < cells; ++i) {
    if(cell_labels[i] == 1) {
      in_H += in_cell[i];
      cases_in_H += cases[i];
      controls_in_H += controls[i];
    } else if(cell_labels[i] == -1) {
      in_L += in_cell[i];
      cases_in_L += cases[i];
      controls_in_L += controls[i];
    } else {
      in_O += in_cell[i];
      cases_in_O += cases[i];
      controls_in_O += controls[i];
    }
  }

  double case_prob_in_H = (double)cases_in_H / (double)in_H;
  double case_prob_out_H = (double)(cases_in_L + cases_in_O) / (double)(in_L + in_O);
  double case_prob_in_L = (double)cases_in_L / (double)in_L;
  double case_prob_out_L = (double)(cases_in_H + cases_in_O) / (double)(in_H + in_O);

  double model_statistic_H = 0;
  double model_statistic_L = 0;

  model_statistic_H = 1 / (mu*(1-mu)) *
    (in_H * pow(case_prob_in_H - mu, 2) +
    (in_L + in_O) * pow(case_prob_out_H - mu, 2));
  model_statistic_L = 1 / (mu*(1-mu)) *
    (in_L * pow(case_prob_in_L - mu, 2) +
    (in_H + in_O) * pow(case_prob_out_L - mu, 2));

  if(std::isnan(model_statistic_H)) {
    if(std::isnan(model_statistic_L)) {
      statistic = 0;
    } else {
      statistic = model_statistic_L;
    }
  } else if(std::isnan(model_statistic_L)){
    statistic = model_statistic_H;
  } else {
    statistic = std::max(model_statistic_H, model_statistic_L);
  }
}
