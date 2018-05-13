#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "GxEDataset.h"

unsigned int MinimumCellCount = 10;

// Patch needed for some Linux compilers
namespace patch
{
template < typename T> std::string to_string(const T& n)
{
  std::ostringstream stm;
  stm << n;
  return stm.str();
}
}

// ***************************************************************************//
//                     Contstructor/Destructor                                //
// ***************************************************************************//

// Constructor
CGxEDataset::CGxEDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
                         const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
                         const bool *_filter) {
  
  m_initialized = false;
  m_good = false;
  m_errorString = "Dataset has not been initialized";
  
  m_numSubjects = _numSub;
  
  m_outcome = _outcome;
  m_missingOutcome = _missingOutcome;
  
  m_sourceCovariates = _cov;
  m_sourceMissingCovariates = _missingCov;
  
  m_numCovariates = _numCov + 1;
  m_numParameters = m_numCovariates + 2;
  m_maxParameters = m_numParameters;
  m_gxeNumber = _gxeNumber;
  
  m_covInt = NULL;
  m_missingCovariates = NULL;
  m_mean = NULL;
  m_std = NULL;
  
  m_gene = NULL;
  m_missingGene = NULL;
  m_geneFrequency = 0.;
  m_minMaf = 0.05;
  m_allelesSwapped = false;
  m_bDosages = false;
  m_bProbabilities = false;
  
  m_include = _filter;
  
  m_complete = NULL;
  m_numCompleteData = 0;
  m_used = NULL;
  m_numUsed = 0;
  
  m_xtransposex = NULL;
  
  m_beta = NULL;
  m_betaNew = NULL;
  
  m_score = NULL;
  
  m_information = NULL;
  m_informationUL = NULL;
  m_inverseInformation = NULL;
}
// Destructor
CGxEDataset::~CGxEDataset() {
#ifdef DEBUG_DELETE
  std::cout << "Entering GxEDataset Destructor" << std::endl;
#endif
  if (m_covInt)
    delete[] m_covInt;
  if (m_missingCovariates)
    delete[] m_missingCovariates;
  if (m_mean)
    delete[] m_mean;
  if (m_std)
    delete[] m_std;
  if (m_complete)
    delete[] m_complete;
  if (m_used)
    delete[] m_used;
  if (m_xtransposex)
    delete[] m_xtransposex;
  if (m_beta)
    delete[] m_beta;
  if (m_betaNew)
    delete[] m_betaNew;
  if (m_score)
    delete[] m_score;
  if (m_information)
    delete[] m_information;
  if (m_informationUL)
    delete[] m_informationUL;
  if (m_inverseInformation)
    delete[] m_inverseInformation;
}

// ***************************************************************************//
//                        Memory management                                   //
// ***************************************************************************//

// Clear allocated memory
void CGxEDataset::DeallocateMemory() {
#ifdef DATASET_DEBUG
  std::cout << "Deallocating GxEDataset memory" << std::endl;
#endif
  if (m_covInt)
    delete[] m_covInt;
  if (m_mean)
    delete[] m_mean;
  if (m_std)
    delete[] m_std;
  
  if (m_missingCovariates)
    delete[] m_missingCovariates;
  if (m_complete)
    delete[] m_complete;
  if (m_used)
    delete[] m_used;
  
  if (m_xtransposex)
    delete[] m_xtransposex;
  
  if (m_beta)
    delete[] m_beta;
  if (m_betaNew)
    delete[] m_betaNew;
  if (m_score)
    delete[] m_score;
  if (m_information)
    delete[] m_information;
  if (m_informationUL)
    delete[] m_informationUL;
  if (m_inverseInformation)
    delete[] m_inverseInformation;
  
  m_covInt = NULL;
  m_mean = NULL;
  m_std = NULL;
  
  m_missingCovariates = NULL;
  m_complete = NULL;
  m_used = NULL;
  
  m_xtransposex = NULL;
  
  m_beta = NULL;
  m_betaNew = NULL;
  m_score = NULL;
  m_information = NULL;
  m_informationUL = NULL;
  m_inverseInformation = NULL;
}
// Allocate memory needed for the score
bool CGxEDataset::AllocateScoreMemory() {
  m_score = new double[m_maxParameters];
  if (m_score == NULL)
    return false;
  return true;
}
// Allocate memory needed for the information
bool CGxEDataset::AllocateInformationMemory() {
  m_information = new double[m_maxParameters * m_maxParameters];
  m_informationUL = new double[m_maxParameters * m_maxParameters];
  if (m_information == NULL || m_informationUL == NULL)
    return false;
  
  return true;
}
// Allocate all the memory needed to do the analyses
bool CGxEDataset::AllocateMemory() {
#ifdef DATASET_DEBUG
  std::cout << "Allocating GxEDataset memory" << std::endl;
#endif
  DeallocateMemory();
  
  m_covInt = new double[m_numParameters * m_numSubjects];
  m_mean = new double[m_numParameters];
  m_std = new double[m_numParameters];
  if (m_covInt == NULL || m_mean == NULL || m_std == NULL)
    return false;
  
  m_missingCovariates = new bool[m_numSubjects];
  m_complete = new bool[m_numSubjects];
  m_used = new bool[m_numSubjects];
  if (m_missingCovariates == NULL || m_complete == NULL || m_used == NULL)
    return false;
  
  m_xtransposex = new double[m_numParameters * m_numParameters];
  if (m_xtransposex == NULL)
    return false;
  
  if (AllocateBetaMemory() == false)
    return false;
  
  if (AllocateScoreMemory() == false)
    return false;
  
  if (AllocateInformationMemory() == false)
    return false;
  
#ifdef DATASET_DEBUG
  std::cout << "Successfully allocated GxEDataset memory" << std::endl;
#endif
  
  return true;
}

// ***************************************************************************//
//                         Data management                                    //
// ***************************************************************************//

// Initialize the dataset
bool CGxEDataset::Initialize() {
#ifdef DATASET_DEBUG
  std::cout << "Entering GxEDataset intialize" << std::endl;
#endif
  m_initialized = true;
  if (m_numSubjects == 0) {
    m_good = false;
    m_errorString = "No subjects in dataset";
    return false;
  }
  
  m_good = AllocateMemory();
#ifdef DATASET_DEBUG
  std::cout << "GxEDataset allocated memory" << std::endl;
#endif
  
  if (m_good == false) {
    m_errorString = "Error allocating memory";
    return false;
  }
  
  AssignCovariateValues();
#ifdef DATASET_DEBUG
  std::cout << "GxEDataset assigned covariates" << std::endl;
#endif
  AssignComplete();
  Standardize();
  RunTests();
  return m_good;
}

// Copy the covariate values and create interaction covariate
void CGxEDataset::AssignCovariateValues() {
  unsigned int ui, uj;
  const double *covPtr;
  double *covIntPtr;
  const bool *miss1, *miss2;
  bool *miss3;
  unsigned int ncm1, ncsd;
  
  ncm1 = m_numCovariates - 1;
  ncsd = ncm1 * sizeof(double);
  memset(m_covInt, 0, m_numSubjects * m_numParameters * sizeof(double));
  memset(m_missingCovariates, false, m_numSubjects * sizeof(bool));
  
  covPtr = m_sourceCovariates;
  covIntPtr = m_covInt;
  miss1 = m_sourceMissingCovariates;
  miss3 = m_missingCovariates;
  for (ui = 0; ui < m_numSubjects; ++ui, covPtr += ncm1, miss1 += ncm1, ++miss3) {
    miss2 = miss1;
    for (uj = 0; uj < ncm1; ++uj, ++miss2) {
      if (*miss2 == true) {
        *miss3 = true;
        break;
      }
    }
    if (*miss3 == true) {
      covIntPtr += m_numParameters;
    }
    else {
      *covIntPtr = 1.;
      ++covIntPtr;
      memmove(covIntPtr, covPtr, ncsd);
      covIntPtr += ncm1;
      *covIntPtr = 1.;
      ++covIntPtr;
      *covIntPtr = covPtr[m_gxeNumber - 1];
      ++covIntPtr;
    }
  }
}
// Assign indicators to use subjects in analysis
void CGxEDataset::AssignComplete() {
  unsigned int ui;
  
  if (m_include != NULL)
    memmove(m_complete, m_include, m_numSubjects * sizeof(bool));
  else
    memset(m_complete, true, m_numSubjects * sizeof(bool));
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if (m_missingCovariates[ui] == true || m_missingOutcome[ui] == true)
      m_complete[ui] = false;
  }
}
// Standardize the covariates
void CGxEDataset::Standardize() {
  unsigned int ui, uj;
  unsigned int count;
  double *x;
  
  memset(m_mean, 0, m_numParameters * sizeof(double));
  memset(m_std, 0, m_numParameters * sizeof(double));
  // Disable standardization for testing purposes
  //	for (ui = 0; ui < m_numParameters; ++ui)
  //		m_std[ui] = 1;
  //	return;
  // Calculate the sum of X and X squared
  x = m_covInt;
  count = 0;
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if (m_complete[ui] == false) {
      x += m_numParameters;
      continue;
    }
    ++count;
    for (uj = 0; uj < m_numParameters; ++uj, ++x){
      m_mean[uj] += *x;
      m_std[uj] += *x * *x;
    }
  }
  
  // Calculate mean and std from the sums of X and X squared
  for (ui = 0; ui < m_numParameters; ++ui) {
    m_mean[ui] /= count;
    m_std[ui] = sqrt(m_std[ui] / count - m_mean[ui] * m_mean[ui]);
  }
  
  x = m_covInt;
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if (m_complete[ui] == false) {
      x += m_numParameters;
      continue;
    }
    for (uj = 0; uj < m_numParameters; ++uj, ++x) {
      if (m_std[uj] != 0)
        *x = (*x - m_mean[uj]) / m_std[uj];
    }
  }
}

// Assign gene to dataset
void CGxEDataset::AssignGene(const double *_gene, const bool *_missingGene, bool _dosages, bool _probabilities) {
  m_gene = _gene;
  m_missingGene = _missingGene;
  m_bDosages = _dosages;
  m_bProbabilities = _probabilities;
}
// Update the genetic values
void CGxEDataset::UpdateGene() {
  unsigned int ui;
  const double *g;
  const double *e;
  const bool *notMissing;
  double *gxe;
  
  memmove(m_used, m_complete, m_numSubjects *sizeof(bool));
  m_numUsed = 0;
  m_geneFrequency = 0.;
  m_allelesSwapped = false;
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if (m_missingGene[ui] == true)
      m_used[ui] = false;
    if (m_used[ui] == true) {
      ++m_numUsed;
      m_geneFrequency += m_gene[ui];
    }
  }
  if (m_numUsed == 0)
    m_geneFrequency = 0;
  else
    m_geneFrequency /= (m_numUsed + m_numUsed);
  if (m_geneFrequency > 0.5) {
    m_geneFrequency = 1. - m_geneFrequency;
    m_allelesSwapped = true;
  }
  
  e = m_covInt + m_gxeNumber;
  gxe = m_covInt + m_numParameters - 2;
  g = m_gene;
  notMissing = m_used;
  m_numUsed = 0;
  for (ui = 0; ui < m_numSubjects; ++ui, e += m_numParameters, gxe += m_numParameters, ++g, ++notMissing) {
    //		std::cout << ui << std::endl;
//    if (ui < 5)
//      std::cout << *g << '\t';
    if (*notMissing == false)
      continue;
    ++m_numUsed;
    if (m_allelesSwapped) {
      *gxe = (2 - *g);
      *(gxe + 1) = (2 - *g) * *e;
    }
    else {
      *gxe = *g;
      *(gxe + 1) = *g * *e;
    }
  }
//  std::cout << std::endl;
  //	std::cout << "Before X transpose X" << std::endl;
  //	XTransposeX(m_covInt, m_xtransposex, m_used, m_numParameters, m_numSubjects);
}

// ***************************************************************************//
//                           Data testing                                     //
// ***************************************************************************//

// Make sure at least one subject has complete data
void CGxEDataset::CompletenessTest() {
  unsigned int ui;
  
  m_numCompleteData = 0;
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if (m_complete[ui] == true)
      ++m_numCompleteData;
  }
  if (m_numCompleteData == 0) {
    m_errorString = "No subjects have complete phenotype and covariate data";
    m_good = false;
  }
}
// Make sure the covariates are linearly independent
int CGxEDataset::LinearDependenceTest() {
  m_good = true;
  XTransposeX(m_covInt, m_xtransposex, m_complete, m_numCovariates, m_numSubjects, m_numParameters);
  if (Invert(m_xtransposex, m_informationUL, m_inverseInformation, m_numCovariates, m_numUsed * 1.0e-12)) {
    m_errorString = "There is collinearity between at least two covariates in the model";
    m_good = false;
    return 1;
  }
  return 0;
}
// Run Tests for valid data
void CGxEDataset::RunTests() {
  CompletenessTest();
}

void CGxEDataset::Print(std::ostream &outfile, int n) {
  int i;
  unsigned int ui;
  double *d;
  bool *x;
  
  if (n == 0)
    n = m_numSubjects;
  d = m_covInt;
  x = (bool *)m_outcome;
  for (i = 0; i < n; ++i, ++x) {
    outfile << (*x ? 1 : 0) << '\t';
    for (ui = 0; ui < m_numCovariates + 2; ++ui, ++d) {
      outfile << *d << '\t';
    }
    outfile << std::endl;
  }
  
}

// ******************************************** mplink logistic regression ********************************************
// This does all the logistic regression, setting assign memory for all possible tests in mplink
// There should be a class for just a single logistic regression model.
// This model uses starting values from previous regressions as the starting point.
// This should be a subclass of a single logistic regresssion class. ???
// ********************************************************************************************************************

// ***************************************************************************//
//                     Contstructor/Destructor                                //
// ***************************************************************************//

// Constructor
CGxELogisticDataset::CGxELogisticDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
                                         const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
                                         const bool *_filter) : CGxEDataset(_numSub, _outcome, _missingOutcome, _numCov, _cov, _missingCov, _gxeNumber, _filter) {
  m_numCasesUsed = 0;
  m_numControlsUsed = 0;
  
  m_numParamD_E = 0;
  m_numParamD_G = 0;
  m_numParamD_GxE = 0;
  
  m_numParamSquaredD_E = 0;
  m_numParamSquaredD_G = 0;
  m_numParamSquaredD_GxE = 0;
  
  m_betaSize = 0;
  m_betaSizeLessDE = 0;
  m_invInfoSize = 0;
  
  m_numParamD_E = m_numCovariates;
  m_numParamD_G = m_numCovariates + 1;
  m_numParamD_GxE = m_numCovariates + 2;
  m_maxParameters = m_numParamD_GxE;
  
  m_numParamSquaredD_E = m_numParamD_E * m_numParamD_E;
  m_numParamSquaredD_G = m_numParamD_G * m_numParamD_G;
  m_numParamSquaredD_GxE = m_numParamD_GxE * m_numParamD_GxE;
  
  m_betaSize = m_numParamD_E + m_numParamD_G + m_numParamD_GxE;
  m_betaSizeLessDE = m_betaSize - m_numParamD_E;
  m_invInfoSize = m_numParamSquaredD_E + m_numParamSquaredD_G + m_numParamSquaredD_GxE;
  
  m_betaD_E = NULL;
  m_betaD_GE = NULL;
  m_betaD_GxE = NULL;
  
  m_scoreConstants = NULL;
  m_maxScore = 0.;
  
  m_inverseInformationD_E = NULL;
  m_inverseInformationD_GE = NULL;
  m_inverseInformationD_GxE = NULL;
}
// Destructor
CGxELogisticDataset::~CGxELogisticDataset() {
#ifdef DEBUG_DELETE
  std::cout << "Entering GxELogisticDataset Destructor" << std::endl;
#endif
  if (m_scoreConstants)
    delete[] m_scoreConstants;
  m_scoreConstants = NULL;
}

// ***************************************************************************//
//                        Memory management                                   //
// ***************************************************************************//

// Clear out allocated memory
void CGxELogisticDataset::DeallocateMemory() {
  CGxEDataset::DeallocateMemory();
  // The following were cleared out in CGxEDataset::Deallocate()
  m_betaD_E = NULL;
  m_betaD_GE = NULL;
  m_betaD_GxE = NULL;
  
  m_inverseInformationD_E = NULL;
  m_inverseInformationD_GE = NULL;
  m_inverseInformationD_GxE = NULL;
  // Clearing out what was added
  if (m_scoreConstants)
    delete[] m_scoreConstants;
  m_scoreConstants = NULL;
  m_maxScore = 0.;
}
// Allocate memory for betas
bool CGxELogisticDataset::AllocateBetaMemory() {
  
  m_beta = new double[m_betaSize];
  m_betaNew = new double[m_maxParameters];
  if (m_beta == NULL || m_betaNew == NULL)
    return false;
  
  m_betaD_E = m_beta;
  m_betaD_GE = m_betaD_E + m_numParamD_E;
  m_betaD_GxE = m_betaD_GE + m_numParamD_G;
  
  return true;
}
// Allocate memory for score
bool CGxELogisticDataset::AllocateScoreMemory() {
  m_score = new double[m_maxParameters];
  m_scoreConstants = new double[m_maxParameters];
  if (m_score == NULL || m_scoreConstants == NULL)
    return false;
  return true;
}
// Allocate memory for information and inverse information
bool CGxELogisticDataset::AllocateInformationMemory() {
  if (CGxEDataset::AllocateInformationMemory() == false)
    return false;
  
  m_inverseInformation = new double[m_invInfoSize];
  if (m_inverseInformation == NULL)
    return false;
  
  m_inverseInformationD_E = m_inverseInformation;
  m_inverseInformationD_GE = m_inverseInformationD_E + m_numParamSquaredD_E;
  m_inverseInformationD_GxE = m_inverseInformationD_GE + m_numParamSquaredD_G;
  
  return true;
}

// ***************************************************************************//
//                           Data testing                                     //
// ***************************************************************************//

// Test if there are there enough cases and controls to run the regression
void CGxELogisticDataset::CaseControlCountTest() {
  unsigned int ui;
  int numCases, numControls;
  const bool *outcome;
  
  numCases = 0;
  numControls = 0;
  outcome = (const bool *)m_outcome;
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if (m_complete[ui] == true) {
      if (outcome[ui] == true)
        ++numCases;
      else
        ++numControls;
    }
  }
  if (numCases < 10) {
    m_errorString = "There must be at least 10 cases in the dataset\nThere are only " + patch::to_string(numCases) + " cases in the current dataset";
    m_good = false;
  }
  else if (numControls < 10) {
    m_errorString = "There must be at least 10 controls in the dataset\nThere are only " + patch::to_string(numControls) + " controls in the current dataset";
    m_good = false;
  }
}
// Test if outcome can be fully explained by one covariate
void CGxELogisticDataset::FullyDeterminedTest() {
  double minCase, maxCase;
  double minControl, maxControl;
  unsigned int ui, uj;
  const double *cov;
  const bool *status;
  
  status = (const bool *)m_outcome;
  for (ui = 0; ui < m_numParameters; ++ui) {
    cov = m_covInt + ui;
    minCase = 1;
    maxCase = 0;
    minControl = 1;
    maxControl = 0;
    for (uj = 0; uj < m_numSubjects; ++uj, cov += m_numParameters) {
      if (m_complete[uj] == true) {
        if (status[uj] == true) {
          if (minCase > maxCase) {
            minCase = *cov;
            maxCase = *cov;
          }
          else if (*cov < minCase)
            minCase = *cov;
          else if (*cov > maxCase)
            maxCase = *cov;
        }
        else {
          if (minControl > maxControl) {
            minControl = *cov;
            maxControl = *cov;
          }
          else if (*cov < minControl)
            minControl = *cov;
          else if (*cov > maxControl)
            maxControl = *cov;
        }
      }
    }
    if (maxCase < minControl || maxControl < minCase) {
      m_errorString = "Disease status can be fully determined by covariate " + patch::to_string(ui);
      m_good = false;
      return;
    }
  }
}
// Peform tests if model is good - id est, convergence is possible
void CGxELogisticDataset::RunTests() {
  CompletenessTest();
  if (m_good == false)
    return;
  
  CaseControlCountTest();
  if (m_good == false)
    return;
  
  FullyDeterminedTest();
  if (m_good == false)
    return;
  
  LinearDependenceTest();
  if (m_good == false)
    return;
  
  InitializeBeta();
  if (m_good == false)
    return;
  
  memmove(m_used, m_complete, m_numSubjects * sizeof(bool));
  CalculateScoreConstants((bool *)m_outcome, m_numCovariates);
//  std::cout << "Before D|E" << std::endl;
  if (Logistic((bool *)m_outcome, m_numCovariates, m_betaD_E, m_inverseInformationD_E)) {
    m_errorString = "Unable to maximum model with no genes in model";
    m_good = false;
    std::cout << "D|E did not converged" << std::endl;
  }
//  else {
//    std::cout << m_betaD_E[0] << '\t' << m_betaD_E[1] << std::endl;
//    std::cout << "D|E converged" << std::endl;
//  }
}

// ***************************************************************************//
//                       Maximization routines                                //
// ***************************************************************************//

// Initialize betas to zero
int CGxELogisticDataset::InitializeBeta() {
  memset(m_beta, 0, (m_betaSize) * sizeof(double));
  
  return 0;
}
// Calculate the constant part of the score when disease status is the outcome
void CGxELogisticDataset::CalculateScoreConstants(const bool *outcome, unsigned int numCovar) {
  unsigned int ui, uj;
  const bool *used, *d;
  const double *x;
  double *sc;
  
  memset(m_scoreConstants, 0, m_numParameters*sizeof(double));
  used = m_used;
  d = outcome;
  x = m_covInt;
  for (ui = 0; ui < m_numSubjects; ++ui, ++used, ++d) {
    if (*used == false || *d == false) {
      x += m_numParameters;
      continue;
    }
    sc = m_scoreConstants;
    for (uj = 0; uj < numCovar; ++uj, ++x, ++sc)
      *sc += *x;
    x += m_numParameters - numCovar;
  }
}
// Calculate the constant part of the score when the gene is the outcome
void CGxELogisticDataset::CalculateScoreConstants(const double *outcome, unsigned int numCovar) {
  unsigned int ui, uj;
  const bool *used;
  const double *d;
  const double *x;
  double *sc;
  
  memset(m_scoreConstants, 0, m_numParameters*sizeof(double));
  used = m_used;
  d = outcome;
  x = m_covInt;
  for (ui = 0; ui < m_numSubjects; ++ui, ++used, ++d) {
    if (*used == false) {
      x += m_numParameters;
      continue;
    }
    sc = m_scoreConstants;
    for (uj = 0; uj < numCovar; ++uj, ++x, ++sc)
      *sc += *x * *d;
    x += m_numParameters - numCovar;
  }
  sc = m_scoreConstants;
  for (ui = 0; ui < numCovar; ++ui, ++sc)
    *sc /= 2.;
}
// Calculate the constant part of the score when the gene is the outcome selected by case control status
void CGxELogisticDataset::CalculateScoreConstants(const double *outcome, unsigned int numCovar, bool casesOnly) {
  unsigned int ui;
  const bool *complete;
  const bool *missingGene;
  bool *used;
  const bool *caseControl;
  unsigned int count;
  
  memset(m_used, true, m_numSubjects*sizeof(bool));
  
  complete = m_complete;
  missingGene = m_missingGene;
  caseControl = (bool *)m_outcome;
  used = m_used;
  count = m_numSubjects;
  for (ui = 0; ui < m_numSubjects; ++ui, ++complete, ++missingGene, ++caseControl, ++used) {
    if (*complete == false || *missingGene == true || *caseControl != casesOnly) {
      *used = false;
      --count;
    }
  }
  
  if (casesOnly == true)
    m_numCasesUsed = count;
  else
    m_numControlsUsed = count;
  
  CalculateScoreConstants(outcome, numCovar);
}
// Calculate the score and information - same for all models only difference is the score constants
void CGxELogisticDataset::CalculateScoreAndInformation(const bool *outcome, unsigned int numCovar, const double *beta) {
  unsigned int ui, uj, uk;
  double *score;
  double *info1, *info2;
  double expit, expit2;
  const bool *used;
  const double *beta1;
  const double *x1, *x2, *x3;
  const double *g;
  double ll;
  
  ll = 0.;
  memmove(m_score, m_scoreConstants, m_numParameters * sizeof(double));
  memset(m_information, 0, m_numParameters * m_numParameters * sizeof(double));
  x1 = m_covInt;
  used = m_used;
  g = m_gene;
  for (ui = 0; ui < m_numSubjects; ++ui, x1 += m_numParameters, ++g, ++used) {
    if (*used == false)
      continue;
    expit = 0;
    x2 = x1;
    beta1 = beta;
    for (uj = 0; uj < numCovar; ++uj, ++x2, ++beta1)
      expit += *beta1 * *x2;
    if (m_gene != NULL)
      ll += *g * expit;
    expit = exp(expit);
    if (m_gene != NULL)
      ll -= 2 * log(1 + expit);
    expit /= (1 + expit);
    expit2 = expit - expit * expit;
    x2 = x1;
    x3 = x1;
    score = m_score;
    info1 = m_information;
    for (uj = 0; uj < numCovar; ++uj, ++x2, ++score) {
      info1 += uj;
      *score -= *x2 * expit;
      x3 = x2;
      for (uk = uj; uk < numCovar; ++uk, ++info1, ++x3)
        *info1 += *x2 * *x3 * expit2;
    }
  }
  info1 = m_information;
  score = m_score;
  m_maxScore = 0;
  for (ui = 0; ui < numCovar; ++ui, ++score) {
    if (fabs(*score) > m_maxScore)
      m_maxScore = fabs(*score);
    info1 += ui;
    info2 = info1 + numCovar;
    ++info1;
    for (uj = ui + 1; uj < numCovar; ++uj, ++info1, info2 += numCovar)
      *info2 = *info1;
  }
}
// Update each beta parameter one by one until information can be inverted
int CGxELogisticDataset::UpdateOneByOne(const bool *outcome, unsigned int numCovar, double *beta, double *inverseInfo) {
  unsigned int ui, uj;
  double *b1;
  double *invInfo;
  double *sc;
  
  for (ui = 0; ui < 25; ++ui) {
    b1 = beta;
    invInfo = m_information;
    sc = m_score;
    for (uj = 0; uj < numCovar; ++uj, ++b1, ++sc, invInfo += numCovar + 1) {
      *b1 += *sc / *invInfo;
      CalculateScoreAndInformation(outcome, numCovar, beta);
      if (Invert(m_information, m_informationUL, inverseInfo, numCovar, m_numUsed * 1.0e-12) == 0)
        return 0;
    }
  }
  return 1;
}
// Fit a logistic model
int CGxELogisticDataset::Logistic(const bool *outcome, unsigned int numCovar, double *beta, double *inverseInfo) {
  unsigned int ui, uj;
  unsigned int fc;
  double *b1, *sc, *invInfo;
  bool improve = false;
  double oldMax;
  
  m_maxScore = 1.;
  CalculateScoreAndInformation(outcome, numCovar, beta);
  fc = 0;
  while (m_maxScore > 1e-6 && fc < 100) {
    improve = false;
    if (Invert(m_information, m_informationUL, inverseInfo, numCovar, m_numUsed * 1.0e-12) == 0) {
      oldMax = m_maxScore;
      memmove(m_betaNew, beta, numCovar * sizeof(double));
      b1 = m_betaNew;
      invInfo = inverseInfo;
      for (ui = 0; ui < numCovar; ++ui, ++b1) {
        sc = m_score;
        for (uj = 0; uj < numCovar; ++uj, ++sc, ++invInfo)
          *b1 += *sc * *invInfo;
      }
      CalculateScoreAndInformation(outcome, numCovar, m_betaNew);
      if (oldMax > m_maxScore) {
        memmove(beta, m_betaNew, numCovar * sizeof(double));
        improve = true;
      }
      ++fc;
    }
    if (improve == false) {
      if (UpdateOneByOne(outcome, numCovar, beta, inverseInfo))
        return 1;
      /*
      b1 = beta;
      invInfo = m_information;
      sc = m_score;
      for (ui = 0; ui < numCovar; ++ui, ++b1, ++sc, invInfo += numCovar + 1)
      *b1 += *sc / *invInfo;
      CalculateScoreAndInformation(outcome, numCovar, beta);
      */
    }
  }
  if (m_maxScore > 1e-6)
    return 1;
  if (Invert(m_information, m_informationUL, inverseInfo, numCovar, m_numUsed * 1.0e-12))
    return 1;
  
  return 0;
}

// Fit the models
int CGxELogisticDataset::FitModels() {
  if (m_initialized == false) {
    m_errorString = "Model not initialized";
    return 0x1fff;
  }
  
  if (m_geneFrequency < m_minMaf) {
    m_errorString = "Minor allele frequency under minimum value";
    return 0x2fff;
  }
  // No subjects with complete data
  if (m_numUsed == 0)
    return 0x03;
  
  // Zero all the betas
  std::fill(m_betaD_GE, m_betaD_GE + m_betaSizeLessDE, 0);
  
  // Initialize beta for D|GE to D|E betas with betaG = 0
  memmove(m_betaD_GE, m_betaD_E, m_numParamD_E*sizeof(double));
  CalculateScoreConstants((bool *)m_outcome, m_numParamD_G);
  if (Logistic((bool *)m_outcome, m_numParamD_G, m_betaD_GE, m_inverseInformationD_GE))
    return 0x03;
  
  // Initialize beta for D|GxE to D|G betas with betaGxE = 0
  memmove(m_betaD_GxE, m_betaD_GE, m_numParamD_G*sizeof(double));
  CalculateScoreConstants((bool *)m_outcome, m_numCovariates + 2);
  if (Logistic((bool *)m_outcome, m_numCovariates + 2, m_betaD_GxE, m_inverseInformationD_GxE))
    return 0x02;
  
  return 0;
}

// ******************************************** mplink polytomous regression ******************************************
// This does all the logistic and polytomous regression, setting aside memory for all possible tests in mplink
// This is derived from the gxe logistic regression. This adds the polytomous regression
// Comments about rearranging the classes still apply.
// ********************************************************************************************************************

// ***************************************************************************//
//                     Contstructor/Destructor                                //
// ***************************************************************************//

// Constructor
CGxEPolytomousDataset::CGxEPolytomousDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
                                             const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
                                             const bool *_filter) : CGxELogisticDataset(_numSub, _outcome, _missingOutcome, _numCov, _cov, _missingCov, _gxeNumber, _filter) {
  m_geCutoff = 1.96;
  m_numParamD_E = m_numCovariates;
  m_numParamD_G = m_numCovariates + 1;
  m_numParamD_GxE = m_numCovariates + 2;
  m_numParamHW = m_numCovariates;
  m_numParamPolytomous = m_numCovariates + m_numCovariates - 1;
  m_numParamRestrictedPolytomous = m_numCovariates + 1;
  if (m_numParamPolytomous > m_numParamD_GxE)
    m_maxParameters = m_numParamPolytomous;
  else
    m_maxParameters = m_numParamD_GxE;
  
  m_numParamSquaredD_E = m_numParamD_E * m_numParamD_E;
  m_numParamSquaredD_G = m_numParamD_G * m_numParamD_G;
  m_numParamSquaredD_GxE = m_numParamD_GxE * m_numParamD_GxE;
  m_numParamSquaredHW = m_numParamHW * m_numParamHW;
  m_numParamSquaredPolytomous = m_numParamPolytomous * m_numParamPolytomous;
  m_numParamSquaredRestrictedPolytomous = m_numParamRestrictedPolytomous * m_numParamRestrictedPolytomous;
  
  m_betaSize = m_numParamD_E + m_numParamD_G + m_numParamD_GxE + 3 * m_numParamHW + 3 * m_numParamPolytomous + 3 * m_numParamRestrictedPolytomous;
  m_betaSizeLessDE = m_betaSize - m_numParamD_E;
  m_invInfoSize = m_numParamSquaredD_E + m_numParamSquaredD_G + m_numParamSquaredD_GxE + 3 * m_numParamSquaredHW + 3 * m_numParamSquaredPolytomous + 3 * m_numParamSquaredRestrictedPolytomous;
  
  memset(m_geneCount, 0, sizeof(m_geneCount));
  m_polyScoreConstants = NULL;
  m_subjectXTX = NULL;
  
  m_betaG_E = NULL;
  m_betaCaseOnly = NULL;
  m_betaCntlOnly = NULL;
  m_betaPolytomousG_E = NULL;
  m_betaPolytomousCaseOnly = NULL;
  m_betaPolytomousControlOnly = NULL;
  m_betaRestrictedPolytomousG_E = NULL;
  m_betaRestrictedPolytomousCaseOnly = NULL;
  m_betaRestrictedPolytomousControlOnly = NULL;
  
  m_inverseInformationG_E = NULL;
  m_inverseInformationCaseOnly = NULL;
  m_inverseInformationCntlOnly = NULL;
  m_inverseInformationPolytomousG_E = NULL;
  m_inverseInformationPolytomousCaseOnly = NULL;
  m_inverseInformationPolytomousControlOnly = NULL;
  m_inverseInformationRestrictedPolytomousG_E = NULL;
  m_inverseInformationRestrictedPolytomousCaseOnly = NULL;
  m_inverseInformationRestrictedPolytomousControlOnly = NULL;
}
// Destructor
CGxEPolytomousDataset::~CGxEPolytomousDataset() {
#ifdef DEBUG_DELETE
  std::cout << "Entering GxEPolytomousDataset Destructor" << std::endl;
#endif
  if (m_polyScoreConstants)
    delete[] m_polyScoreConstants;
  if (m_subjectXTX)
    delete[] m_subjectXTX;
}

// ***************************************************************************//
//                        Memory management                                   //
// ***************************************************************************//

// Clear out allocated memory
void CGxEPolytomousDataset::DeallocateMemory() {
  CGxELogisticDataset::DeallocateMemory();
  // The following were cleared out in CGxELogisticDataset::Deallocate()
  m_betaG_E = NULL;
  m_betaCaseOnly = NULL;
  m_betaCntlOnly = NULL;
  m_betaPolytomousG_E = NULL;
  m_betaPolytomousCaseOnly = NULL;
  m_betaPolytomousControlOnly = NULL;
  m_betaRestrictedPolytomousG_E = NULL;
  m_betaRestrictedPolytomousCaseOnly = NULL;
  m_betaRestrictedPolytomousControlOnly = NULL;
  
  m_inverseInformationG_E = NULL;
  m_inverseInformationCaseOnly = NULL;
  m_inverseInformationCntlOnly = NULL;
  m_inverseInformationPolytomousG_E = NULL;
  m_inverseInformationPolytomousCaseOnly = NULL;
  m_inverseInformationPolytomousControlOnly = NULL;
  m_inverseInformationRestrictedPolytomousG_E = NULL;
  m_inverseInformationRestrictedPolytomousCaseOnly = NULL;
  m_inverseInformationRestrictedPolytomousControlOnly = NULL;
  // Clearing out what was added
  if (m_polyScoreConstants)
    delete[] m_polyScoreConstants;
  if (m_subjectXTX)
    delete[] m_subjectXTX;
  
  m_polyScoreConstants = NULL;
  m_subjectXTX = NULL;
}
// Allocate memory needed for betas
bool CGxEPolytomousDataset::AllocateBetaMemory() {
  if (CGxELogisticDataset::AllocateBetaMemory() == false)
    return false;
  m_betaG_E = m_betaD_GxE + m_numParamD_GxE;
  m_betaCaseOnly = m_betaG_E + m_numParamHW;
  m_betaCntlOnly = m_betaCaseOnly + m_numParamHW;
  m_betaPolytomousG_E = m_betaCntlOnly + m_numParamHW;
  m_betaPolytomousCaseOnly = m_betaPolytomousG_E + m_numParamPolytomous;
  m_betaPolytomousControlOnly = m_betaPolytomousCaseOnly + m_numParamPolytomous;
  m_betaRestrictedPolytomousG_E = m_betaPolytomousControlOnly + m_numParamPolytomous;
  m_betaRestrictedPolytomousCaseOnly = m_betaRestrictedPolytomousG_E + m_numParamRestrictedPolytomous;
  m_betaRestrictedPolytomousControlOnly = m_betaRestrictedPolytomousCaseOnly + m_numParamRestrictedPolytomous;
  
  return true;
}
// Allocate memory needed for score
bool CGxEPolytomousDataset::AllocateScoreMemory() {
  if (CGxELogisticDataset::AllocateScoreMemory() == false)
    return false;
  m_polyScoreConstants = new double[m_numParamPolytomous];
  if (m_polyScoreConstants == NULL)
    return false;
  
  return true;
}
// Allocate memory needed for information and inverse information
bool CGxEPolytomousDataset::AllocateInformationMemory() {
  if (CGxELogisticDataset::AllocateInformationMemory() == false)
    return false;
  
  m_inverseInformationG_E = m_inverseInformationD_GxE + m_numParamSquaredD_GxE;
  m_inverseInformationCaseOnly = m_inverseInformationG_E + m_numParamSquaredHW;
  m_inverseInformationCntlOnly = m_inverseInformationCaseOnly + m_numParamSquaredHW;
  m_inverseInformationPolytomousG_E = m_inverseInformationCntlOnly + m_numParamSquaredHW;
  m_inverseInformationPolytomousCaseOnly = m_inverseInformationPolytomousG_E + m_numParamSquaredPolytomous;
  m_inverseInformationPolytomousControlOnly = m_inverseInformationPolytomousCaseOnly + m_numParamSquaredPolytomous;
  m_inverseInformationRestrictedPolytomousG_E = m_inverseInformationPolytomousControlOnly + m_numParamSquaredPolytomous;
  m_inverseInformationRestrictedPolytomousCaseOnly = m_inverseInformationRestrictedPolytomousG_E + m_numParamSquaredRestrictedPolytomous;
  m_inverseInformationRestrictedPolytomousControlOnly = m_inverseInformationRestrictedPolytomousCaseOnly + m_numParamSquaredRestrictedPolytomous;
  
  return true;
}
// Allocate all memory needed for analysis
bool CGxEPolytomousDataset::AllocateMemory() {
  unsigned int xtxSize;
  
  if (CGxEDataset::AllocateMemory() == false)
    return false;
  
  xtxSize = m_numSubjects * (((m_numCovariates + 1) * m_numCovariates) / 2);
  m_subjectXTX = new double[xtxSize];
  if (m_subjectXTX == NULL)
    return false;
  
  return true;
}

// ***************************************************************************//
//                       Maximization routines                                //
// ***************************************************************************//

//		Polytomous logistic regression where beta1 = 2*beta2 only for covariate interacting with gene

// Calculate the constant part of the score when genes are measured
void CGxEPolytomousDataset::CalculatePolytomousMeasuredScoreConstants() {
  unsigned int ui, uj;
  unsigned int ncm1;
  const double *g, *x1, *x2;
  double *sc1, *sc2, *sci;
  const bool *used;
  
  memset(m_geneCount, 0, sizeof(m_geneCount));
  memset(m_polyScoreConstants, 0, m_numParamPolytomous * sizeof(double));
  g = m_gene;
  x1 = m_covInt;
  used = m_used;
  ncm1 = m_numCovariates - 1;
  sc2 = m_polyScoreConstants + ncm1;
  sci = m_polyScoreConstants + m_numParamPolytomous - 1;
  for (ui = 0; ui < m_numSubjects; ++ui, ++g, x1 += m_numParameters, ++used) {
    if (*used == false)
      continue;
    
    x2 = x1;
    
    if (*g == 1) {
      ++m_geneCount[1];
      sc1 = m_polyScoreConstants;
      for (uj = 0; uj < m_gxeNumber; ++uj, ++x2, ++sc1)
        *sc1 += *x2;
      *sci += *x2;
      ++uj;
      ++x2;
      for (; uj < m_numCovariates; ++uj, ++x2, ++sc1)
        *sc1 += *x2;
    }
    else if (*g == 2) {
      ++m_geneCount[2];
      sc1 = sc2;
      for (uj = 0; uj < m_gxeNumber; ++uj, ++x2, ++sc1)
        *sc1 += *x2;
      *sci += *x2 + *x2;
      ++uj;
      ++x2;
      for (; uj < m_numCovariates; ++uj, ++x2, ++sc1)
        *sc1 += *x2;
    }
    else {
      ++m_geneCount[0];
    }
  }
}
// Calculate the constant part of the score when using dosage values
void CGxEPolytomousDataset::CalculatePolytomousDosageScoreConstants() {
  unsigned int ui, uj;
  unsigned int ncm1;
  const double *x1, *x2;
  const double *p0, *p1, *p2;
  double *sc1, *sc2, *sci;
  const bool *used;
  
  memset(m_geneCount, 0, sizeof(m_geneCount));
  memset(m_polyScoreConstants, 0, m_numParamPolytomous * sizeof(double));
  p0 = m_gene + m_numSubjects;
  p1 = p0 + m_numSubjects;
  p2 = p1 + m_numSubjects;
  ncm1 = m_numCovariates - 1;
  sci = m_polyScoreConstants + m_numParamPolytomous - 1;
  x1 = m_covInt;
  used = m_used;
  for (ui = 0; ui < m_numSubjects; ++ui, ++p0, ++p1, ++p2, x1 += m_numParameters, ++used) {
    if (*used == false)
      continue;
    x2 = x1;
    m_geneCount[0] += *p0;
    m_geneCount[1] += *p1;
    m_geneCount[2] += *p2;
    sc1 = m_polyScoreConstants;
    sc2 = m_polyScoreConstants + ncm1;
    for (uj = 0; uj < m_gxeNumber; ++uj, ++x2, ++sc1, ++sc2) {
      *sc1 += *p1 * *x2;
      *sc2 += *p2 * *x2;
    }
    *sci += (*p1 + *p2 + *p2) * *x2;
    ++x2;
    ++uj;
    for (; uj < m_numCovariates; ++uj, ++x2, ++sc1, ++sc2) {
      *sc1 += *p1 * *x2;
      *sc2 += *p2 * *x2;
    }
  }
}
// Calculate the score and information
void CGxEPolytomousDataset::CalculatePolytomousScoreAndInfo(const double *_beta) {
  const double *beta1, *beta2;
  double betai, betai2;
  double betax1, betax2;
  double spart1, spart2;
  double *sc1, *sc2, *sci;
  double *info1, *info2, *info3, *info4, *info5, *info6;
  double *info1a, *info2a, *info3a, *info4a, *info5a, *info6a;
  double *info1b, *info2b, *info3b, *info4b;
  double ipart1, ipart2, ipart3, ipart4, ipart5;
  double denom;
  const double *x1, *x2, *xtx1, *xtx2;
  const bool *used;
  unsigned int xtxSize;
  unsigned int nc2, nc2m1;
  unsigned int ui, uj, uk;
  
  xtxSize = (m_numCovariates * (m_numCovariates + 1)) / 2;
  nc2 = m_numCovariates + m_numCovariates;
  nc2m1 = nc2 - 1;
  sci = m_score + nc2m1 - 1;
  //	betai = 0.1;
  betai = *(_beta + nc2m1 - 1);
  betai2 = betai + betai;
  x1 = m_covInt;
  xtx1 = m_subjectXTX;
  used = m_used;
  memset(m_score, 0, m_numParamPolytomous * sizeof(double));
  memmove(m_score, m_polyScoreConstants, m_numParamPolytomous * sizeof(double)); // ??? correct but commented out for testing
  memset(m_information, 0, m_numParamSquaredPolytomous*sizeof(double));
  for (ui = 0; ui < m_numSubjects; ++ui, x1 += m_numParameters, xtx1 += xtxSize, ++used) {
    if (*used == false)
      continue;
    
    x2 = x1;
    beta1 = _beta;
    beta2 = beta1 + m_numCovariates - 1;
    betax1 = 0;
    betax2 = 0;
    for (uj = 0; uj < m_gxeNumber; ++uj, ++x2, ++beta1, ++beta2) {
      betax1 += *beta1 * *x2;
      betax2 += *beta2 * *x2;
    }
    betax1 += betai * *x2;
    betax2 += betai2 * *x2;
    ++uj;
    ++x2;
    for (; uj < m_numCovariates; ++uj, ++x2, ++beta1, ++beta2) {
      betax1 += *beta1 * *x2;
      betax2 += *beta2 * *x2;
    }
    spart1 = exp(betax1);
    spart2 = exp(betax2);
    denom = 1 + spart1 + spart2;
    spart1 /= denom;
    spart2 /= denom;
    ipart1 = spart1 - spart1 * spart1;
    ipart2 = -spart1 * spart2;
    ipart3 = spart2 - spart2 * spart2;
    ipart4 = ipart1 + ipart2 + ipart2;
    ipart5 = ipart3 + ipart3 + ipart2;
    
    x2 = x1;
    sc1 = m_score;
    sc2 = m_score + m_numCovariates - 1;
    info1 = m_information;
    info2 = m_information + m_numCovariates - 1;
    info3 = m_information + nc2 * (m_numCovariates - 1);
    xtx2 = xtx1;
    
    for (uj = 0; uj < m_gxeNumber; ++uj, ++x2, ++sc1, ++sc2, info1 += nc2, info2 += nc2, info3 += nc2) {
      *sc1 -= spart1 * *x2;
      *sc2 -= spart2 * *x2;
      info1a = info1;
      info2a = info2;
      info3a = info3;
      for (uk = uj; uk < m_numCovariates - 1; ++uk, ++info1a, ++info2a, ++info3a, ++xtx2) {
        *info1a += ipart1 * *xtx2;
        *info2a += ipart2 * *xtx2;
        *info3a += ipart3 * *xtx2;
      }
    }
    
    *sci -= (spart1 + spart2 + spart2) * *x2;
    ++x2;
    
    for (; uj < m_numCovariates - 1; ++uj, ++x2, ++sc1, ++sc2, info1 += nc2, info2 += nc2, info3 += nc2) {
      *sc1 -= spart1 * *x2;
      *sc2 -= spart2 * *x2;
      info1a = info1;
      info2a = info2;
      info3a = info3;
      for (uk = uj; uk < m_numCovariates - 1; ++uk, ++info1a, ++info2a, ++info3a, ++xtx2) {
        *info1a += ipart1 * *xtx2;
        *info2a += ipart2 * *xtx2;
        *info3a += ipart3 * *xtx2;
      }
    }
    
    info1 = m_information + nc2m1 - 1;
    info2 = info1 + nc2m1 * (m_numCovariates - 1);
    for (uj = 0; uj < m_numCovariates - 1; ++uj, ++xtx2, info1 += nc2m1, info2 += nc2m1) {
      *info1 += ipart4 * *xtx2;
      *info2 += ipart5 * *xtx2;
    }
    *info2 += (ipart4 + ipart5 + ipart5) * *xtx2;
    
  }
  
  info1 = m_information;
  info2 = m_information + m_numCovariates - 1;
  info3 = m_information + nc2m1 * (m_numCovariates - 1);
  info4 = m_information + nc2 * (m_numCovariates - 1);
  info5 = m_information + nc2m1 - 1;
  info5a = info5 + nc2m1 * (m_numCovariates - 1);
  info6 = m_information + nc2m1 * (nc2m1 - 1);
  info6a = info6 + (m_numCovariates - 1);
  for (uj = 0; uj < m_numCovariates - 1; ++uj, info1 += nc2, info2 += nc2, info3 += nc2, info4 += nc2, info5 += nc2m1, info5a += nc2m1, ++info6, ++info6a) {
    info1a = info1;
    info2a = info2;
    info3a = info3;
    info4a = info4;
    info1b = info1;
    info2b = info2;
    info3b = info3;
    info4b = info4;
    for (uk = uj; uk < m_numCovariates - 1; ++uk, ++info1a, ++info2a, ++info3a, ++info4a, info1b += nc2m1, info2b += nc2m1, info3b += nc2m1, info4b += nc2m1) {
      *info1b = *info1a;
      *info2b = *info2a;
      *info3a = *info2a;
      *info3b = *info3a;
      *info4b = *info4a;
    }
    *info6 = *info5;
    *info6a = *info5a;
  }
  sc1 = m_score;
  m_maxScore = 0;
  for (ui = 0; ui < nc2m1; ++ui, ++sc1) {
    if (fabs(*sc1) > m_maxScore)
      m_maxScore = fabs(*sc1);
  }
}
// Fit the polytomous logistic model
int CGxEPolytomousDataset::PolytomousLogistic(double *_beta, double *_inverseInfo) {
  unsigned int ui, uj;
  unsigned int fc;
  unsigned int nc2m1;
  double *b1, *sc, *invInfo;
  bool improve = false;
  double oldMax;
  
  nc2m1 = m_numCovariates + m_numCovariates - 1;
  m_maxScore = 1.;
  CalculatePolytomousScoreAndInfo(_beta);
  fc = 0;
  while (m_maxScore > 1e-6 && fc < 100) {
    improve = false;
    if (Invert(m_information, m_informationUL, _inverseInfo, nc2m1, m_numUsed * 1.0e-12) == 0) {
      oldMax = m_maxScore;
      memmove(m_betaNew, _beta, nc2m1 * sizeof(double));
      b1 = m_betaNew;
      invInfo = _inverseInfo;
      for (ui = 0; ui < nc2m1; ++ui, ++b1) {
        sc = m_score;
        for (uj = 0; uj < nc2m1; ++uj, ++sc, ++invInfo)
          *b1 += *sc * *invInfo;
      }
      CalculatePolytomousScoreAndInfo(m_betaNew);
      if (oldMax > m_maxScore) {
        memmove(_beta, m_betaNew, nc2m1 * sizeof(double));
        improve = true;
      }
      ++fc;
    }
    if (improve == false) {
      // Need to come up with an idea on what to do if information is not invertable ???
      return 1;
    }
  }
  if (m_maxScore > 1e-6)
    return 1;
  if (Invert(m_information, m_informationUL, _inverseInfo, nc2m1, m_numUsed * 1.0e-12))
    return 1;
  
  return 0;
}

//		Polytomous logistic regression where beta1 = 2*beta2 for all covariates

// Calculate the constant part of the score when genes are measured
void CGxEPolytomousDataset::CalculateRestrictedPolytomousMeasuredScoreConstants() {
  unsigned int ui, uj;
  const double *g, *x1, *x2;
  double *sc1, *sc2, *sci;
  const bool *used;
  
  memset(m_geneCount, 0, sizeof(m_geneCount));
  memset(m_polyScoreConstants, 0, m_numParamRestrictedPolytomous * sizeof(double));
  g = m_gene;
  x1 = m_covInt;
  used = m_used;
  sc2 = m_polyScoreConstants + m_numCovariates - 1;
  sci = sc2 + 1;
  for (ui = 0; ui < m_numSubjects; ++ui, ++g, x1 += m_numParameters, ++used) {
    if (*used == false)
      continue;
    
    x2 = x1;
    
    if (*g == 1) {
      ++m_geneCount[1];
      sc1 = m_polyScoreConstants;
      for (uj = 0; uj < m_gxeNumber; ++uj, ++x2, ++sc1)
        *sc1 += *x2;
      *sci += *x2;
      ++uj;
      ++x2;
      for (; uj < m_numCovariates; ++uj, ++x2, ++sc1)
        *sc1 += *x2;
    }
    else if (*g == 2) {
      ++m_geneCount[2];
      *sc2 += *x2;
      ++x2;
      sc1 = m_polyScoreConstants + 1;
      for (uj = 1; uj < m_gxeNumber; ++uj, ++x2, ++sc1)
        *sc1 += *x2 + *x2;
      *sci += *x2 + *x2;
      ++uj;
      ++x2;
      for (; uj < m_numCovariates; ++uj, ++x2, ++sc1)
        *sc1 += *x2 + *x2;
    }
    else {
      ++m_geneCount[0];
    }
  }
}
// Calculate the constant part of the score when using dosage values
void CGxEPolytomousDataset::CalculateRestrictedPolytomousDosageScoreConstants() {
  unsigned int ui, uj;
  unsigned int ncm1;
  const double *x1, *x2;
  const double *p0, *p1, *p2;
  const double *d;
  double *sc1, *sc2, *sci;
  const bool *used;
  
  memset(m_geneCount, 0, sizeof(m_geneCount));
  memset(m_polyScoreConstants, 0, m_numParamPolytomous * sizeof(double));
  d = m_gene;
  p0 = m_gene + m_numSubjects;
  p1 = p0 + m_numSubjects;
  p2 = p1 + m_numSubjects;
  ncm1 = m_numCovariates - 1;
  sc2 = m_polyScoreConstants + ncm1;
  sci = sc2 + 1;
  x1 = m_covInt;
  used = m_used;
  for (ui = 0; ui < m_numSubjects; ++ui, ++d, ++p0, ++p1, ++p2, x1 += m_numParameters, ++used) {
    if (*used == false)
      continue;
    x2 = x1;
    m_geneCount[0] += *p0;
    m_geneCount[1] += *p1;
    m_geneCount[2] += *p2;
    *m_polyScoreConstants += *p1 * *x2;
    *sc2 += *p2 * *x2;
    ++x2;
    sc1 = m_polyScoreConstants + 1;
    for (uj = 1; uj < m_gxeNumber; ++uj, ++x2, ++sc1)
      *sc1 += *d * *x2;
    *sci += *d * *x2;
    ++x2;
    ++uj;
    for (; uj < m_numCovariates; ++uj, ++x2, ++sc1)
      *sc1 += *d * *x2;
  }
}
// Calculate the score and information
void CGxEPolytomousDataset::CalculateRestrictedPolytomousScoreAndInfo(const double *_beta) {
  const double *beta1;
  double beta2, bx;
  double betai;
  double betax1, betax2;
  double spart1, spart2, spart3;
  double *sc1, *sc2, *sci;
  double *info1, *info2, *info3;
  double *info1a, *info2a;
  double ipart1, ipart2, ipart3, ipart4, ipart5, ipart6;
  double denom;
  const double *x1, *x2, *xtx1, *xtx2;
  const bool *used;
  unsigned int xtxSize;
  unsigned int ui, uj, uk;
  
  xtxSize = (m_numCovariates * (m_numCovariates + 1)) / 2;
  sci = m_score + m_numCovariates;
  //	betai = 0.1;
  beta2 = *(_beta + (m_numCovariates - 1));
  betai = *(_beta + m_numCovariates);
//  betai2 = betai + betai;
  x1 = m_covInt;
  xtx1 = m_subjectXTX;
  used = m_used;
  memset(m_score, 0, m_numParamRestrictedPolytomous * sizeof(double));
  memmove(m_score, m_polyScoreConstants, (m_numCovariates + 1) * sizeof(double)); // ??? correct but commented out for testing
  memset(m_information, 0, m_numParamSquaredRestrictedPolytomous * sizeof(double));
  for (ui = 0; ui < m_numSubjects; ++ui, x1 += m_numParameters, xtx1 += xtxSize, ++used) {
    if (*used == false)
      continue;
    
    x2 = x1;
    beta1 = _beta;
    betax1 = *beta1 * *x2;
    betax2 = beta2 * *x2;
    ++x2;
    ++beta1;
    for (uj = 1; uj < m_gxeNumber; ++uj, ++x2, ++beta1) {
      bx = *beta1 * *x2;
      betax1 += bx;
      betax2 += bx + bx;
    }
    bx = betai * *x2;
    betax1 += bx;
    betax2 += bx + bx;
    ++uj;
    ++x2;
    for (; uj < m_numCovariates; ++uj, ++x2, ++beta1) {
      bx = *beta1 * *x2;
      betax1 += bx;
      betax2 += bx + bx;
    }
    spart1 = exp(betax1);
    spart2 = exp(betax2);
    denom = 1 + spart1 + spart2;
    spart1 /= denom;
    spart2 /= denom;
    spart3 = spart1 + spart2 + spart2;
    ipart1 = spart1 - spart1 * spart1;
    ipart2 = -spart1 * spart2;
    ipart3 = spart2 - spart2 * spart2;
    ipart4 = ipart1 + ipart2 + ipart2;
    ipart5 = ipart3 + ipart3 + ipart2;
    ipart6 = ipart4 + ipart5 + ipart5;
    
    x2 = x1;
    sc1 = m_score;
    sc2 = m_score + m_numCovariates - 1;
    info1 = m_information;
    info2 = m_information + m_numCovariates - 1;
    info3 = m_information + (m_numCovariates + 2) * (m_numCovariates - 1);
    xtx2 = xtx1;
    
    *sc2 -= spart2 * *x2;
    *sc1 -= spart1 * *x2;
    *info3 += ipart3 * *xtx2;
    *info1 += ipart1 * *xtx2;
    *info2 += ipart2 * *xtx2;
    
    info1a = info1 + 1;
    info2a = info2 + m_numCovariates + 1;
    ++xtx2;
    for (uj = 1; uj < m_numCovariates - 1; ++uj, ++info1a, info2a += m_numCovariates + 1, ++xtx2) {
      *info1a += ipart4 * *xtx2;
      *info2a += ipart5 * *xtx2;
    }
    
    ++x2;
    ++sc1;
    info1 += m_numCovariates + 2;
    info2 += m_numCovariates + 1;
    
    for (uj = 1; uj < m_gxeNumber; ++uj, ++x2, ++sc1, info1 += m_numCovariates + 2) {
      *sc1 -= spart3 * *x2;
      info1a = info1;
      for (uk = uj; uk < m_numCovariates - 1; ++uk, ++info1a, ++xtx2)
        *info1a += ipart6 * *xtx2;
    }
    
    *sci -= spart3 * *x2;
    ++x2;
    
    for (; uj < m_numCovariates - 1; ++uj, ++x2, ++sc1, info1 += m_numCovariates + 2, info2 += m_numCovariates + 1) {
      *sc1 -= spart3 * *x2;
      info1a = info1;
      for (uk = uj; uk < m_numCovariates - 1; ++uk, ++info1a, ++xtx2)
        *info1a += ipart6 * *xtx2;
    }
    
    info1 = m_information + m_numCovariates;
    *info1 += ipart4 * *xtx2;
    *(m_information + (m_numCovariates + 1) * m_numCovariates - 1) += ipart5 * *xtx2;
    info1 += m_numCovariates + 1;
    ++xtx2;
    for (uj = 1; uj < m_numCovariates - 1; ++uj, ++xtx2, info1 += m_numCovariates + 1)
      *info1 += ipart6 * *xtx2;
    info1 += m_numCovariates + 1;
    *info1 += ipart6 * *xtx2;
  }
  
  info1 = m_information + 1;
  info2 = m_information + m_numCovariates + 1;
  for (uj = 0; uj < m_numCovariates; ++uj, info1 += m_numCovariates + 2, info2 += m_numCovariates + 2) {
    info1a = info1;
    info2a = info2;
    for (uk = uj; uk < m_numCovariates; ++uk, ++info1a, info2a += m_numCovariates + 1)
      *info2a = *info1a;
  }
  
  sc1 = m_score;
  m_maxScore = 0;
  for (ui = 0; ui < m_numCovariates + 1; ++ui, ++sc1) {
    if (fabs(*sc1) > m_maxScore)
      m_maxScore = fabs(*sc1);
  }
  
}
// Fit the restricted polytomous logistic model
int CGxEPolytomousDataset::RestrictedPolytomousLogistic(double *_beta, double *_inverseInfo) {
  unsigned int ui, uj;
  unsigned int fc;
  double *b1, *sc, *invInfo;
  bool improve = false;
  double oldMax;
  
  m_maxScore = 1.;
  CalculateRestrictedPolytomousScoreAndInfo(_beta);
  fc = 0;
  while (m_maxScore > 1e-6 && fc < 100) {
    improve = false;
    if (Invert(m_information, m_informationUL, _inverseInfo, m_numParamRestrictedPolytomous, m_numUsed * 1.0e-12) == 0) {
      oldMax = m_maxScore;
      memmove(m_betaNew, _beta, m_numParamRestrictedPolytomous * sizeof(double));
      b1 = m_betaNew;
      invInfo = _inverseInfo;
      for (ui = 0; ui < m_numParamRestrictedPolytomous; ++ui, ++b1) {
        sc = m_score;
        for (uj = 0; uj < m_numParamRestrictedPolytomous; ++uj, ++sc, ++invInfo)
          *b1 += *sc * *invInfo;
      }
      CalculateRestrictedPolytomousScoreAndInfo(m_betaNew);
      if (oldMax > m_maxScore) {
        memmove(_beta, m_betaNew, m_numParamRestrictedPolytomous * sizeof(double));
        improve = true;
      }
      ++fc;
    }
    if (improve == false) {
      // Need to come up routine to handle situation when information cannot be inverted
      return 1;
    }
  }
  if (m_maxScore > 1e-6)
    return 1;
  if (Invert(m_information, m_informationUL, _inverseInfo, m_numParamRestrictedPolytomous, m_numUsed * 1.0e-12))
    return 1;
  
  return 0;
}

// Doubles the inverse information -- need for logistic regression assuming HWE
void CGxEPolytomousDataset::DoubleInverseInformation(double *invInfo, double numCovariates){
  unsigned int ui;
  double nc2;
  double *inverseInfo;
  
  inverseInfo = invInfo;
  nc2 = numCovariates * numCovariates;
  for (ui = 0; ui < nc2; ++ui, ++inverseInfo)
    *inverseInfo /= 2;
  inverseInfo = invInfo;
  nc2 = 1;
}

// Subset the data to either cases or controls -- all subjects are selected when UpdateGene is called
void CGxEPolytomousDataset::SelectCaseControl(bool _case) {
  unsigned int ui;
  const bool *complete;
  const bool *missingGene;
  bool *used;
  const bool *caseControl;
  unsigned int count;
  
  memset(m_used, true, m_numSubjects*sizeof(bool));
  
  complete = m_complete;
  missingGene = m_missingGene;
  caseControl = (bool *)m_outcome;
  used = m_used;
  count = m_numSubjects;
  for (ui = 0; ui < m_numSubjects; ++ui, ++complete, ++missingGene, ++caseControl, ++used) {
    if (*complete == false || *missingGene == true || *caseControl != _case) {
      *used = false;
      --count;
    }
  }
  
  if (_case == true)
    m_numCasesUsed = count;
  else
    m_numControlsUsed = count;
}

// This sets the XtX values after standardizing
// This seems a little odd. Should consider adding standardizing in the AssignCovariateValues routine
void CGxEPolytomousDataset::Standardize() {
  unsigned int ui, uj, uk;
  unsigned int xtxSize;
  double *xtx, *xtxi;
  const double *x1, *x2, *x3;
  const bool *miss;
  CGxEDataset::Standardize();
  
  xtxSize = (m_numCovariates * (m_numCovariates + 1)) / 2;
  std::fill(m_subjectXTX, m_subjectXTX + xtxSize * m_numSubjects, 0);
  xtx = m_subjectXTX;
  xtxi = (xtx + xtxSize) - m_numCovariates;
  x1 = m_covInt;
  miss = m_missingCovariates;
  for (ui = 0; ui < m_numSubjects; ++ui, x1 += m_numParameters, ++miss) {
    if (*miss == true) {
      xtx += xtxSize;
      continue;
    }
    xtxi = (xtx + xtxSize) - m_numCovariates;
    x2 = x1;
    for (uj = 0; uj < m_numCovariates; ++uj, ++x2) {
      if (uj == m_gxeNumber) {
        x3 = x1;
        for (uk = 0; uk < m_numCovariates; ++uk, ++x3) {
          if (uk == m_gxeNumber)
            continue;
          *xtxi = *x2 * *x3;
          ++xtxi;
        }
        *xtxi = *x2 * *x2;
        ++xtxi;
      }
      else {
        x3 = x2;
        for (uk = uj; uk < m_numCovariates; ++uk, ++x3) {
          if (uk == m_gxeNumber)
            continue;
          *xtx = *x2 * *x3;
          ++xtx;
        }
      }
    }
    xtx += m_numCovariates;
  }
}

// Fit all the models
int CGxEPolytomousDataset::FitModels() {
  int retval = 0;
  double zval1 = 0, zval2 = 0;
  
  if (m_geneFrequency < m_minMaf) {
    m_errorString = "Minor allele frequency under minimum value";
    return 0x2fff;
  }
  
  retval = CGxELogisticDataset::FitModels();
  
  memset(m_betaG_E, 0, m_numCovariates * sizeof(double));
  m_betaG_E[0] = log(m_geneFrequency / (1. - m_geneFrequency));
  CalculateScoreConstants(m_gene, m_numCovariates);
  if (Logistic((bool *)m_outcome, m_numCovariates, m_betaG_E, m_inverseInformationG_E)) {
    retval |= 0x0124;
  } else {
    DoubleInverseInformation(m_inverseInformationG_E, m_numCovariates);
    zval1 = fabs(m_betaG_E[m_numCovariates - 1] / sqrt(m_inverseInformationG_E[(m_numCovariates - 1) * (m_numCovariates + 1)]));
    memmove(m_betaCaseOnly, m_betaG_E, m_numCovariates * sizeof(double));
    memmove(m_betaCntlOnly, m_betaG_E, m_numCovariates * sizeof(double));

    if (zval1 > m_geCutoff) {
      if (m_bProbabilities == true || m_bDosages == false) {
        if (m_bProbabilities == true)
          CalculateRestrictedPolytomousDosageScoreConstants();
        else
          CalculateRestrictedPolytomousMeasuredScoreConstants();
        if (m_geneCount[0] > MinimumCellCount && m_geneCount[2] > MinimumCellCount) {
          memset(m_betaRestrictedPolytomousG_E, 0, m_numParamRestrictedPolytomous * sizeof(double));
          memmove(m_betaRestrictedPolytomousG_E, m_betaG_E, (m_numCovariates - 1) * sizeof(double));
          m_betaRestrictedPolytomousG_E[0] += log(2);
          m_betaRestrictedPolytomousG_E[m_numParamRestrictedPolytomous - 2] = m_betaG_E[0] + m_betaG_E[0];
          m_betaRestrictedPolytomousG_E[m_numParamRestrictedPolytomous - 1] = m_betaG_E[m_numCovariates - 1];
          if (RestrictedPolytomousLogistic(m_betaRestrictedPolytomousG_E, m_inverseInformationRestrictedPolytomousG_E) != 0) {
            retval |= 0x0120;
          } else {
            if (m_bProbabilities == true)
              CalculatePolytomousDosageScoreConstants();
            else
              CalculatePolytomousMeasuredScoreConstants();
            memset(m_betaPolytomousG_E, 0, m_numParamPolytomous * sizeof(double));
            memmove(m_betaPolytomousG_E, m_betaRestrictedPolytomousG_E, (m_numParamRestrictedPolytomous - 1) * sizeof(double));
            memmove(m_betaPolytomousG_E + m_numParamRestrictedPolytomous - 1, m_betaRestrictedPolytomousG_E + 1, (m_numParamRestrictedPolytomous - 3) * sizeof(double));
            m_betaPolytomousG_E[m_numParamPolytomous - 1] = m_betaPolytomousG_E[m_numParamRestrictedPolytomous - 1];
            if (PolytomousLogistic(m_betaPolytomousG_E, m_inverseInformationPolytomousG_E) != 0)
              retval |= 0x0020;
          }
        } else {
          retval |= 0x0120;
        }
      } else {
        retval |= 0x0120;
      }
    } else {
      retval |= 0x0120;
    }
  }

  // Case-only HW model
  SelectCaseControl(true);
  CalculateScoreConstants(m_gene, m_numCovariates, true);
  if (Logistic((bool *)m_outcome, m_numCovariates, m_betaCaseOnly, m_inverseInformationCaseOnly)) {
    retval |= 0x0248;
  } else {
    DoubleInverseInformation(m_inverseInformationCaseOnly, m_numCovariates);
    zval1 = fabs(m_betaCaseOnly[m_numCovariates - 1] / sqrt(m_inverseInformationCaseOnly[(m_numCovariates - 1) * (m_numCovariates + 1)]));
  }
  
  // Control only models
  SelectCaseControl(false);
  CalculateScoreConstants(m_gene, m_numCovariates, false);
  if (Logistic((bool *)m_outcome, m_numCovariates, m_betaCntlOnly, m_inverseInformationCntlOnly)) {
    retval |= 0x0490;
  } else {
    DoubleInverseInformation(m_inverseInformationCntlOnly, m_numCovariates);
    zval2 = fabs(m_betaCntlOnly[m_numCovariates - 1] / sqrt(m_inverseInformationCntlOnly[(m_numCovariates - 1) * (m_numCovariates + 1)]));
  }
  
  SelectCaseControl(true);
  if (((retval & 0x0008) == 0) && (zval1 > m_geCutoff || zval2 > m_geCutoff)) {
    if (m_bProbabilities == true || m_bDosages == false) {
      if (m_bProbabilities == true)
        CalculateRestrictedPolytomousDosageScoreConstants();
      else
        CalculateRestrictedPolytomousMeasuredScoreConstants();
      if (m_geneCount[0] > MinimumCellCount && m_geneCount[2] > MinimumCellCount) {
        memset(m_betaRestrictedPolytomousCaseOnly, 0, m_numParamRestrictedPolytomous * sizeof(double));
        memmove(m_betaRestrictedPolytomousCaseOnly, m_betaCaseOnly, (m_numCovariates - 1) * sizeof(double));
        m_betaRestrictedPolytomousCaseOnly[0] += log(2);
        m_betaRestrictedPolytomousCaseOnly[m_numParamRestrictedPolytomous - 2] = m_betaCaseOnly[0] + m_betaCaseOnly[0];
        m_betaRestrictedPolytomousCaseOnly[m_numParamRestrictedPolytomous - 1] = m_betaCaseOnly[m_numCovariates - 1];
        if (RestrictedPolytomousLogistic(m_betaRestrictedPolytomousCaseOnly, m_inverseInformationRestrictedPolytomousCaseOnly) != 0) {
          retval |= 0x0240;
        } else {
          if (m_bProbabilities == true)
            CalculatePolytomousDosageScoreConstants();
          else
            CalculatePolytomousMeasuredScoreConstants();
          memset(m_betaPolytomousCaseOnly, 0, m_numParamPolytomous * sizeof(double));
          memmove(m_betaPolytomousCaseOnly, m_betaRestrictedPolytomousCaseOnly, (m_numParamRestrictedPolytomous - 1) * sizeof(double));
          memmove(m_betaPolytomousCaseOnly + m_numParamRestrictedPolytomous - 1, m_betaRestrictedPolytomousCaseOnly + 1, (m_numParamRestrictedPolytomous - 3) * sizeof(double));
          m_betaPolytomousCaseOnly[m_numParamPolytomous - 1] = m_betaPolytomousCaseOnly[m_numParamRestrictedPolytomous - 1];
          if (PolytomousLogistic(m_betaPolytomousCaseOnly, m_inverseInformationPolytomousCaseOnly) != 0)
            retval |= 0x0040;
        }
      } else {
        retval |= 0x0240;
      }
    } else {
      retval |= 0x0240;
    }
  } else {
    retval |= 0x240;
  }
  
  // Control only models
  SelectCaseControl(false);
  if (((retval & 0x0010) == 0) && (zval1 > m_geCutoff || zval2 > m_geCutoff)) {
    if (m_bProbabilities == true || m_bDosages == false) {
      if (m_bProbabilities == true)
        CalculateRestrictedPolytomousDosageScoreConstants();
      else
        CalculateRestrictedPolytomousMeasuredScoreConstants();
      if (m_geneCount[0] > MinimumCellCount && m_geneCount[2] > MinimumCellCount) {
        memset(m_betaRestrictedPolytomousControlOnly, 0, m_numParamRestrictedPolytomous * sizeof(double));
        memmove(m_betaRestrictedPolytomousControlOnly, m_betaCntlOnly, (m_numCovariates - 1) * sizeof(double));
        m_betaRestrictedPolytomousControlOnly[0] += log(2);
        m_betaRestrictedPolytomousControlOnly[m_numParamRestrictedPolytomous - 2] = m_betaCntlOnly[0] + m_betaCntlOnly[0];
        m_betaRestrictedPolytomousControlOnly[m_numParamRestrictedPolytomous - 1] = m_betaCntlOnly[m_numCovariates - 1];
        if (RestrictedPolytomousLogistic(m_betaRestrictedPolytomousControlOnly, m_inverseInformationRestrictedPolytomousControlOnly) != 0) {
          retval |= 0x0480;
        } else {
          if (m_bProbabilities == true)
            CalculatePolytomousDosageScoreConstants();
          else
            CalculatePolytomousMeasuredScoreConstants();
          memset(m_betaPolytomousControlOnly, 0, m_numParamPolytomous * sizeof(double));
          memmove(m_betaPolytomousControlOnly, m_betaRestrictedPolytomousControlOnly, (m_numParamRestrictedPolytomous - 1) * sizeof(double));
          memmove(m_betaPolytomousControlOnly + m_numParamRestrictedPolytomous - 1, m_betaRestrictedPolytomousControlOnly + 1, (m_numParamRestrictedPolytomous - 3) * sizeof(double));
          m_betaPolytomousControlOnly[m_numParamPolytomous - 1] = m_betaPolytomousControlOnly[m_numParamRestrictedPolytomous - 1];
          if (PolytomousLogistic(m_betaPolytomousControlOnly, m_inverseInformationPolytomousControlOnly) != 0)
            retval |= 0x0080;
        }
      } else {
        retval |= 0x0480;
      }
    } else {
      retval |= 0x0480;
    }
  } else {
    retval |= 0x0480;
  }
//  retval |= 0x080;
/*  
  if (m_bProbabilities == true) {
    CalculatePolytomousDosageScoreConstants();
  }
  else if (m_bDosages == false) {
    CalculatePolytomousMeasuredScoreConstants();
  }
  if (m_geneCount[0] > MinimumCellCount && m_geneCount[2] > MinimumCellCount) {
    if ((retval & 0x20) != 0) {
      retval |= 0x80;
    }
    else {
      memset(m_betaPolytomousControlOnly, 0, m_numParamPolytomous * sizeof(double));
      if (PolytomousLogistic(m_betaPolytomousControlOnly, m_inverseInformationPolytomousControlOnly) != 0)
        retval |= 0x80;
    }
    if ((retval & 0x0100) != 0) {
      retval |= 0x400;
    }
    else {
      if (m_bProbabilities == true)
        CalculateRestrictedPolytomousDosageScoreConstants();
      else if (m_bDosages == false)
        CalculateRestrictedPolytomousMeasuredScoreConstants();
      memset(m_betaRestrictedPolytomousControlOnly, 0, m_numParamRestrictedPolytomous * sizeof(double));
      if (RestrictedPolytomousLogistic(m_betaRestrictedPolytomousControlOnly, m_inverseInformationRestrictedPolytomousControlOnly) != 0)
        retval |= 0x0400;
    }
  }
  else {
    retval |= 0x0480;
  }
*/
  return retval;
}
/*
// ******************************************** mplink linear regression ********************************************
// This does all the linear regression, assigning memory for all possible tests in mplink
// There should be a class for just a single linear regression model.
// This should be a subclass of a single linear regresssion class. ???
// ********************************************************************************************************************

// *************************************************************************** //
//                     Contstructor/Destructor                                 //
// *************************************************************************** //

// Constructor

CGxELinearDataset::CGxELinearDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
                                     const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
                                     const bool *_filter) : CGxEDataset(_numSub, _outcome, _missingOutcome, _numCov, _cov, _missingCov, _gxeNumber, _filter) {
#ifdef DATASET_DEBUG
  std::cout << "Entering GxELinearDataset constructor" << std::endl;
#endif
  
  m_betaD_E = NULL;
  m_betaD_GE = NULL;
  m_betaD_GxE = NULL;
  m_betaLevene = NULL;
  
  m_inverseInformationD_E = NULL;
  m_inverseInformationD_GE = NULL;
  m_inverseInformationD_GxE = NULL;
  m_inverseInformationLevene = NULL;
  
  m_residuals = NULL;
  m_xty = NULL;
  
  m_Levene = NULL;
}
// Destructor
CGxELinearDataset::~CGxELinearDataset() {
#ifdef DEBUG_DELETE
  std::cout << "Entering GxELogisticDataset Destructor" << std::endl;
#endif
  if (m_residuals)
    delete[] m_residuals;
  if (m_xty)
    delete[] m_xty;
  if (m_Levene)
    delete m_Levene;
}

// *************************************************************************** //
//                        Memory management                                    //
// *************************************************************************** //

// Clear out allocated memory
void CGxELinearDataset::DeallocateMemory() {
#ifdef DATASET_DEBUG
  std::cout << "Deallocating GxELinearDataset memory" << std::endl;
#endif
  if (m_residuals)
    delete[] m_residuals;
  if (m_xty)
    delete[] m_xty;
  
  m_residuals = NULL;
  m_xty = NULL;
}
// Allocate memory needed for betas
// Need to consider using GxEDataset routine
bool CGxELinearDataset::AllocateBetaMemory() {
  m_beta = new double[m_numParameters + m_numCovariates + m_numCovariates + 2];
  if (m_beta == NULL)
    return false;
  m_betaD_E = m_beta;
  m_betaD_GE = m_betaD_E + m_numCovariates;
  m_betaD_GxE = m_betaD_GE + m_numCovariates + 1;
  m_betaLevene = m_betaD_GxE + m_numParameters;
  
  return true;
}
// Allocate memory needed for score test - currently score not used
bool CGxELinearDataset::AllocateScoreMemory() {
  return true;
}
// Allocate memory needed for information and inverse information
// Again, I need to consider using GxEDataset routine
bool CGxELinearDataset::AllocateInformationMemory() {
  m_information = new double[m_numParameters * m_numParameters];
  m_informationUL = new double[m_numParameters * m_numParameters];
  m_inverseInformation = new double[2 * m_numCovariates * m_numCovariates + 2 * m_numCovariates + m_numParameters * m_numParameters + 2];
  
  if (m_information == NULL || m_inverseInformation == NULL || m_informationUL == NULL)
    return false;
  
  m_inverseInformationD_E = m_inverseInformation;
  m_inverseInformationD_GE = m_inverseInformationD_E + m_numCovariates * m_numCovariates;
  m_inverseInformationD_GxE = m_inverseInformationD_GE + m_numCovariates * m_numCovariates + 2 * m_numCovariates + 1;
  m_inverseInformationLevene = m_inverseInformationD_GxE + m_numParameters * m_numParameters;
  
  return true;
}
// Allocate all memory needed to do analyses
bool CGxELinearDataset::AllocateMemory() {
#ifdef DATASET_DEBUG
  std::cout << "Allocating linear dataset memory" << std::endl;
#endif
  if (CGxEDataset::AllocateMemory() == false)
    return false;
  
  m_residuals = new double[m_numSubjects];
  if (m_residuals == NULL)
    return false;
  
  m_xty = new double[m_numParameters];
  if (m_xty == NULL)
    return false;
  
  return true;
}

// *************************************************************************** //
//                         Data management                                     //
// *************************************************************************** //

// Update the gene
void CGxELinearDataset::UpdateGene() {
  CGxEDataset::UpdateGene();
  XTransposeY(m_covInt, (double *)m_outcome, m_xty, m_used, m_numParameters, m_numSubjects);
  XTransposeX(m_covInt, m_xtransposex, m_used, m_numParameters, m_numSubjects);
}

// *************************************************************************** //
//                           Data testing                                      //
// *************************************************************************** //

// Test data to see if convergence is possible
void CGxELinearDataset::RunTests() {
  unsigned int ui;
  double *residual;
  const double *x1;
  const double *y;
  const bool *used;
  double denom;
  double alpha, beta;
  double m1, m2, m3;
  
  CompletenessTest();
  if (m_good == false)
    return;
  
  LinearDependenceTest();
  if (m_good == false)
    return;
  
  InitializeBeta();
  if (m_good == false)
    return;
  
  memmove(m_used, m_complete, m_numSubjects * sizeof(bool));
  m_numUsed = m_numCompleteData;
  
  XTransposeY(m_covInt, (double *)m_outcome, m_xty, m_used, m_numParameters, m_numSubjects);
  XTransposeX(m_covInt, m_xtransposex, m_used, m_numParameters, m_numSubjects);
  
  m1 = *(m_xtransposex + (m_numParameters - 1) * (m_numParameters - 1) - 1);
  m2 = *(m_xtransposex + (m_numParameters - 1) * (m_numParameters - 1));
  m3 = *(m_xtransposex + m_numParameters * m_numParameters - 1);
  denom = m1*m3 - m2*m2;
  alpha = (m3*m_xty[m_numParameters - 2] - m2*m_xty[m_numParameters - 1]) / denom;
  beta = (-m2*m_xty[m_numParameters - 2] + m1*m_xty[m_numParameters - 1]) / denom;
  
  if (Linear(m_numCovariates, m_betaD_E, m_inverseInformationD_E, m_s2[0])) {
    m_errorString = "Unable to maximum model with no genes in model";
    m_good = false;
    return;
  }
  
  x1 = m_covInt + m_numParameters - 1;
  used = m_used;
  y = (double *)m_outcome;
  residual = m_residuals;
  for (ui = 0; ui < m_numSubjects; ++ui, x1 += m_numParameters, ++y, ++residual, ++used) {
    if (*used == false) {
      *residual = 1e10;
      continue;
    }
    *residual = *y - alpha - beta * *x1;
  }
  
  if (m_Levene != NULL)
    delete m_Levene;
  m_Levene = NULL;
}

// *************************************************************************** //
//                       Maximization routines                                 //
// *************************************************************************** //

// Select measured of dosage Levene's test
void CGxELinearDataset::SetLevene(bool _dosages, bool _probabilities) {
  m_bDosages = _dosages;
  m_bProbabilities = _probabilities;
  
  if (m_Levene != NULL)
    delete m_Levene;
  m_Levene = NULL;
  
  if (m_bDosages == true) {
    if (m_bProbabilities == true)
      m_Levene = new CLeveneDosageMedian(m_residuals, m_used, m_numSubjects);
    else
      m_Levene = NULL;
  }
  else {
    m_Levene = new CLeveneGeneMeasuredMedian(m_residuals, m_used, m_numSubjects);
  }
}
// Linear regression
int CGxELinearDataset::Linear(unsigned int numCovar, double *beta, double *inverseInfo, double &s2) {
  double *ii;
  double *b;
  const double *x1, *x2, *xtx, *y;
  double *info, *xty;
  bool *used;
  double residual;
  unsigned int ui, uj;
  
  memset(beta, 0, numCovar * sizeof(double));
  memset(m_information, 0, m_numParameters * m_numParameters * sizeof(double));
  memset(inverseInfo, 0, numCovar * numCovar * sizeof(double));
  
  xtx = m_xtransposex;
  info = m_information;
  for (ui = 0; ui < numCovar; ++ui, info += numCovar, xtx += m_numParameters)
    memmove(info, xtx, numCovar * sizeof(double));
  
  if (Invert(m_information, m_informationUL, inverseInfo, numCovar, m_numUsed * 1.0e-12)) {
    return 1;
  }
  
  b = beta;
  ii = inverseInfo;
  for (ui = 0; ui < numCovar; ++ui, ++b) {
    xty = m_xty;
    for (uj = 0; uj < numCovar; ++uj, ++ii, ++xty)
      *b += *xty * *ii;
  }
  
  s2 = 0;
  x1 = m_covInt;
  used = m_used;
  y = (double *)m_outcome;
  for (ui = 0; ui < m_numSubjects; ++ui, x1 += m_numParameters, ++y, ++used) {
    if (*used == false)
      continue;
    x2 = x1;
    residual = *y;
    b = beta;
    for (uj = 0; uj < numCovar; ++uj, ++x2, ++b)
      residual -= *b * *x2;
    s2 += residual * residual;
  }
  s2 /= (m_numUsed - numCovar);
  
  return 0;
}
// Fit all the models
int CGxELinearDataset::FitModels() {
  int retval;
  
  retval = 0;
  if (m_initialized == false) {
    m_errorString = "Model not initialized";
    return -1;
  }
  
  if (m_numUsed == 0)
    return 7;
  
  // Should be numCovariates + 1, + numParameters, the extra 1 is for Levene' test
  // This should be removed
  //	memset(m_betaD_GE, 0, (m_numCovariates + m_numParameters + 2) * sizeof(double));
  //	memmove(m_betaD_GE, m_betaD_E, m_numCovariates*sizeof(double));
  
  if (Linear(m_numCovariates + 1, m_betaD_GE, m_inverseInformationD_GE, m_s2[1]))
    retval += 1;
  
  if (Linear(m_numParameters, m_betaD_GxE, m_inverseInformationD_GxE, m_s2[2]))
    retval += 2;
  
  if (m_Levene == NULL || m_Levene->Levene(m_gene))
    retval += 4;
  
  return retval;
}
*/

