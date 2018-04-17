#ifndef GXEDATASET_H
#define GXEDATASET_H 1

#ifndef MATRIXFUNCTIONS_H
#include "MatrixFunctions.h"
#endif
/*
#ifndef LEVENE_H
#include "Levene.h"
#endif
*/
class CGxEDataset {
protected:
  // Has the dataset been initialized?
  bool m_initialized;
  // Does the dataset have no errors?
  bool m_good;
  // Error message
  std::string m_errorString;
  
  // Number of subjects (records) in dataset
  unsigned int m_numSubjects;
  // Number of covariates in dataset - including intercept
  unsigned int m_numCovariates;
  // Number of genetic terms in model - removed -- should always be 1 ??? 
  // unsigned int m_numInteractions;
  // Number of parameters in D|G,E,GxE model
  unsigned int m_numParameters;
  // Maximum number of parameters in model
  unsigned int m_maxParameters;
  // GxE covariate number
  unsigned int m_gxeNumber;
  
  // Outcome values
  // Outcome - void because can be binary or continuous
  const void *m_outcome;
  // Missing outcome value
  const bool *m_missingOutcome;
  
  // Source list of covariates
  // Number of covariates in source list -- not needed since all are used ???
  // unsigned int m_numSourceCovariates;
  // Array of source covariates
  const double *m_sourceCovariates;
  // Array of missing source covariates
  const bool *m_sourceMissingCovariates;
  // Array of covariates to be used -- removed since all covariates must be used ??? may change
  // const bool *m_usedSourceCovariates;
  
  // Covariate values - a covariate can be in both groups
  // Covariates not in interactions
  // double *m_covariates;
  // Covariates in interaction
  // double *m_interactions; 
  // All covariates and interaction values
  double *m_covInt;
  // Missing one or more covariate values
  bool *m_missingCovariates;
  // Mean of the covariates
  double *m_mean;
  // Standard deviation of the covariates
  double *m_std;
  
  // Genetic values
  const double *m_gene;
  // Missing indicator for genetic values
  const bool *m_missingGene;
  // Gene frequency in subjects with complete data
  double m_geneFrequency;
  // Indicator if alleles 1 and 2 were swapped - This occurs when the allele frequency is greater than 0.5
  bool m_allelesSwapped;
  // Is this dosage data
  bool m_bDosages;
  // Are the dosage probabilities available
  bool m_bProbabilities;
  // Excluded from analysis -- filter, used for case-only and control-only analyses
  const bool *m_include;
  
  // Record has complete data and is included in the analysis
  bool *m_complete;
  // Number of subjects with complete data
  unsigned int m_numCompleteData;
  
  // Indicator if subject is used in the analysis
  bool *m_used;
  // Number of subjects used in the analysis
  unsigned int m_numUsed;
  
  // X transpose X
  double *m_xtransposex;
  
  // beta values - array contains values for all model
  double *m_beta;
  // Proposed values for beta
  double *m_betaNew;
  // Score values - not kept test to test because they should be 0
  double *m_score;
  // Information matrix - Only one copy needed - not kept test to test
  double *m_information;
  // Upper lower decompostion of information matrix - needed to calculate inverse
  double *m_informationUL;
  // Inverse information array - contains inverse information array for each model - needed for testing
  double *m_inverseInformation;
  
  virtual bool AllocateBetaMemory() = 0;
  virtual bool AllocateScoreMemory();
  virtual bool AllocateInformationMemory();
  virtual void DeallocateMemory();
  virtual bool AllocateMemory();
  
  // Copy the covariate vales
  virtual void AssignCovariateValues();
  // Assign complete data indicators for subjects
  virtual void AssignComplete();
  // Standard the covariates
  virtual void Standardize();
  
  // Make sure at least one subject has complete data
  void CompletenessTest();
  // Test if covariates are linearly independent
  int LinearDependenceTest();
  // Run tests if data is valid
  virtual void RunTests();
  
  virtual int InitializeBeta() { return 0; }
public:
  CGxEDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
              const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
              const bool *_filter);
  virtual ~CGxEDataset();
  
  bool Initialize();
  virtual int FitModels() = 0;
  
  void AssignGene(const double *_gene, const bool *_missingGene, bool _dosages, bool _probabilities);
  virtual void UpdateGene();
  
  bool isGood() const { return m_good; }
  const std::string ErrorString() const { return m_errorString; }
  unsigned int NumSubjects() const { return m_numSubjects; }
  unsigned int NumCovariates() const { return m_numCovariates; }
  //	unsigned int NumInteractions() const { return m_numInteractions; }
  unsigned int NumParameters() const { return m_numParameters; }
  unsigned int MaxParameters() const { return m_maxParameters; }
  //	const double *Covaraiates() const { return m_covariates; }
  //	const double *Interactions() const { return m_interactions; }
  const double *CompleteCovariates() const { return m_covInt; }
  const double *Mean() const { return m_mean; }
  const double *StandardDeviation() const { return m_std; }
  const bool *MissingCovariates() const { return m_missingCovariates; }
  const void *Outcome() const { return m_outcome; }
  const bool *MissingOutcome() const { return m_missingOutcome; }
  const bool *Included() const { return m_include; }
  const bool *Complete() const { return m_complete; }
  const bool *Used() const { return m_used; }
  unsigned int NumUsed() const { return m_numUsed; }
  unsigned int NumCompleteData() const { return m_numCompleteData; }
  double GeneFrequency() const { return m_geneFrequency; }
  bool AllelesSwapped() const { return m_allelesSwapped; }
  const double *Beta() const { return m_beta; }
  const double *InverseInformation() const { return m_inverseInformation; }
  
  void Print(std::ostream &outfile, int n = 0);
};

class CGxELogisticDataset : public CGxEDataset {
protected:
  unsigned int m_numCasesUsed;
  unsigned int m_numControlsUsed;
  
  // Number of parameters in model
  unsigned int m_numParamD_E;
  unsigned int m_numParamD_G;
  unsigned int m_numParamD_GxE;
  // Number of parameters squared for each model - used to allocate memory for inverse information
  unsigned int m_numParamSquaredD_E;
  unsigned int m_numParamSquaredD_G;
  unsigned int m_numParamSquaredD_GxE;
  // Number of betas summed over all models - used for allocating memory
  unsigned int m_betaSize;
  // Number of betas summed ove all model less the number of betas for the D|E model. Used for zeroing.
  unsigned int m_betaSizeLessDE;
  // Sum of the squares of the number of parmaters for each model - size needed to store all the inverse information matrices
  unsigned int m_invInfoSize;
  // Size needed to store all inverse information matrices less the D|E model. Used for zeroing.
  unsigned int m_invInfoSizeLessDE;
  
  // pointers to beta location for each model results in m_beta array
  double *m_betaD_E;
  double *m_betaD_GE;
  double *m_betaD_GxE;
  
  // Maximum absolute value of score
  double m_maxScore;
  // The constant part of the score calculation - Score values are set to these values at start of each iteration
  double *m_scoreConstants;
  
  // pointers to inverse information location for each model results in m_inverseInformation array
  double *m_inverseInformationD_E;
  double *m_inverseInformationD_GE;
  double *m_inverseInformationD_GxE;
  
  virtual void DeallocateMemory();
  virtual bool AllocateBetaMemory();
  virtual bool AllocateScoreMemory();
  virtual bool AllocateInformationMemory();
  
  void CaseControlCountTest();
  void FullyDeterminedTest();
  
  virtual void RunTests();
  virtual int InitializeBeta();
  
  void CalculateScoreConstants(const bool *outcome, unsigned int numCovar);
  void CalculateScoreConstants(const double *outcome, unsigned int numCovar);
  void CalculateScoreConstants(const double *outcome, unsigned int numCovar, bool casesOnly);
  void CalculateScoreAndInformation(const bool *outcome, unsigned int numCovar, const double *beta);
  int UpdateOneByOne(const bool *outcome, unsigned int numCovar, double *beta, double *inverseInfo);
  int Logistic(const bool *outcome, unsigned int numCovar, double *beta, double *inverseInfo);
public:
  CGxELogisticDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
                      const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
                      const bool *_filter);
  virtual ~CGxELogisticDataset();
  
  virtual int FitModels();
  
  unsigned int NumCasesUsed() const { return m_numCasesUsed; }
  unsigned int NumControlsUsed() const { return m_numControlsUsed; }
  const double *BetaD_E() const { return m_betaD_E; }
  const double *BetaD_GE() const { return m_betaD_GE; }
  const double *BetaD_GxE() const { return m_betaD_GxE; }
  
  const double *InverseInformationD_E() const { return m_inverseInformationD_E; }
  const double *InverseInformationD_GE() const { return m_inverseInformationD_GE; }
  const double *InverseInformationD_GxE() const { return m_inverseInformationD_GxE; }
};

class CGxEPolytomousDataset : public CGxELogisticDataset {
protected:
  unsigned int m_interactionCovariate;
  double m_geneCount[3];
  double *m_polyScoreConstants;
  double *m_subjectXTX;
  
  // Number of parameters in each model
  unsigned int m_numParamHW;
  unsigned int m_numParamPolytomous;
  unsigned int m_numParamRestrictedPolytomous;
  // Number of parameters squared in each model
  unsigned int m_numParamSquaredHW;
  unsigned int m_numParamSquaredPolytomous;
  unsigned int m_numParamSquaredRestrictedPolytomous;
  
  double *m_betaG_E;
  double *m_betaCaseOnly;
  double *m_betaCntlOnly;
  double *m_betaPolytomousG_E;
  double *m_betaPolytomousCaseOnly;
  double *m_betaPolytomousControlOnly;
  double *m_betaRestrictedPolytomousG_E;
  double *m_betaRestrictedPolytomousCaseOnly;
  double *m_betaRestrictedPolytomousControlOnly;
  
  double *m_inverseInformationG_E;
  double *m_inverseInformationCaseOnly;
  double *m_inverseInformationCntlOnly;
  double *m_inverseInformationPolytomousG_E;
  double *m_inverseInformationPolytomousCaseOnly;
  double *m_inverseInformationPolytomousControlOnly;
  double *m_inverseInformationRestrictedPolytomousG_E;
  double *m_inverseInformationRestrictedPolytomousCaseOnly;
  double *m_inverseInformationRestrictedPolytomousControlOnly;
  
  virtual void DeallocateMemory();
  virtual bool AllocateBetaMemory();
  virtual bool AllocateScoreMemory();
  virtual bool AllocateInformationMemory();
  virtual bool AllocateMemory();
  
  virtual void Standardize();
  void SelectCaseControl(bool _case);
  void CalculatePolytomousMeasuredScoreConstants();
  void CalculatePolytomousDosageScoreConstants();
  void CalculatePolytomousScoreAndInfo(const double *_beta);
  void CalculateRestrictedPolytomousMeasuredScoreConstants();
  void CalculateRestrictedPolytomousDosageScoreConstants();
  void CalculateRestrictedPolytomousScoreAndInfo(const double *_beta);
  void DoubleInverseInformation(double *inverseInfo, double numParam);
  // Polytomous logistic regression when beta2 = 2 beta1 for only the interacting covariate
  int PolytomousLogistic(double *beta, double *inverseInfo);
  int RestrictedPolytomousLogistic(double *beta, double *inverseInfo);
public:
  CGxEPolytomousDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
                        const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
                        const bool *_filter);
  virtual ~CGxEPolytomousDataset();
  
  virtual int FitModels();
  
  const double *BetaG_E() const { return m_betaG_E; }
  const double *BetaCaseOnly() const { return m_betaCaseOnly; }
  const double *BetaCntlOnly() const { return m_betaCntlOnly; }
  const double *BetaPolytomousG_E() const { return m_betaPolytomousG_E; }
  const double *BetaPolytomousCaseOnly() const { return m_betaPolytomousCaseOnly; }
  const double *BetaPolytomousControlOnly() const { return m_betaPolytomousControlOnly; }
  const double *BetaRestrictedPolytomousG_E() const { return m_betaRestrictedPolytomousG_E; }
  const double *BetaRestrictedPolytomousCaseOnly() const { return m_betaRestrictedPolytomousCaseOnly; }
  const double *BetaRestrictedPolytomousControlOnly() const { return m_betaRestrictedPolytomousControlOnly; }
  const double *InverseInformationG_E() const { return m_inverseInformationG_E; }
  const double *InverseInformationCO() const { return m_inverseInformationCaseOnly; }
  const double *InverseInformationCntlOnly() const { return m_inverseInformationCntlOnly; }
  const double *InverseInformationPolytomousG_E() const { return m_inverseInformationPolytomousG_E; }
  const double *InverseInformationPolytomousCaseOnly() const { return m_inverseInformationPolytomousCaseOnly; }
  const double *InverseInformationPolytomousControlOnly() const { return m_inverseInformationPolytomousControlOnly; }
  const double *InverseInformationRestrictedPolytomousG_E() const { return m_inverseInformationRestrictedPolytomousG_E; }
  const double *InverseInformationRestrictedPolytomousCaseOnly() const { return m_inverseInformationRestrictedPolytomousCaseOnly; }
  const double *InverseInformationRestrictedPolytomousControlOnly() const { return m_inverseInformationRestrictedPolytomousControlOnly; }
};
/*
class CGxELinearDataset : public CGxEDataset {
protected:
  // residuals needed for levene's test
  double *m_residuals;
  
  // pointers to beta location for each model results in m_beta array
  double *m_betaD_E;
  double *m_betaD_GE;
  double *m_betaD_GxE;
  double *m_betaLevene;
  
  // inverse information matrices
  double *m_inverseInformationD_E;
  double *m_inverseInformationD_GE;
  double *m_inverseInformationD_GxE;
  double *m_inverseInformationLevene;
  
  // XTY - easier to calculate this up front
  double *m_xty;
  
  // estimate of S squared for the models;
  double m_s2[4];
  
  CLeveneGene *m_Levene;
  
  virtual void DeallocateMemory();
  virtual bool AllocateBetaMemory();
  virtual bool AllocateScoreMemory();
  virtual bool AllocateInformationMemory();
  virtual bool AllocateMemory();
  
  virtual void RunTests();
  
  int Linear(unsigned int numCovar, double *beta, double *inverseInfo, double &s2);
  
public:
  CGxELinearDataset(const unsigned int _numSub, const void *_outcome, const bool *_missingOutcome,
                    const unsigned int _numCov, const double *_cov, const bool *_missingCov, unsigned int _gxeNumber,
                    const bool *_filter);
  virtual ~CGxELinearDataset();
  
  virtual int FitModels();
  void SetLevene(bool _dosages, bool _probabilities);
  virtual void UpdateGene();
  
  const double *BetaD_E() const { return m_betaD_E; }
  const double *BetaD_GE() const { return m_betaD_GE; }
  const double *BetaD_GxE() const { return m_betaD_GxE; }
  const double *InverseInformationD_E() const { return m_inverseInformationD_E; }
  const double *InverseInformationD_GE() const { return m_inverseInformationD_GE; }
  const double *InverseInformationD_GxE() const { return m_inverseInformationD_GxE; }
  const double *S2() const { return m_s2; }
  const double W() const { return m_Levene->W(); }
  const unsigned int DF1() const { return m_Levene->DF1(); }
  const unsigned int DF2() const { return m_Levene->DF2(); }
};
*/
#endif
