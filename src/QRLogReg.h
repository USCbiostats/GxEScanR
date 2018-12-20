#ifndef QRLOGREG
#define QRLOGREG 1

class CQRLogReg {
protected:
  int m_numRow, m_numCol;
  
  arma::vec m_abx, m_expabx, m_expabxp1, m_expitabx;
  arma::vec m_w;
  arma::vec m_wInv;
  arma::mat m_ql, m_rtl;
  arma::vec m_yp0;
  arma::vec m_zt0;
  arma::vec m_k0;
  arma::vec m_beta0;
  
  arma::vec m_yp;
  arma::vec m_zt;
  arma::vec m_k;
  arma::mat m_rtr;
  arma::mat m_t;
  
  arma::mat m_qr;
  arma::mat m_rbr;
  arma::mat m_h;
  arma::vec m_zb;
  
  arma::vec m_bb, m_bt;
  arma::vec m_beta;
  arma::mat m_xr;
  arma::vec m_score, m_score0;
  
  arma::mat m_xl;
  arma::vec m_y;
  
  double m_logLike;
public:
  CQRLogReg(int numRow, int numCol, arma::vec &y, arma::mat &x, arma::vec &beta0, std::vector<std::vector<double> > &space);
  ~CQRLogReg() {}
  
  int FitModel(arma::mat &xr);
  //  int FitModel2();
  
  const arma::vec &Beta() const { return m_beta; }
  const arma::mat &Score() const { return m_score; }
  double LogLike() const { return m_logLike; }
};

#endif
