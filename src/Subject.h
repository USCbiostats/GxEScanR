#ifndef SUBJECT_H
#define SUBJECT_H 1

int TestSubjectData(Rcpp::DataFrame &subjectData);

class CSubject {
protected:
  std::string m_fid;
  std::string m_iid;
  
  void print(std::ostream &_stream) const { _stream << m_fid << '\t' << m_iid; }
public:
  CSubject() { m_fid = ""; m_iid = ""; }
  CSubject(const CSubject &_src) { m_fid = _src.m_fid; m_iid = _src.m_iid; }
  CSubject(const std::string &_fid, const std::string &_iid) { m_fid = _fid; m_iid = _iid; }
  ~CSubject() {}
  
  CSubject &operator=(const CSubject &_rhs);
  bool operator==(const CSubject &_rhs) const;
  
  CSubject &Assign(const std::string &_fid, const std::string &_iid) {
    m_fid = _fid;
    m_iid = _iid;
    return *this;
  }
  void FamilyID(const std::string &_fid) { m_fid = _fid; }
  void SubjectID(const std::string &_iid) { m_iid = _iid; }
  const std::string &FamilyID() const { return m_fid; }
  const std::string &SubjectID() const { return m_iid; }
  
  friend std::ostream &operator<<(std::ostream &_stream, const CSubject &_rhs) {
    _rhs.print(_stream);
    return _stream;
  }
};

#endif