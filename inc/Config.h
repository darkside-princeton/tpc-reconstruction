#ifndef Config_h
#define Config_h

#include <iostream>
#include <vector>
#include <map>

using namespace std;

typedef map<string, int> paramsI;
typedef map<string, double> paramsD;
typedef map<string, string> paramsS;
typedef map<string, vector<int> > paramsIvec;
typedef map<string, vector<double> > paramsDvec;

class Config {

  static Config *config;

  Config() {}

  Config(string fname){
    _config_file_name_vec.push_back(fname);
    Load();
  }
  
  int Load();

  vector<string> _config_file_name_vec;
  paramsI _parsI;
  paramsD _parsD;
  paramsS _parsS;
  paramsIvec     _parsIvec;
  paramsDvec     _parsDvec;


 public:

  static Config *Get() {
    if(!config)
      config = new Config;
    return config;
  }

  int AddFile(string fname){
    _config_file_name_vec.push_back(fname);
    return Load();
  }

  int SetParameter(string pname, int val) {
    bool overwritten = _parsI.find(pname) != _parsI.end();
    _parsI[pname] = val;
    return int(overwritten);
  }
  int SetParameter(string pname, double val) {
    bool overwritten = _parsD.find(pname) != _parsD.end();
    _parsD[pname] = val;
    return int(overwritten);
  }

  int SetParameter(string pname, string val) {
    bool overwritten = _parsS.find(pname) != _parsS.end();
    _parsS[pname] = val;
    return int(overwritten);
  }

  int SetParameter(string pname, vector<int> val) {
    bool overwritten = _parsIvec.find(pname) != _parsIvec.end();
    _parsIvec[pname] = val;
    return int(overwritten);
  }

  int SetParameter(string pname, vector<double> val) {
    bool overwritten = _parsDvec.find(pname) != _parsDvec.end();
    _parsDvec[pname] = val;
    return int(overwritten);
  }

  int GetParameterI(string pname){
    bool valid = _parsI.find(pname) != _parsI.end();
    if(!valid){
      cout << "ERROR: Invalid parameter " << pname << endl;
      return -999;
    }
    return _parsI[pname];
  }
  double GetParameterD(string pname){
    bool valid = _parsD.find(pname) != _parsD.end();
    if(!valid){
      cout << "ERROR: Invalid parameter " << pname << endl;
      return -999;
    }
    return _parsD[pname];
  }

  string GetParameterS(string pname){
    bool valid = _parsS.find(pname) != _parsS.end();
    if(!valid){
      cout << "ERROR: Invalid parameter " << pname << endl;
      return "";
    }
    return _parsS[pname];
  }
  vector<int> GetParameterIvec(string pname){
    bool valid = _parsIvec.find(pname) != _parsIvec.end();
    if(!valid){
      cout << "ERROR: Invalid parameter " << pname << endl;
      vector<int> emptyvec;
      return emptyvec;
    }
    return _parsIvec[pname];
  }
  vector<double> GetParameterDvec(string pname){
    bool valid = _parsDvec.find(pname) != _parsDvec.end();
    if(!valid){
      cout << "ERROR: Invalid parameter " << pname << endl;
      vector<double> emptyvec;
      return emptyvec;
    }
    return _parsDvec[pname];
  }


};

#endif
