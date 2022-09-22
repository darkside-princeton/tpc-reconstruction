#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Config.h"

using namespace std;

Config *Config::config = NULL;

int Config::Load(){
  for(int i_f = 0; i_f < int(_config_file_name_vec.size()); i_f++){
    string config_file_name = _config_file_name_vec[i_f];
    cout << config_file_name << endl;
    ifstream cfile(config_file_name.c_str());
    if(!cfile.is_open()){
      cout << "ERROR: Could not open file " << config_file_name << endl;
      return 1;
    } // if !cfile.is_open()                                                                                                                                             

    string buff;
    while(getline(cfile,buff)){
      if(buff == "") continue;
      stringstream ss(buff);
      string type, pname, val_str, other;
      ss >> type >> pname >> val_str >> other;
      if(type[0] == '#' || pname[0] == '#' || val_str[0] == '#') continue;

      stringstream val_ss(val_str);
      if(type == "int"){
        int val = 0;
        val_ss >> val;
        SetParameter(pname,val);
      }
      else if(type == "double"){
        double val = 0;
        val_ss >> val;
        SetParameter(pname,val);
      }
      else if(type == "string"){
        string val = "";
        val_ss >> val;
        SetParameter(pname,val);
      }
      else if(type == "vector"){
        string parname = val_str;
        if(other != "("){
          cout << "ERROR: In config file, vector " << parname << " must enumerate items between parentheses surrounded by spaces" << endl;
          exit(0);
        }
        type = pname;
	if(type == "int"){
          string newval_str = "(";
          vector<int> val_vec;
          while(newval_str != ")"){
            ss >> newval_str;
            if(newval_str == ")") break;
            stringstream newval_ss(newval_str);
            int val = 0;
            newval_ss >> val;
            val_vec.push_back(val);
          } // while newval_str != )
          SetParameter(parname, val_vec);
        } // if type == int 
	else if(type == "double"){
          string newval_str = "(";
          vector<double> val_vec;
          while(newval_str != ")"){
            ss >> newval_str;
            if(newval_str == ")") break;
            stringstream newval_ss(newval_str);
            double val = 0;
            newval_ss >> val;
            val_vec.push_back(val);
          } 
          SetParameter(parname, val_vec);
        } // if type == double            
	else{
          cout << "WARNING: Unrecognized vector type for line: " << buff << endl;
        }
      }
      else{
        cout << "WARNING: Unrecognized type for line: " << buff << endl;
      }
    }
  }
  return 0;
}

