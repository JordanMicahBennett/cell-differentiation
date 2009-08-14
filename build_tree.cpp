#include <deque>
#include <vector>
#include <set>
#include <multiset>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

bool read_rules(string filename, 
		vector<string>& symbols,
		vector<vector<size_t> >& rules) {

  // Initial structures
  set<string> sorted_symbols;
  vector<multiset<string> > raw_rules;

  // Open a file
  ifstream infile;
  infile.open(filename.c_str());
  
  if (infile.eof() || infile.bad())
    return false;

  // Read each line
  string line;
  getline(infile,line);
  while (!infile.eof()) {

    // Parse by spaces
    stringstream parser;
    parser << line;
    
    string symbol;
    symbol << parser;

    if (!parser.eof())
      raw_rules.push_back(multiset<string>)

    while (!parser.eof()) {
      sorted_symbols.insert(symbol);
      raw_rules[raw_rules.length()-1].insert(symbol);
    }
			    
    getline(infile,line)
  }

  for (set<string>::iterator itr = sorted_symbols.begin();
       itr != sorted_symbols.begin();
       
  }

  for (size_t x = 0; x < raw_rules.length
  return true;
}

int main(int argc, 
	 char* argv[]) {

}
