#include <deque>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;

bool read_rules(string filename, 
		vector<string>& symbols,
		vector<pair<size_t, vector<pair<size_t,size_t> > > >& rules) {
  
  // Initial structures
  set<string> sorted_symbols;
  map<string,size_t> invert_symbols;
  vector<pair<string, multiset<string> > > raw_rules;

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

    // Grab a rule
    string symbol;
    parser >> symbol;

    if (!parser.eof()) {
      symbols.push_back(symbol);
      raw_rules.push_back(pair<string,multiset<string> >(symbol,multiset<string>()));
      parser >> symbol;
    }
    
    while (!parser.fail()) {
      symbols.push_back(symbol);
      raw_rules[raw_rules.size()-1].second.insert(symbol);
      parser >> symbol;
    }
    
    getline(infile,line);
  }

  sort(symbols.begin(),symbols.end());
  vector<string>::iterator sym_itr = unique(symbols.begin(),symbols.end());
  symbols.resize(sym_itr - symbols.begin());

  // Put symbols into structures in order
  for (int x = 0; x < symbols.size(); x++)
    invert_symbols.insert(pair<string,size_t>(symbols[x],symbols.size()-1));

  cout << "Symbols: " << symbols.size() << endl;

  for (size_t x = 0; x < raw_rules.size(); x++) {
    cout << "Rule " << x << ": " << raw_rules[x].first << " - ";

    rules.push_back(pair<size_t,vector<pair<size_t,size_t> > >(invert_symbols.find(raw_rules[x].first)->second,vector<pair<size_t,size_t> >()));
    
    for (size_t y = 0; y < symbols.size(); y++) {
      size_t count = raw_rules[x].second.count(symbols[y]);
      cout << count << "*" << symbols[y] << " ";
      if (count > 0)
	rules[x].second.push_back(pair<size_t,size_t>(y,count));
    }
    cout << endl;
  }
  
  if (symbols.size() == 0 || rules.size() == 0)
    return false;
  
  return true;
}

void print_rules(ostream& os,
		 vector<string>& symbols,
		 vector<pair<size_t, vector<pair<size_t,size_t> > > >& rules) {

  for (size_t x = 0; x < rules.size(); x++) {
    cout << symbols[rules[x].first] << " ";
    for (size_t y = 0; y < rules[x].second.size(); y++) {
      for (size_t w = 0; w < rules[x].second[y].second; w++)
	cout << symbols[rules[x].second[y].first] << " ";
    }
    cout << endl;
  }

  return;
}

int main(int argc,
	 char* argv[]) {

  vector<string> symbols;
  vector<pair<size_t, vector<pair<size_t,size_t> > > > rules;

  cout << "Reading:" << endl;
  if (argc > 1)
    cout << read_rules(string(argv[1]),symbols,rules) << endl;

  cout << "Writing:" << endl;
  print_rules(cout,symbols,rules);
  
  return 0;
}
