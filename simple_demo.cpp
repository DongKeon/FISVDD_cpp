#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

#include "fisvdd.cpp"

using namespace std;

vector<vector<double>> read_csv(string filename){
  // Reads a CSV file into a vector of <string, vector<int>> pairs where
  // each pair represents <column name, column values>

  // Create a vector of <string, int vector> pairs to store the result
  vector<vector<double>> result;

  // Create an input filestream
  ifstream myFile(filename);

  // Make sure the file is open
  if(!myFile.is_open()) throw runtime_error("Could not open file");

  // Helper vars
  string line, colname;
  double val;

  // Read the column names
  if(myFile.good())
  {
    // Extract the first line in the file
    getline(myFile, line);

    // Create a stringstream from line
    stringstream ss(line);

    // Extract each column name
    while(getline(ss, colname, ',')){}
  }

  // Read data, line by line
  while(getline(myFile, line))
  {
    // Create a stringstream of the current line
    stringstream ss(line);
    // Keep track of the current column index
    vector<double> row_vec;
    // Extract each integer
    while(ss >> val){
      row_vec.push_back(val);
      
      // If the next token is a comma, ignore it and move on
      if(ss.peek() == ',') ss.ignore();
      
      // Increment the column index
    }
    result.push_back(row_vec);
  }

  // Close file
  myFile.close();

  return result;
}


int main() {
  vector<vector<double>> data = read_csv("sample_input.csv");

  tuple<vector<double>, vector<vector<double>>> fisvdd_output;
  vector<double> alpha; vector<vector<double>> sv;
  // initialization with a random number between 0~1
  double s = 0.8;

  fisvdd *model = new fisvdd(data[0], s);
  fisvdd_output = model->find_sv(data);

  alpha = get<0>(fisvdd_output); sv = get<1>(fisvdd_output);
  printf("alpha -------\n");
  print_1d_vector(alpha);
  printf("\nsupport vector -------\n");
  print_2d_vector(sv);

}
