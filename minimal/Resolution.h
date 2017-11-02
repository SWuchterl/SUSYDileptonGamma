#include <sstream>
#include <string>

#include <iostream>
using std::cout;
using std::cin;
#include <fstream>
using std::ifstream;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <iterator>
using std::istream_iterator;
#include <algorithm>
using std::copy;

class Resolution {
 public:
  Resolution(){};
  ~Resolution(){};
  Resolution(const string& filename);
  float get(float pt, float eta, float rho);
 private:
  std::vector<float> numbers;
};

Resolution::Resolution(const string& filename)
{
  ifstream myfile(filename);
  copy(istream_iterator<float>(myfile), istream_iterator<float>(), back_inserter(numbers));
}

float Resolution::get(float pt, float eta, float rho)
{
  float res = -1;
  if (eta<-4.7) eta= -4.69999;
  if (eta>4.7) eta = 4.69999;
  if (rho>44.31) rho = 44.3;
  for (unsigned int i=0;i<numbers.size();i += 11) {
    if (numbers.at(i) < eta && eta < numbers.at(i+1) && numbers.at(i+2) < rho && rho < numbers.at(i+3)) {
      if (pt<numbers.at(i+5)) pt = numbers.at(i+5);
      if (pt>numbers.at(i+6)) pt = numbers.at(i+6);
      float a=numbers.at(i+7), b=numbers.at(i+8), c=numbers.at(i+9), d=numbers.at(i+10);
      res = sqrt(a*abs(a)/(pt*pt)+b*b*pow(pt,d)+c*c);
      break;
    }
  }
  if (res<0) { cout << "Resolution::get could not find resolution for pt = " << pt << " eta = " << eta << " rho = " << rho << endl; }
  return res;
}
