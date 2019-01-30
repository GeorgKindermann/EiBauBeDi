#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <valarray>
#include <algorithm>
#include <memory>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

using namespace std;

namespace {
  struct sTree {
    string tree;
    double gew;
    double d;
    double h;
  };

  class wgFun {
  public:
    virtual double operator()(sTree) = 0;
    virtual ~wgFun() {};
  };

  class gFun : public wgFun {
  public:
    double operator()(sTree);
  };
  double gFun::operator()(sTree x) {return pow(x.d/200.,2) * M_PI * x.gew;}

  class vFun : public wgFun {
  public:
    double operator()(sTree);
  };
  double vFun::operator()(sTree x) {return pow(x.d/200.,2) * M_PI * x.h/2. * x.gew;}

}

int main(int argc, char** argv) {
  if(argc < 3) {
    cout << "Usage: " << argv[0] << " filename plotsize [minG/ha] [function]" << endl;
    exit(0);
  }
  vector<sTree> dat;
  const double plotSize = atof(argv[2]);
  const double minGHa = argc > 3 ? atof(argv[3]) : 0.;
  
  unique_ptr<wgFun> wg;
  if(argc > 4) {
    switch(atoi(argv[4])) {
    case 1 : wg = make_unique<vFun>();
      break;
    default : wg = make_unique<gFun>();
    }
  } else {wg = make_unique<gFun>();}

  {
    ifstream istrm(argv[1]);
    if (!istrm.is_open()) {
      std::cout << "failed to open " << argv[1] << '\n';
      exit(0);
    }
    string line;
    getline(istrm,line);
    sTree x;
    string discard;
    while(istrm >> x.tree >> discard >> discard >> discard >> x.gew >> discard >> x.d >> x.h >> discard >> discard) {
      if(x.d > 0.) {x.d += 1.3;
      } else {x.d += x.h > 1.3 ? 1.3 : x.h;}
      dat.push_back(x);
    }
  }
  double sumTreeInfluence = 0.;
  for(auto const& tree: dat) {
    sumTreeInfluence += (*wg)(tree);
  }
  valarray<double> share(dat.size());
  valarray<bool> fix(dat.size());
  for(std::vector<sTree>::size_type i = 0; i != dat.size(); ++i) {
    share[i] = (*wg)(dat[i]) / sumTreeInfluence;
    fix[i] = false;
  }
  bool updateShare = true;
  double freeShare = 0.;
  while(updateShare) {
    updateShare = false;
    for(std::vector<sTree>::size_type i = 0; i != dat.size(); ++i) {
      double growingArea = plotSize * share[i];
      double density = dat[i].gew * pow(dat[i].d/2., 2) * M_PI / growingArea;
      if(density < minGHa) {
	updateShare = true;
	fix[i] = true;
	freeShare += share[i] * (1. - density / minGHa);
	share[i] *= density / minGHa;
      }
    }
    if(updateShare && freeShare > 0.) {
      double varShare = valarray(share[!fix]).sum();
      share[!fix] = valarray(share[!fix]) * (varShare + freeShare) / varShare;
      freeShare = 0.;
    }
  }
  cout << "#tree g/ha growingArea" << endl;
  for(std::vector<sTree>::size_type i = 0; i != dat.size(); ++i) {
    double growingArea = plotSize * share[i];
    double density = dat[i].gew * pow(dat[i].d/2., 2) * M_PI / growingArea;
    cout << dat[i].tree << " "
	 << " " << std::scientific << density
	 << " " << std::scientific << growingArea
	 << "\n";
  }
  
  return 0;
}
