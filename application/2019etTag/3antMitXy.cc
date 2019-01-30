#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <valarray>
#include <algorithm>
#include <sstream>
#include <utility>

#include "../../eiBauBeDi.h"

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
    virtual double operator()(EiBauBeDi::tree) = 0;
    virtual ~wgFun() {};
  };

  class gFun : public wgFun {
  public:
    double operator()(EiBauBeDi::tree);
  };
  double gFun::operator()(EiBauBeDi::tree x) {return pow(x.d/200.,2) * M_PI * x.wgt;}

  class vFun : public wgFun {
  public:
    double operator()(EiBauBeDi::tree);
  };
  double vFun::operator()(EiBauBeDi::tree x) {return pow(x.d/200.,2) * M_PI * x.h/2. * x.wgt;}

  //define influence function
  double funCircle(const double &px, const double &py, const EiBauBeDi::tree &tree) {
    double weight = -1.;
    double dx = fabs(tree.x - px);
    double dy = fabs(tree.y - py);
    if(max(dx, dy) < tree.influence0) {
      double distance = sqrt(pow(dx, 2) + pow(dy, 2));
      if(distance <= tree.influence0) {
	weight = 1. / (pow(tree.influence0, 2) * M_PI);
	if(distance == tree.influence0) {weight /= 2.;}
      }
    }
    return(weight);
  }

  double funCircleWgt(const double &px, const double &py, const EiBauBeDi::tree &tree) {
    double weight = -1.;
    double dx = fabs(tree.x - px);
    double dy = fabs(tree.y - py);
    if(max(dx, dy) < tree.influence0) {
      double distance = sqrt(pow(dx, 2) + pow(dy, 2));
      if(distance <= tree.influence0) {
	weight = 1. / (pow(tree.influence0, 2) * M_PI);
	weight *= 2. * (1. - pow(distance, 2)/pow(tree.influence0, 2));
      }
    }
    return(weight);
  }

  double funCircleWgt2(const double &px, const double &py, const EiBauBeDi::tree &tree) {
    double weight = 0.;
    double dx = tree.x - px;
    double dy = tree.y - py;
    if(max(dx, dy) < tree.influence0) {
      double distance = sqrt(pow(dx, 2) + pow(dy, 2));
      if(distance < tree.influence0) {
	weight = tree.h * (1. - pow(distance / tree.influence0, tree.influence1));
      }
    }
    return(weight);
  }

}

int main(int argc, char** argv) {
  if(argc < 3) {
    cout << "Usage: " << argv[0] << " filenameTrees filenameCorners [function]" << endl;
    exit(0);
  }
  
  wgFun *wg;
  if(argc > 3) {
    switch(atoi(argv[3])) {
    case 1 : wg = new (vFun);
      break;
    default : wg = new (gFun);
    }
  } else {wg = new (gFun);}
  
  //Read in the plot corners
  std::vector<std::pair<double, double> > plotCorners;
  {
    vector<array<double, 3> > tmpPlotCorners; //cornerNumber x y
    ifstream infile(argv[2]);
    if (!infile.is_open()) {
      std::cout << "failed to open " << argv[2] << '\n';
      exit(0);
    }
    string line;
    while (getline(infile, line)) {
      if(line[0] == '#') continue;
      istringstream iss(line);
      double cornerNumber, x, y;
      if (!(iss >> cornerNumber >> x >> y)) {break;} //Read error
      tmpPlotCorners.push_back(array<double, 3>{cornerNumber, x, y});
    }
    sort(tmpPlotCorners.begin(), tmpPlotCorners.end());
    for(auto&& i : tmpPlotCorners) {plotCorners.push_back(make_pair(i[1],i[2]));}
  }
  EiBauBeDi::forestStand stand(plotCorners);
  
  //Read in the tree data
  {
    ifstream infile(argv[1]);
    if (!infile.is_open()) {
      std::cout << "failed to open " << argv[1] << '\n';
      exit(0);
    }
    string line;
    while (getline(infile, line)) {
      if(line[0] == '#') continue;
      istringstream iss(line);
      vector<string> tokens{istream_iterator<string>{iss},
	  istream_iterator<string>{}};
      //Use only trees which are inside the plot or on the border
      bool in = stand.poly.pointInPolyCN(stod(tokens[1]),stod(tokens[2]));
      bool on = stand.poly.pointOnPoly(stod(tokens[1]),stod(tokens[2]));
      if(in || on) {
	EiBauBeDi::tree tree;
	tree.nr = tokens[0];
	tree.x = stod(tokens[1]);
	tree.y = stod(tokens[2]);
	tree.z = stod(tokens[3]);
	tree.wgt = stod(tokens[4]);
	tree.d = stod(tokens[6]);
	tree.h = stod(tokens[7]);
	if(tree.d > 0.) {tree.d += 1.3;
	} else {tree.d += tree.h > 1.3 ? 1.3 : tree.h;}
	stand.trees.push_back(tree);
      }
    }
  }

  //Set the tree impact
  for(auto&& i : stand.trees) {i.impact = (*wg)(i);}
  
  //Homogenous Basal area on the plot
  cout << "#homogenous Basal area on the plot" << endl;
  cout << "#type tree density growingArea" << endl;
  {
    double density = stand.getImpactSum() / stand.poly.getArea();
    for(auto&& i : stand.trees) {
      double growingArea = i.impact / density;
      cout << "hba " << i.nr << " " << density * 10000. << " " << growingArea << "\n";
    }
  }

  //Fix sample plot at the position of each tree
  cout << "\n#Basal area from fixed sample plot around each tree" << endl;
  //Set influence0 to 7 (influence radius)
  for(auto&& i : stand.trees) {
    i.influence0 = 7.;
    if(stand.poly.shareInside(i.x, i.y, i.influence0) <= 0.) {
      cout << "#Sample plot size of " << i.influence0 << " is to large for tree " << i.nr << endl;;
    }
  }
  cout << "#type tree density growingArea" << endl;
  for(auto&& i : stand.trees) {
    double density = stand.subsamplePoint(i.x, i.y, funCircle, true);
    double growingArea = i.impact / density;
    cout << "fspt " << i.nr << " " << density * 10000. << " " << growingArea << "\n";
  }

  //Fix sample plot at the position of each tree with size of tree-height
  cout << "\n#Basal area from fixed sample plot = h around each tree" << endl;
  //Set influence0 to h (influence radius)
  for(auto&& i : stand.trees) {
    i.influence0 = i.h;
    if(stand.poly.shareInside(i.x, i.y, i.influence0) <= 0.) {
      cout << "#Sample plot size of " << i.influence0 << " is to large for tree " << i.nr << endl;;
    }
  }
  cout << "#type tree density growingArea" << endl;
  for(auto&& i : stand.trees) {
    double density = stand.subsamplePoint(i.x, i.y, funCircle, true);
    double growingArea = i.impact / density;
    cout << "vspht " << i.nr << " " << density * 10000. << " " << growingArea << "\n";
  }

  //Angle count at the position of each tree
  cout << "\n#Basal area from angle count sample around each tree" << endl;
  //Set influence0 to d/4
  for(auto&& i : stand.trees) {
    i.influence0 = i.d/4.;
    if(stand.poly.shareInside(i.x, i.y, i.influence0) <= 0.) {
      cout << "#Sample plot size of " << i.influence0 << " is to large for tree " << i.nr << endl;;
    }
  }
  cout << "#type tree density growingArea" << endl;
  for(auto&& i : stand.trees) {
    double density = stand.subsamplePoint(i.x, i.y, funCircle, true);
    double growingArea = i.impact / density;
    cout << "wzpt " << i.nr << " " << density * 10000. << " " << growingArea << "\n";
  }
  //Weighted angle count at the position of each tree
  cout << "\n#Basal area from variable angle count sample around each tree" << endl;
  cout << "#type tree density growingArea" << endl;
  for(auto&& i : stand.trees) {
    double density = stand.subsamplePoint(i.x, i.y, funCircleWgt, true);
    double growingArea = i.impact / density;
    cout << "wwzpt " << i.nr << " " << density * 10000. << " " << growingArea << "\n";
  }

  //Overlapping circles at the position of each tree
  cout << "\n#Basal area from overlaping circles around each tree" << endl;
  cout << "\n#Circle size = 100*d/4" << endl;
  {
    //Set influence0 to d/4 -> impacton this area is 4m2/ha
    double maxInf = 0.;
    for(auto&& i : stand.trees) {
      i.influence0 = i.d/4.;
      i.impact = (*wg)(i) / (pow(i.influence0,2) * M_PI);
      if(maxInf < i.influence0) {maxInf = i.influence0;}
    }
    for(auto&& i : stand.trees) {
      if(stand.poly.shareInside(i.x, i.y, i.influence0 + maxInf) <= 0.) {
	cout << "#Sample plot size of " << i.influence0 << " is to large for tree " << i.nr << endl;;
      }
    }
  }
  cout << "#type tree density growingArea" << endl;
  for(auto&& i : stand.trees) {
    double density = stand.subsampleCircle(i.x, i.y, i.d/4., true);
    double growingArea = (*wg)(i) / density;
    cout << "ovlap " << i.nr << " " << density  * 10000. << " " << growingArea << "\n";
  }

  {
    //Fix sample plot at systematic rasterpoints
    cout << "\n#Basal area from fixed sample plot on raster" << endl;
    cout << "#type tree density growingArea" << endl;
    double r = 7.;
    for(auto&& i : stand.trees) {i.impact = (*wg)(i);}
    for(auto&& i : stand.trees) {i.influence0 = r;}
    std::array<double, 4> ext = stand.poly.getExtends();
    double dxy = 0.1; //Distance between raster point
    std::valarray<double> tmp(0., stand.trees.size());
    std::valarray<double> density(0., stand.trees.size());
    std::valarray<unsigned int> count(0U, stand.trees.size());
    for(double x = ext[0] - r; x <= ext[1] + r; x += dxy) {
      for(double y = ext[2] - r; y <= ext[3] + r; y += dxy) {
	tmp = stand.subsamplePointTree(x, y, funCircle, false);
	density += tmp;
	if(stand.poly.pointInPolyCN(x, y)) {
	  count[tmp>0.] = valarray(count[tmp>0.]) + 1U;
	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      if(count[i] > 0.) {density[i] /= count[i];
      } else {density[i] = 0.;}
      double growingArea = stand.trees[i].impact / density[i];
      cout << "fspr " << stand.trees[i].nr << " " << density[i] * 10000. << " " << growingArea << "\n";
    }
    //The same but sample only inside the stand
    cout << "\n#Basal area from fixed sample plot on raster inside stand" <<endl;
    cout << "#type tree density growingArea" << endl;
    density = 0.;
    count = 0U;
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.subsamplePointTree(x, y, funCircle, true);
	  density += tmp;
	  count[tmp>0.] = valarray(count[tmp>0.]) + 1U;
 	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      if(count[i] > 0.) {density[i] /= count[i];
      } else {density[i] = 0.;}
      double growingArea = stand.trees[i].impact / density[i];
      cout << "fspri " << stand.trees[i].nr << " " << density[i] * 10000. << " " << growingArea << "\n";
    }
    //Weighted with the distance between tree and sample center
    cout << "\n#Basal area from weightened fixed sample plot on raster" << endl;
    cout << "#type tree density growingArea" << endl;
    density = 0.;
    count = 0U;
    for(double x = ext[0] - r; x <= ext[1] + r; x += dxy) {
      for(double y = ext[2] - r; y <= ext[3] + r; y += dxy) {
	tmp = stand.subsamplePointTree(x, y, funCircleWgt, false);
	density += tmp;
	if(stand.poly.pointInPolyCN(x, y)) {
	  count[tmp>0.] = valarray(count[tmp>0.]) + 1U;
	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      if(count[i] > 0.) {density[i] /= count[i];
      } else {density[i] = 0.;}
      double growingArea = stand.trees[i].impact / density[i];
      cout << "fsprw " << stand.trees[i].nr << " " << density[i] * 10000. << " " << growingArea << "\n";
    }
  }

  //assign plot area to trees - Winner takes it all
  {
    cout << "\n#Voronoi" << endl;
    cout << "#type tree density growingArea" << endl;
    for(auto&& i : stand.trees) {i.impact = 99.;}
    for(auto&& i : stand.trees) {i.influence0 = 99.;}
    double dxy = 0.1; //Distance between raster point
    std::array<double, 4> ext = stand.poly.getExtends();
    ext[0] -= dxy/2.; ext[1] += dxy/2.;
    ext[2] -= dxy/2.; ext[3] += dxy/2.;
    std::valarray<double> tmp(0., stand.trees.size());
    std::valarray<unsigned int> count(0U, stand.trees.size());
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.influencePoint(x, y, funCircleWgt, false);
	  if(tmp.sum() > 0.) {
	    ++count[distance(begin(tmp), max_element(begin(tmp), end(tmp)))];
	  }
	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      double growingArea = count[i]*dxy*dxy;
      double density = 0.;
      if(growingArea > 0.) {density = (*wg)(stand.trees[i]) / (growingArea / 10000.);}
      cout << "wtaVor " << stand.trees[i].nr << " " << density << " " << growingArea << "\n";
    }

    cout << "\n#Standflaeche WTA Kegel" << endl;
    cout << "#type tree density growingArea" << endl;
    count=0;
    for(auto&& i : stand.trees) {i.impact = (*wg)(i);}
    for(auto&& i : stand.trees) {i.influence0 = 100.*sqrt(i.impact / M_PI) / 2.;}
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.influencePoint(x, y, funCircleWgt, false);
	  if(tmp.sum() > 0.) {
	    ++count[distance(begin(tmp), max_element(begin(tmp), end(tmp)))];
	  }
	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      double growingArea = count[i]*dxy*dxy;
      double density = 0.;
      if(growingArea > 0.) {density = stand.trees[i].impact / (growingArea / 10000.);}
      cout << "wtaK " << stand.trees[i].nr << " " << density << " " << growingArea << "\n";
    }

    //Assign plot area to trees - Share pixel between trees
    cout << "\n#Standflaeche Share Kegel" << endl;
    cout << "#type tree density growingArea" << endl;
    for(auto&& i : stand.trees) {i.impact = (*wg)(i);}
    for(auto&& i : stand.trees) {i.influence0 = 100.*sqrt(i.impact / M_PI) / 2.;}
    std::valarray<double> area(0., stand.trees.size());
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.influencePoint(x, y, funCircleWgt, false);
	  double total = tmp.sum();
	  if(total > 0.) {area += tmp/total;}
	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      double growingArea = area[i]*dxy*dxy;
      double density = 0.;
      if(growingArea > 0.) {density = stand.trees[i].impact / (growingArea / 10000.);}
      cout << "shK " << stand.trees[i].nr << " " << density << " " << growingArea << "\n";
    }

    //Assign plot area to trees - Share pixel between trees
    cout << "\n#Standflaeche Share Krone" << endl;
    cout << "#type tree density growingArea" << endl;
    area = 0;
    for(auto&& i : stand.trees) {i.impact = 1.;}
    for(auto&& i : stand.trees) {i.influence0 = i.d/4.;}
    for(auto&& i : stand.trees) {i.influence1 = 0.5;}
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.influencePoint(x, y, funCircleWgt2, false);
	  double tmpMax = tmp.max();
	  if(tmpMax > 0.) {
	    tmp /= tmpMax;
	    tmp = pow(tmp, 2.);
	    double total = tmp.sum();
	    if(total > 0.) {area += tmp/total;}
	  }
	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      double growingArea = area[i]*dxy*dxy;
      double density = 0.;
      if(growingArea > 0.) {density = (*wg)(stand.trees[i]) / (growingArea / 10000.);}
      cout << "shC " << stand.trees[i].nr << " " << density << " " << growingArea << "\n";
    }
    
  }
  
  return 0;
}
