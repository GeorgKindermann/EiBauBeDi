// eiBauBeDi 0.1a, EInzelBAUmBEstandesDIchte - Single Tree Stand Density
// Programm to calculte the stand density of singe trees in a forest stand
// Copyright (C) 2017-2019 Georg Kindermann
// Home: https://github.com/GeorgKindermann/EiBauBeDi

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3 of the License.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <valarray>
#include <chrono>

#include "eiBauBeDi.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

using namespace std;

//define influence function
double funCircle(const double &px, const double &py, const EiBauBeDi::tree &tree) {
  double weight = -1.;
  double dx = tree.x - px;
  double dy = tree.y - py;
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
  double dx = tree.x - px;
  double dy = tree.y - py;
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
double fun0Rast(const EiBauBeDi::tree &ctree, const double &dist2, const double &maxDist2, const double &c1, const double &c2) {
  return(pow(ctree.h * (1. - pow(dist2 / maxDist2, c1)),c2));
  //return(ctree.h * (1. - dist2 / maxDist2));
}
double fun0RastWta(const EiBauBeDi::tree &ctree __attribute__((unused)), const double &dist2, const double &maxDist2, const double &c1 __attribute__((unused)), const double &c2 __attribute__((unused))) {
  return(dist2 / maxDist2);
}
  

int main(int argc, char *argv[]) {
  string FileInCorners = "./data/corners.txt";
  string FileInTrees = "./data/trees.txt";
  string FileOutTrees = "./result/standDensity.txt";
  if(argc == 4) {
    FileInCorners = argv[1];
    FileInTrees = argv[2];
    FileOutTrees = argv[3];
  }

  //Read in the plot corners
  std::vector<std::pair<double, double> > plotCorners;
  {
    vector<array<double, 3> > tmpPlotCorners; //cornerNumber x y
    ifstream infile(FileInCorners);
    string line;
    while (getline(infile, line)) {
      if(line[0] == '#') continue;
      istringstream iss(line);
      string plotNr;
      double cornerNumber, x, y;
      if (!(iss >> plotNr >> cornerNumber >> x >> y)) {break;} //Read error
      tmpPlotCorners.push_back(array<double, 3>{cornerNumber, x, y});
    }
    sort(tmpPlotCorners.begin(), tmpPlotCorners.end());
    for(auto&& i : tmpPlotCorners) {plotCorners.push_back(make_pair(i[1],i[2]));}
  }
  EiBauBeDi::forestStand stand(plotCorners);
  
  //Read in the tree data
  {
    ifstream infile(FileInTrees);
    string line;
    while (getline(infile, line)) {
      if(line[0] == '#') continue;
      istringstream iss(line);
      vector<string> tokens{istream_iterator<string>{iss},
	  istream_iterator<string>{}};
      //Use only trees which are inside the plot
      if(stand.poly.pointInPolyCN(stod(tokens[2]),stod(tokens[3]))) {
	EiBauBeDi::tree tree;
	tree.nr = tokens[1];
	tree.x = stod(tokens[2]);
	tree.y = stod(tokens[3]);
	tree.z = stod(tokens[4]);
	tree.d = stod(tokens[6]);
	tree.h = stod(tokens[7]);
	tree.hcr = stod(tokens[8]);
	stand.trees.push_back(tree);
      }
    }
  }

  
  //Set the tree impact to its basal area
  for(auto&& i : stand.trees) {i.impact = pow(i.d/2.,2) * M_PI;}
  //Basal area / hectare of the plot
  cout << "Basal area [m2/ha]: " << stand.getImpactSum() / stand.poly.getArea() << endl;

  //Fix sample plot at the position of each tree
  cout << "\nBasal area from fixed sample plot around each tree" << endl;
  cout << "tree g/ha" << endl;
  //Set influence0 to 7 (influence radius)
  for(auto&& i : stand.trees) {i.influence0 = 7.;}
  for(auto&& i : stand.trees) {
    double gha = stand.subsamplePoint(i.x, i.y, funCircle, true);
    cout << i.nr << " " << gha << endl;
  }

  //Angle count at the position of each tree
  cout << "\nBasal area from angle count sample around each tree" << endl;
  cout << "tree g/ha" << endl;
  //Set influence0 to d/4
  for(auto&& i : stand.trees) {i.influence0 = i.d/4.;}
  for(auto&& i : stand.trees) {
    double gha = stand.subsamplePoint(i.x, i.y, funCircle, true);
    cout << i.nr << " " << gha << endl;
  }

  //Weighted angle count at the position of each tree
  cout << "\nBasal area from variable angle count sample around each tree" << endl;
  cout << "tree g/ha" << endl;
  for(auto&& i : stand.trees) {
    double gha = stand.subsamplePoint(i.x, i.y, funCircleWgt, true);
    cout << i.nr << " " << gha << endl;
  }

  //Overlapping circles at the position of each tree
  cout << "\nBasal area from overlaping circles around each tree" << endl;
  cout << "tree g/ha" << endl;
  for(auto&& i : stand.trees) {i.impact = 4.;}
  for(auto&& i : stand.trees) {
    double gha = stand.subsampleCircle(i.x, i.y, i.d/4., true);
    cout << i.nr << " " << gha << endl;
  }

  //Fix sample plot at systematic rasterpoints
  cout << "\nBasal area from fixed sample plot on raster" << endl;
  cout << "tree g/ha" << endl;
  {
    double r = 7.;
    for(auto&& i : stand.trees) {i.impact = pow(i.d/2.,2) * M_PI;}
    for(auto&& i : stand.trees) {i.influence0 = r;}
    std::array<double, 4> ext = stand.poly.getExtends();
    double dxy = 0.1; //Distance between raster point
    std::valarray<double> tmp(0., stand.trees.size());
    std::valarray<double> gha(0., stand.trees.size());
    std::valarray<unsigned int> count(0U, stand.trees.size());
    for(double x = ext[0] - r; x <= ext[1] + r; x += dxy) {
      for(double y = ext[2] - r; y <= ext[3] + r; y += dxy) {
	tmp = stand.subsamplePointTree(x, y, funCircle, false);
	gha += tmp;
	for(size_t i = 0; i<tmp.size(); ++i) {if(tmp[i] > 0.) {++count[i];}}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      if(count[i] > 0.) {gha[i] /= count[i];
      } else {gha[i] = 0.;}
      cout << stand.trees[i].nr << " " << gha[i] << endl;
    }
    //The same but sample only inside the stand
    cout << "\nBasal area from fixed sample plot on raster inside stand" <<endl;
    gha = 0.;
    count = 0U;
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.subsamplePointTree(x, y, funCircle, true);
	  gha += tmp;
	  for(size_t i = 0; i<tmp.size(); ++i) {if(tmp[i] > 0.) {++count[i];}}
 	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      if(count[i] > 0.) {gha[i] /= count[i];
      } else {gha[i] = 0.;}
      cout << stand.trees[i].nr << " " << gha[i] << endl;
    }
    //Weighted with the distance between tree and sample center
    cout << "\nBasal area from weightened fixed sample plot on raster" << endl;
    cout << "tree g/ha" << endl;
    gha = 0.;
    count = 0U;
    for(double x = ext[0] - r; x <= ext[1] + r; x += dxy) {
      for(double y = ext[2] - r; y <= ext[3] + r; y += dxy) {
	tmp = stand.subsamplePointTree(x, y, funCircleWgt, false);
	gha += tmp;
	for(size_t i = 0; i<tmp.size(); ++i) {if(tmp[i] > 0.) {++count[i];}}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      if(count[i] > 0.) {gha[i] /= count[i];
      } else {gha[i] = 0.;}
      cout << stand.trees[i].nr << " " << gha[i] << endl;
    }
  }

  //assign plot area to trees - Winner takes it all
  cout << "\nBasal area from single tree growing area on raster" << endl;
  cout << "tree g/ha distanceCenter directionCenter uncircularity axis" << endl;
  {
    for(auto&& i : stand.trees) {i.impact = pow(i.d/2.,2) * M_PI;}
    for(auto&& i : stand.trees) {i.influence0 = i.d/4.;}
    double dxy = 0.1; //Distance between raster point
    std::array<double, 4> ext = stand.poly.getExtends();
    ext[0] -= dxy/2.; ext[1] += dxy/2.;
    ext[2] -= dxy/2.; ext[3] += dxy/2.;
    std::valarray<double> tmp(0., stand.trees.size());
    std::valarray<unsigned int> count(0U, stand.trees.size());
    valarray<double> sumX(0., stand.trees.size());
    valarray<double> sumY(0., stand.trees.size());
    valarray<double> sumD(0., stand.trees.size());
    valarray<double> sumXY(0., stand.trees.size());
    valarray<double> sumX2(0., stand.trees.size());
    valarray<double> sumY2(0., stand.trees.size());
    valarray<double> comp(0., stand.trees.size()*stand.trees.size());
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.influencePoint(x, y, funCircleWgt, false);
	  auto result = std::max_element(begin(tmp), end(tmp));
	  size_t idx = distance(begin(tmp), result);
	  ++count[idx];
	  sumX[idx] += x; sumY[idx] += y;
	  comp[slice(idx*stand.trees.size(), stand.trees.size(), 1)] += tmp;
	}
      }
    }
    valarray<double> treeX(0., stand.trees.size());
    valarray<double> treeY(0., stand.trees.size());
    for(size_t i=0; i<count.size(); ++i) { //Schwerpunkt der Standflaeche
      treeX[i] = sumX[i] / count[i];
      treeY[i] = sumY[i] / count[i];
    }
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.influencePoint(x, y, funCircleWgt, false);
	  auto result = std::max_element(begin(tmp), end(tmp));
	  size_t idx = distance(begin(tmp), result);
	  sumD[idx] += sqrt(pow(treeX[idx] - x, 2) + pow(treeY[idx] - y, 2));
	  sumXY[idx] += (x-treeX[idx]) * (y - treeY[idx]);
	  sumX2[idx] += pow(x-treeX[idx], 2);
	  sumY2[idx] += pow(y-treeY[idx], 2);
	}
      }
    }
    for(size_t i = 0; i<tmp.size(); ++i) {
      sumX[i] /= count[i]; sumY[i] /= count[i]; sumD[i] /= count[i];
      double gha = 0.;
      if(count[i] > 0.) {gha = stand.trees[i].impact / (count[i]*dxy*dxy);}
      cout << stand.trees[i].nr << " " << gha;
      cout << " " << sqrt(pow(sumX[i] - stand.trees[i].x,2) + pow(sumY[i] - stand.trees[i].y,2));
      cout << " " << atan2(sumY[i] - stand.trees[i].y, sumX[i] - stand.trees[i].x);
      cout << " " << sumD[i] / (2./3. * sqrt(count[i]*dxy*dxy) / M_PI);
      cout << " " << (sumXY[i] / sumX2[i] + sumY2[i] / sumXY[i]) / 2.;
      cout << endl;
    }
    cout << "tree Competitor..." << endl;
    cout << "X";
    for(size_t i = 0; i<tmp.size(); ++i) {cout << " " << stand.trees[i].nr;}
    cout << endl;
    for(size_t i = 0; i<tmp.size(); ++i) {
      cout << stand.trees[i].nr;
      valarray<double> competitors = comp[slice(i*stand.trees.size(), stand.trees.size(), 1)];
      //competitors[i] = 0.; //In case the tree by itself should not be included
      competitors /= competitors.sum();
      for(size_t j = 0; j<stand.trees.size(); ++j) {
	cout << " " << competitors[j];
      }
      cout << endl;
    }
  }
  
  //Assign plot area to trees - Share pixel between trees
  cout << "\nBasal area from single tree growing shared area on raster" << endl;
  cout << "tree g/ha" << endl;
  {
    for(auto&& i : stand.trees) {i.impact = pow(i.d/2.,2) * M_PI;}
    for(auto&& i : stand.trees) {i.influence0 = i.d/4.;}
    std::array<double, 4> ext = stand.poly.getExtends();
    double dxy = 0.1; //Distance between raster point
    std::valarray<double> tmp(0., stand.trees.size());
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
      double gha = 0.;
      if(area[i] > 0.) {gha = stand.trees[i].impact / (area[i]*dxy*dxy);}
      cout << stand.trees[i].nr << " " << gha << endl;
    }
  }

  //Rasterize example
  cout << "\nBasal area from single tree growing shared area using rasterize" << endl;
  cout << "tree g/ha" << endl;
  for(auto&& i : stand.trees) {i.impact = pow(i.d/2.,2) * M_PI;}
  for(auto&& i : stand.trees) {i.influence0 = i.d/4.;}
  {
    std::array<double, 4> ext = stand.poly.getExtends();
    std::valarray<double> area(0., stand.trees.size());
    area = stand.rasterize(ext[0],ext[1],ext[2],ext[3], 0.1, 0.1, 0.5, 2, fun0Rast);
    for(size_t i = 0; i<area.size(); ++i) {
      double gha = 0.;
      if(area[i] > 0.) {gha = stand.trees[i].impact / area[i];}
      cout << stand.trees[i].nr << " " << gha << endl;
    }
  }


  //Speed comparrison
  cout << "\nBasal area from single tree growing shared area on raster" << endl;
  cout << "tree area0 area1 area2" << endl;
   {
    for(auto&& i : stand.trees) {i.impact = 1.;}
    for(auto&& i : stand.trees) {i.influence0 = i.d/4.;}
    for(auto&& i : stand.trees) {i.influence1 = 0.5;}
    std::array<double, 4> ext = stand.poly.getExtends();
    double dxy = 0.1; //Distance between raster point
    auto t0 = std::chrono::high_resolution_clock::now();
    std::valarray<double> tmp(0., stand.trees.size());
    std::valarray<double> area(0., stand.trees.size());
    for(double x = ext[0]; x <= ext[1]; x += dxy) {
      for(double y = ext[2]; y <= ext[3]; y += dxy) {
	if(stand.poly.pointInPolyCN(x, y)) {
	  tmp = stand.influencePoint(x, y, funCircleWgt2, false);
	  tmp /= tmp.max();
	  tmp = pow(tmp, 2.);
	  double total = tmp.sum();
	  if(total > 0.) {area += tmp/total;}
	}
      }
    }
    area *= dxy*dxy;
    auto t1 = std::chrono::high_resolution_clock::now();
    std::valarray<double> area1(0., stand.trees.size());
    area1 = stand.rasterize(ext[0],ext[1],ext[2],ext[3], dxy, dxy, 0.25, 2, fun0Rast);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::valarray<double> area2(0., stand.trees.size());
    area2 = stand.rasterizeWta(ext[0],ext[1],ext[2],ext[3], dxy, dxy, 1, 1, fun0RastWta);
    auto t3 = std::chrono::high_resolution_clock::now();
    double sum0=0.; double sum1=0.; double sum2=0.;
    for(size_t i = 0; i<tmp.size(); ++i) {
      cout << stand.trees[i].nr << " " << area[i] << " " << area1[i] << " " << area2[i] << endl;
      sum0 += area[i]; sum1 += area1[i]; sum2 += area2[i];
    }
    std::chrono::duration<double> diff0 = t1 - t0;
    std::chrono::duration<double> diff1 = t2 - t1;
    std::chrono::duration<double> diff2 = t3 - t2;
    cout << "Duration: " << diff0.count() << " " << diff1.count() << " " << diff2.count() << endl;
    cout << "Areasum: " << sum0 << " " << sum1 << " " << sum2 << endl;
  }


  return(0);
}
