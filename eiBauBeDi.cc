// eiBauBeDi 0.1a, EInzelBAUmBEstandesDIchte - Single Tree Stand Density
// Programm to calculte the stand density of singe trees in a forest stand
// Copyright (C) 2017 Georg Kindermann
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
#include <array>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <cmath>
#include <iterator>
#include <valarray>
#include <stack>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

using namespace std;

namespace EiBauBeDi {

  struct point {
    double x{0.};
    double y{0.};
  };
  
  struct tree {
    string nr{""}; //tree indentifier
    double x{0.}; //x Position
    double y{0.}; //y Position
    double z{0.}; //z Position
    double d{0.}; //Diameter
    double h{0.}; //Height
    double hcr{0.}; //Height of the crown base
  };

  double polygonArea(const vector<point>& polygon) {
    double area = 0.;
    vector<point>::const_iterator previous = prev(polygon.end());
    for(vector<point>::const_iterator current = polygon.begin(); current != polygon.end(); ++current) {
      area += previous->x * current->y;
      area -= previous->y * current->x;
      previous = current;
    }
    area /= 2.;
    return(abs(area));
  }

  //tests if a point is Left|On|Right of an infinite line.
  //  >0 for P2 left of the line through P0 and P1, 0 for P2  on the line,
  //  <0 for P2  right of the line
  int isLeft(const point &p0, const point &p1, const point &p2) {
    return ((p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y));}

  //winding number test for a point in a polygon
  //  0 (false) = outside, 1 (true) = inside
  bool pointInPolyWN(const point &pkt, const vector<point> &polygon) {
    int nw = 0; //Number of windings
    vector<point>::const_iterator previous = prev(polygon.end());
    for(vector<point>::const_iterator current = polygon.begin(); current != polygon.end(); ++current) {
      if(previous->y <= pkt.y) {
	if(current->y > pkt.y) {
	  if(isLeft(*previous, *current, pkt) > 0) {++nw;}
	}
      } else {
	if(current->y <= pkt.y) {
	  if(isLeft(*previous, *current, pkt) < 0) {--nw;}
	}
      }
      previous = current;
    }
    return(abs(nw) % 2);
  }
  
  //crossing number test for a point in a polygon
  //  0 (false) = outside, 1 (true) = inside
  bool pointInPolyCN(const point &pkt, const vector<point> &polygon) {
    int nx = 0; //Number of crossings
    vector<point>::const_iterator previous = prev(polygon.end());
    for(vector<point>::const_iterator current = polygon.begin(); current != polygon.end(); ++current) {
      if(((previous->y <= pkt.y) && (current->y > pkt.y))
	 || ((previous->y > pkt.y) && (current->y <= pkt.y))) {
	double xc = (pkt.y - previous->y) / (current->y - previous->y);
	if(pkt.x < previous->x + xc * (current->x - previous->x)) {++nx;}
      }
      previous = current;
    }
    return(nx % 2);
  }

  //Points where line cuts a circle
  vector<point> cutLineCircle(const point &a0, const point &a1,
			      const point &c, const double &r) {
    vector<point> ret;
    double dx = a1.x - a0.x;
    double dy = a1.y - a0.y;
    double dr = sqrt(pow(dx, 2) + pow(dy,2));
    double D = (a0.x - c.x) * (a1.y - c.y) - (a1.x - c.x) * (a0.y - c.y);
    double incidence = pow(r,2) * pow(dr,2) - pow(D,2);
    if(incidence > 0) { //==0 .. tangent is not needed here
      array<double, 2> x{};
      array<double, 2> y{};
      double tt = sqrt(incidence);
      x[0] = c.x + (D*dy + copysign(1.0, dy) * dx * tt) / pow(dr,2);
      x[1] = c.x + (D*dy - copysign(1.0, dy) * dx * tt) / pow(dr,2);
      y[0] = c.y + (-D*dx + abs(dy) * tt) / pow(dr,2);
      y[1] = c.y + (-D*dx - abs(dy) * tt) / pow(dr,2);
      for(int i=0; i<2; ++i) {
	if(x[i] >= min(a0.x, a1.x) && x[i] <= max(a0.x, a1.x) &&
	   y[i] >= min(a0.y, a1.y) && y[i] <= max(a0.y, a1.y)) {
	  ret.push_back(point{x[i], y[i]});}
      }
    }
    return(ret);
  }
  
  //Share of circle circumference inside polygon
  double wgtCircC(const point &pkt, const double &radius, const vector<point> &polygon) {
    double wgt = 1.;
    //polygon is totaly inside the circle
    unsigned int nPointsOutside = 0;
    double r2 = pow(radius,2);
    for(auto&& i : polygon) {
      double dist2 = pow(i.x - pkt.x, 2) + pow(i.y - pkt.y, 2);
      if(dist2 > r2) {++nPointsOutside;}
    }
    if(nPointsOutside > 0) {
      vector<point> cutPoint;
      {
	vector<point>::const_iterator previous = prev(polygon.end());
	for(vector<point>::const_iterator current = polygon.begin(); current != polygon.end(); ++current) {
	  vector<point> tmp = cutLineCircle(*previous, *current, pkt, radius);
	  cutPoint.insert(cutPoint.end(), tmp.begin(), tmp.end());
	  previous = current;
	}
      }
      if(cutPoint.size() > 1) {
	vector<double> rad;
	for(auto&& i : cutPoint) {rad.push_back(atan2(pkt.y-i.y, pkt.x-i.x));}
	sort(rad.begin(), rad.end());
	auto last = unique(rad.begin(), rad.end());
	rad.erase(last, rad.end());
	double sum = 0.;
	double sumIn = 0.;
	{
	  vector<double>::const_iterator previous = prev(rad.end());
	  double between = *previous + (*rad.begin() - *previous)/2.;
	  bool inPlot = EiBauBeDi::pointInPolyCN(EiBauBeDi::point{pkt.x + radius * cos(between), pkt.y + radius * sin(between)}, polygon);
	  for(vector<double>::const_iterator current = rad.begin(); current != rad.end(); ++current) {
	    double tmp = fmod(*current - *previous + M_PI * 2., M_PI * 2.);
	    sum += tmp;
	    if(inPlot) {sumIn += tmp;}
	    inPlot = !inPlot;
	    previous = current;
	  }
	}
	if(sum > 0.) {wgt = sumIn/sum;}  //sum should be M_PI * 2.
      }
    } else {wgt = 0.;}
    return(wgt);
  }
  
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
  vector<EiBauBeDi::point> plotCorners;
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
    for(auto&& i : tmpPlotCorners) {plotCorners.push_back(EiBauBeDi::point{i[1],i[2]});}
  }
  //for(auto&& i : plotCorners) cout << i.x << " " << i.y << endl;

  //Read in the tree data
  vector<EiBauBeDi::tree> trees;
  {
    ifstream infile(FileInTrees);
    string line;
    while (getline(infile, line)) {
      if(line[0] == '#') continue;
      istringstream iss(line);
      vector<string> tokens{istream_iterator<string>{iss},
	  istream_iterator<string>{}};
      if(EiBauBeDi::pointInPolyCN(EiBauBeDi::point{stod(tokens[2]),stod(tokens[3])}, plotCorners)) { //Use only trees which are inside the plot
	EiBauBeDi::tree tree;
	tree.nr = tokens[1];
	tree.x = stod(tokens[2]);
	tree.y = stod(tokens[3]);
	tree.z = stod(tokens[4]);
	tree.d = stod(tokens[6]);
	tree.h = stod(tokens[7]);
	tree.hcr = stod(tokens[8]);
	trees.push_back(tree);
      }
    }
  }
  //for(auto&& i : trees) {cout << i.nr << " " << i.x << endl;}

  //Basal area of the plot
  {
    double sumg = 0.;
    for(auto&& i : trees) {sumg += pow(i.d,2);}
    sumg *= M_PI / 4.;
    cout << "Basal area [m2/ha]: " << sumg / polygonArea(plotCorners) << endl;
  }

  //Fix sample plot at the position of each tree
  {
    cout << "\nBasal area from fixed sample plot around each tree" << endl;
    cout << "tree g/ha" << endl;
    double r = 7.; //Sample Radius
    valarray<double> sumg(0., trees.size());
    double r2 = pow(r, 2);
    size_t line = 0;
    for(auto&& i : trees) {
      for(auto&& j : trees) {
	double dist2 = pow(i.x - j.x, 2) + pow(i.y - j.y, 2);
	if(dist2 <= r2) {
	  //Plot border correction
	  double weight = 1./EiBauBeDi::wgtCircC(EiBauBeDi::point{i.x,i.y}, sqrt(dist2), plotCorners);
	  sumg[line] += pow(j.d,2) * weight;
	}
      }
      ++line;
    }
    sumg /= 4.*r2;
    for(size_t i=0; i < trees.size(); ++i) {
      cout << trees[i].nr << " " << sumg[i] << endl;
    }
  }

  //Angle count at the position of each tree
  {
    cout << "\nBasal area from angle count sample around each tree" << endl;
    cout << "tree g/ha" << endl;
    double k = 4.; //counting factor
    valarray<double> sumg(0., trees.size());
    size_t line = 0;
    for(auto&& i : trees) {
      for(auto&& j : trees) {
	double dist2 = pow(i.x - j.x, 2) + pow(i.y - j.y, 2);
	if(dist2 <= 2500 * pow(j.d/100., 2) / k) {
	  //Plot border correction
	  double weight = 1./EiBauBeDi::wgtCircC(EiBauBeDi::point{i.x,i.y}, sqrt(dist2), plotCorners);
	  sumg[line] += k * weight;
	}
      }
      ++line;
    }
    for(size_t i=0; i < trees.size(); ++i) {
      cout << trees[i].nr << " " << sumg[i] << endl;
    }
  }

  //Fix sample plot at systematic rasterpoints
  {
    cout << "\nBasal area from fixed sample plot on raster" << endl;
    cout << "tree g/ha" << endl;
    //Find plot extends
    EiBauBeDi::point plu = {plotCorners[0].x,plotCorners[0].y};
    EiBauBeDi::point pro = {plotCorners[0].x,plotCorners[0].y};
    for(auto&& i : plotCorners) {
      plu.x = min(plu.x, i.x); plu.y = min(plu.y, i.y);
      pro.x = max(pro.x, i.x); pro.y = max(pro.y, i.y);
    }
    double r = 7.; //Sample Radius
    double r2 = pow(r, 2);
    double dxy = 0.1; //Distance between raster points
    valarray<double> sumg(0., trees.size());
    valarray<unsigned int> nSamples(0u, trees.size());
    for(double x = plu.x - r; x <= pro.x + r; x += dxy) {
      for(double y = plu.y - r; y <= pro.y + r; y += dxy) {
	bool insidePlot =
	  EiBauBeDi::pointInPolyCN(EiBauBeDi::point{x, y}, plotCorners);
	double sg = 0.;  //Basal area of this sample
	stack<size_t> treesInSample;
	size_t line = 0;
	for(auto&& i : trees) {
	  if(pow(i.x - x, 2) + pow(i.y - y, 2) <= r2) {
	    sg += pow(i.d,2);
	    treesInSample.push(line);
	    if(insidePlot) {++nSamples[line];}
	  }
	  ++line;
	}
	while(!treesInSample.empty()) {
	  sumg[treesInSample.top()] += sg;
	  treesInSample.pop();}
      }
    }
    sumg /= 4.*r2;
    //sumg /= nSamples;  //Different types
    for(size_t i=0; i<sumg.size(); ++i) {sumg[i] /= nSamples[i];}
    for(size_t i=0; i < trees.size(); ++i) {
      cout << trees[i].nr << " " << sumg[i] << endl;}
  }

  //Fix sample plot at systematic rasterpoints
  //Weightened with the distance between tree and sample center
  {
    cout << "\nBasal area from weightened fixed sample plot on raster" << endl;
    cout << "tree g/ha" << endl;
    //Find plot extends
    EiBauBeDi::point plu = {plotCorners[0].x,plotCorners[0].y};
    EiBauBeDi::point pro = {plotCorners[0].x,plotCorners[0].y};
    for(auto&& i : plotCorners) {
      plu.x = min(plu.x, i.x); plu.y = min(plu.y, i.y);
      pro.x = max(pro.x, i.x); pro.y = max(pro.y, i.y);
    }
    double r = 7.; //Sample Radius
    double r2 = pow(r, 2);
    double dxy = 0.1; //Distance between raster points
    valarray<double> sumg(0., trees.size());
    valarray<double> sumWeight(0., trees.size());
    for(double x = plu.x - r; x <= pro.x + r; x += dxy) {
      for(double y = plu.y - r; y <= pro.y + r; y += dxy) {
	bool insidePlot =
	  EiBauBeDi::pointInPolyCN(EiBauBeDi::point{x, y}, plotCorners);
	double sg = 0.;  //Basal area of this sample
	stack<pair<size_t, double> > treesInSample; //idx, distance
	size_t line = 0;
	for(auto&& i : trees) {
	  double dist2 = pow(i.x - x, 2) + pow(i.y - y, 2);
	  if(dist2 <= r2) {
	    double wgt = 1. - sqrt(dist2) / r;
	    sg += pow(i.d,2);
	    treesInSample.push(make_pair(line, wgt));
	    if(insidePlot) {sumWeight[line] += wgt;}
	  }
	  ++line;
	}
	while(!treesInSample.empty()) {
	  sumg[treesInSample.top().first] += sg * treesInSample.top().second;
	  treesInSample.pop();}
      }
    }
    sumg /= 4.*r2*sumWeight;
    for(size_t i=0; i < trees.size(); ++i) {
      cout << trees[i].nr << " " << sumg[i] << endl;}
  }

  //assign plot area to trees - Winner takes it all
  {
    cout << "\nBasal area from single tree growing area on raster" << endl;
    cout << "tree g/ha distanceCenter directionCenter uncircularity" << endl;
    //Find plot extends
    EiBauBeDi::point plu = {plotCorners[0].x,plotCorners[0].y};
    EiBauBeDi::point pro = {plotCorners[0].x,plotCorners[0].y};
    for(auto&& i : plotCorners) {
      plu.x = min(plu.x, i.x); plu.y = min(plu.y, i.y);
      pro.x = max(pro.x, i.x); pro.y = max(pro.y, i.y);
    }
    double dxy = 0.1; //Distance between raster points
    valarray<unsigned int> sump(0u, trees.size());
    valarray<double> sumpX(0., trees.size());
    valarray<double> sumpY(0., trees.size());
    valarray<double> sumd(0., trees.size());
    for(double x = plu.x + dxy/2.; x <= pro.x; x += dxy) {
      for(double y = plu.y + dxy/2.; y <= pro.y; y += dxy) {
	if(EiBauBeDi::pointInPolyCN(EiBauBeDi::point{x, y}, plotCorners)) {
	  size_t winner = 0;
	  double minInfl = std::numeric_limits<double>::infinity();
	  size_t line = 0;
	  for(auto&& i : trees) {
	    double infl = (pow(i.x - x, 2) + pow(i.y - y, 2)) / pow(i.d,2);
	    if(minInfl > infl) {
	      minInfl = infl;
	      winner = line;
	    }
	    ++line;
	  }
	  ++sump[winner];
	  sumpX[winner] += x; sumpY[winner] += y;
	}
      }
    }
    for(double x = plu.x + dxy/2.; x <= pro.x; x += dxy) {
      for(double y = plu.y + dxy/2.; y <= pro.y; y += dxy) {
	if(EiBauBeDi::pointInPolyCN(EiBauBeDi::point{x, y}, plotCorners)) {
	  size_t winner = 0;
	  double minInfl = std::numeric_limits<double>::infinity();
	  size_t line = 0;
	  for(auto&& i : trees) {
	    double infl = (pow(i.x - x, 2) + pow(i.y - y, 2)) / pow(i.d,2);
	    if(minInfl > infl) {
	      minInfl = infl;
	      winner = line;
	    }
	    ++line;
	  }
	  double centerX = sumpX[winner] / sump[winner];
	  double centerY = sumpY[winner] / sump[winner];
	  sumd[winner] += sqrt(pow(centerX - x, 2) + pow(centerY - y, 2));
	}
      }
    }
    for(size_t i=0; i < trees.size(); ++i) {
      double area = sump[i]*dxy*dxy;
      cout << trees[i].nr << " " << pow(trees[i].d/2.,2)*M_PI/area
	   << " " << sqrt(pow(sumpX[i] / sump[i] - trees[i].x,2) + pow(sumpY[i] / sump[i] - trees[i].y,2))
	   << " " << atan2(sumpY[i] / sump[i] - trees[i].y, sumpX[i] / sump[i] - trees[i].x)
	   << " " << (sumd[i] / sump[i]) / (2./3. * sqrt(area) / M_PI)
	   << endl;
    }
  }

    //assign plot area to trees - Share pixel between trees
  {
    cout << "\nBasal area from single tree growing shared area on raster" << endl;
    cout << "tree g/ha" << endl;
    //Find plot extends
    EiBauBeDi::point plu = {plotCorners[0].x,plotCorners[0].y};
    EiBauBeDi::point pro = {plotCorners[0].x,plotCorners[0].y};
    for(auto&& i : plotCorners) {
      plu.x = min(plu.x, i.x); plu.y = min(plu.y, i.y);
      pro.x = max(pro.x, i.x); pro.y = max(pro.y, i.y);
    }
    double dxy = 0.1; //Distance between raster points
    double maxInfl = 0.04; //trees with higher infl are not used
    valarray<double> sump(0., trees.size());
    for(double x = plu.x + dxy/2.; x <= pro.x; x += dxy) {
      for(double y = plu.y + dxy/2.; y <= pro.y; y += dxy) {
	if(EiBauBeDi::pointInPolyCN(EiBauBeDi::point{x, y}, plotCorners)) {
	  double sInfl = 0.; //summ of influence on this grid
	  stack<pair<size_t, double> > treesInSample; //idx, infl
	  size_t line = 0;
	  for(auto&& i : trees) {
	    double infl = (pow(i.x - x, 2) + pow(i.y - y, 2)) / pow(i.d,2);
	    if(maxInfl > infl) {
	      sInfl += maxInfl - infl;
	      treesInSample.push(make_pair(line, maxInfl - infl));}
	    ++line;
	  }
	  while(!treesInSample.empty()) {
   sump[treesInSample.top().first] += treesInSample.top().second / sInfl;
	    treesInSample.pop();}
 	}
      }
    }
    sump *= dxy*dxy;
    for(size_t i=0; i < trees.size(); ++i) {
      cout << trees[i].nr << " " << pow(trees[i].d/2.,2)*M_PI/sump[i]
	   << endl;
    }
  }

  return(0);
}
