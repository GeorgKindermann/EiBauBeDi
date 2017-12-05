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
    r = pow(r, 2);
    size_t line = 0;
    for(auto&& i : trees) {
      for(auto&& j : trees) {
	if(pow(i.x - j.x, 2) + pow(i.y - j.y, 2) <= r) {
	  //Plot border correction is currently MISSING
	  double weight = 1.;
	  sumg[line] += pow(j.d,2) * weight;
	}
      }
      ++line;
    }
    sumg /= 4.*r;
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
	if(pow(i.x - j.x, 2) + pow(i.y - j.y, 2) <= 2500 * pow(j.d/100., 2) / k) {
	  //Plot border correction is currently MISSING
	  double weight = 1.;
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

  return(0);
}
