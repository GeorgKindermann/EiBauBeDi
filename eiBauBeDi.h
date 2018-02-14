// eiBauBeDi 0.1a, EInzelBAUmBEstandesDIchte - Single Tree Stand Density
// Programm to calculte the stand density of singe trees in a forest stand
// Copyright (C) 2017-2018 Georg Kindermann
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

#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <valarray>
#include <stack>
#include <numeric>
#include <tuple>
#include <forward_list>

namespace EiBauBeDi {

  class polygon {
  public:
    polygon(const std::vector<std::pair<double, double> > &cornerPoints);
    //~polygon();
    size_t updateCornerPoints(const std::vector<std::pair<double, double> > &cornerPoints);
    double getArea();
    //Test if point is inside polygon
    //  0 (false) = outside, 1 (true) = inside
    bool pointInPolyWN(const double &x, const double &y); //winding number
    bool pointInPolyCN(const double &x, const double &y); //crossing number
    //Share of circle circumference inside polygon
    double shareInside(const double &x, const double &y, const double &radius);
    std::array<double, 4> getExtends();
  private:
    struct point {
      double x;
      double y;
    };
    std::vector<point> corners;
    //tests if a point is Left|On|Right of an infinite line.
    //  >0 for P2 left of the line through P0 and P1, 0 for P2  on the line,
    //  <0 for P2  right of the line
    int isLeft(const point &p0, const point &p1, const point &p2);
    //Points where line cuts a circle
    std::vector<EiBauBeDi::polygon::point> cutLineCircle(const point &a0, const point &a1, const point &c, const double &r);
  };

  class tree {
  public:
    tree();
    tree(const std::string &nr, const double &x, const double &y, const double &z, const double &d, const double &h, const double &hcr, const double &impact, const double &influence0, const double &influence1);
    //~tree();
    std::string nr; //tree indentifier
    double x; //x Position
    double y; //y Position
    double z; //z Position
    double d; //Diameter
    double h; //Height
    double hcr; //Height of the crown base
    double impact; //e.g. the basal area, Volumne ...
    double influence0;
    double influence1;
    //double getWeight(const double &px, const double &py, double f(const double &dx, const double &dy, const double &influence0, const double &influence1));
    double getWeight(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree));
    friend class forestStand;
  protected:
    double pointInfl;
    double maxDist2;
  private:
  };

  //Forest stand
  class forestStand {
  public:
    forestStand(const std::vector<std::pair<double, double> > &cornerPoints);
    //~stand();
    polygon poly;
    double getImpactSum();
    double subsamplePoint(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree), const bool &makeBorderCorrection); //Sample at x, y
    double subsampleCircle(const double &px, const double &py, const double &r, const bool &makeBorderCorrection);
    std::valarray<double> subsamplePointTree(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree), const bool &makeBorderCorrection);
    std::valarray<double> influencePoint(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree), const bool &makeBorderCorrection);
    std::vector<tree> trees;
    std::valarray<double> rasterize(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const double &dx, const double &dy, const double &c1, const double &c2, double f0(const tree &tree, const double &dist2, const double &maxDist2, const double &c1, const double &c2));
    std::valarray<double> rasterizeWta(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const double &dx, const double &dy, const double &c1, const double &c2, double f0(const tree &tree, const double &dist2, const double &maxDist2, const double &c1, const double &c2));
  private:
  };

  
}
