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

#include "eiBauBeDi.h"

#include <iostream>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

namespace EiBauBeDi {

  polygon::polygon(const std::vector<std::pair<double, double> > &cornerPoints){
    updateCornerPoints(cornerPoints);
  }
  
  size_t polygon::updateCornerPoints(const std::vector<std::pair<double, double> > &cornerPoints) {
    corners.clear();
    corners.reserve(cornerPoints.size());
    for(auto&& i : cornerPoints) {
      corners.push_back(point{i.first, i.second});
    }
    return(corners.size());
  }

  double polygon::getArea() {
    double area = 0.;
    std::vector<point>::const_iterator previous = std::prev(corners.end());
    for(std::vector<point>::const_iterator current = corners.begin(); current != corners.end(); ++current) {
      area += previous->x * current->y;
      area -= previous->y * current->x;
      previous = current;
    }
    area /= 2.;
    return(std::abs(area));
  }

  int polygon::isLeft(const point &p0, const point &p1, const point &p2) {
    return ((p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y));}

  bool polygon::pointInPolyWN(const double &x, const double &y) {
    int nw = 0; //Number of windings
    std::vector<point>::const_iterator previous = prev(corners.end());
    for(std::vector<point>::const_iterator current = corners.begin(); current != corners.end(); ++current) {
      if(previous->y <= y) {
	if(current->y > y) {
	  if(isLeft(*previous, *current, point{x, y}) > 0) {++nw;}
	}
      } else {
	if(current->y <= y) {
	  if(isLeft(*previous, *current, point{x, y}) < 0) {--nw;}
	}
      }
      previous = current;
    }
    return(abs(nw) % 2);
  }

  bool polygon::pointInPolyCN(const double &x, const double &y) {
    int nx = 0; //Number of crossings
    std::vector<point>::const_iterator previous = std::prev(corners.end());
    for(std::vector<point>::const_iterator current = corners.begin(); current != corners.end(); ++current) {
      if(((previous->y <= y) && (current->y > y))
	 || ((previous->y > y) && (current->y <= y))) {
	double xc = (y - previous->y) / (current->y - previous->y);
	if(x < previous->x + xc * (current->x - previous->x)) {++nx;}
      }
      previous = current;
    }
    return(nx % 2);
  }
  
  std::vector<EiBauBeDi::polygon::point> polygon::cutLineCircle(const point &a0, const point &a1, const point &c, const double &r) {
    std::vector<point> ret;
    double dx = a1.x - a0.x;
    double dy = a1.y - a0.y;
    double dr = std::sqrt(pow(dx, 2) + pow(dy,2));
    double D = (a0.x - c.x) * (a1.y - c.y) - (a1.x - c.x) * (a0.y - c.y);
    double incidence = std::pow(r,2) * std::pow(dr,2) - std::pow(D,2);
    if(incidence > 0) { //==0 .. tangent is not needed here
      std::array<double, 2> x{};
      std::array<double, 2> y{};
      double tt = std::sqrt(incidence);
      x[0] = c.x + (D*dy + std::copysign(1.0, dy) * dx * tt) / std::pow(dr,2);
      x[1] = c.x + (D*dy - std::copysign(1.0, dy) * dx * tt) / std::pow(dr,2);
      y[0] = c.y + (-D*dx + std::abs(dy) * tt) / std::pow(dr,2);
      y[1] = c.y + (-D*dx - std::abs(dy) * tt) / std::pow(dr,2);
      for(int i=0; i<2; ++i) {
	if(x[i] >= std::min(a0.x, a1.x) && x[i] <= std::max(a0.x, a1.x) &&
	   y[i] >= std::min(a0.y, a1.y) && y[i] <= std::max(a0.y, a1.y)) {
	  ret.push_back(point{x[i], y[i]});}
      }
    }
    return(ret);
  }

    //Share of circle circumference inside polygon
  double polygon::shareInside(const double &x, const double &y, const double &radius) {
    double wgt = 1.;
    //polygon is totaly inside the circle
    unsigned int nPointsOutside = 0;
    double r2 = std::pow(radius,2);
    for(auto&& i : corners) {
      double dist2 = std::pow(i.x - x, 2) + std::pow(i.y - y, 2);
      if(dist2 > r2) {++nPointsOutside;}
    }
    if(nPointsOutside > 0) {
      std::vector<point> cutPoint;
      {
	std::vector<point>::const_iterator previous = std::prev(corners.end());
	for(std::vector<point>::const_iterator current = corners.begin(); current != corners.end(); ++current) {
	  std::vector<point> tmp = cutLineCircle(*previous, *current, point{x, y}, radius);
	  cutPoint.insert(cutPoint.end(), tmp.begin(), tmp.end());
	  previous = current;
	}
      }
      if(cutPoint.size() > 1) {
	std::vector<double> rad;
	for(auto&& i : cutPoint) {rad.push_back(atan2(y-i.y, x-i.x));}
	std::sort(rad.begin(), rad.end());
	auto last = std::unique(rad.begin(), rad.end());
	rad.erase(last, rad.end());
	double sumIn = 0.;
	std::vector<double>::const_iterator previous = std::prev(rad.end());
	double between = *previous + (*rad.begin() - *previous)/2.;
	bool inPlot = pointInPolyCN(x + radius * cos(between), y + radius * sin(between));
	for(std::vector<double>::const_iterator current = rad.begin(); current != rad.end(); ++current) {
	  double tmp = fmod(*current - *previous + M_PI * 2., M_PI * 2.);
	  if(inPlot) {sumIn += tmp;}
	  inPlot = !inPlot;
	  previous = current;
	}
	wgt = sumIn/(2. * M_PI);
      }
    } else {wgt = 0.;}
    return(wgt);
  }

  std::array<double, 4> polygon::getExtends() {
    std::array<double, 4> ret = {corners[0].x, corners[0].x, corners[0].y, corners[0].y};
    for(auto&& i : corners) {
      ret[0] = std::min(ret[0], i.x); ret[1] = std::max(ret[1], i.x);
      ret[2] = std::min(ret[2], i.y); ret[3] = std::max(ret[3], i.y);
    }
    return(ret);
  }

  
  tree::tree() : nr(""), x(0.), y(0.), z(0.), d(0.), h(0.), hcr(0.), impact(0.), influence0(0.), influence1(0.) {}
  
  tree::tree(const std::string &anr, const double &ax, const double &ay, const double &az, const double &ad, const double &ah, const double &ahcr, const double &aimpact, const double &ainfluence0, const double &ainfluence1) : nr(anr), x(ax), y(ay), z(az), d(ad), h(ah), hcr(ahcr), impact(aimpact), influence0(ainfluence0), influence1(ainfluence1) {}

  //double tree::getWeight(const double &px, const double &py, double f(const double &dx, const double &dy, const double &influence0, const double &influence1)) {
  //return(f(px-x, py-y, influence0, influence1));
  //}

  double tree::getWeight(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree)) {
    return(f(px, py, *this));
  }


  
  forestStand::forestStand(const std::vector<std::pair<double, double> > &cornerPoints) :
    poly(cornerPoints) {}

  double forestStand::getImpactSum() {
    double sum = 0.;
    for(auto&& i : trees) {sum += i.impact;}
    return(sum);
  }

  double forestStand::subsamplePoint(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree), const bool &makeBorderCorrection) {
    double sum=0.;
    for(auto&& i : trees) {
      double wgt = f(px, py, i);
      if(wgt > 0.) {
	double borderCor = 1.;
	if(makeBorderCorrection) {
	  borderCor /= poly.shareInside(px, py, sqrt(pow(i.x-px,2) + pow(i.y-py,2)));}
	sum += i.impact * wgt * borderCor;
      }
    }
    return(sum);
  }

  double circleCircleIntersectionArea(const double &distance, const double &r0, const double &r1) {
    double ret = 0.;
    if(distance < r0 + r1) { //overlap
      if(distance <= std::abs(r0 - r1)) { //One totaly insede the other
	ret = M_PI * std::pow(std::min(r0, r1), 2);
      } else { //Partly overlap
	double r02 = std::pow(r0, 2);
	double r12 = std::pow(r1, 2);
	double x = (r02 - r12 + std::pow(distance, 2)) / (2. * distance);
	double z = std::pow(x, 2);
	double y = std::sqrt(r02 - z);
	ret = r02 * std::asin(y / r0) + r12 * std::asin(y / r1) - y * (x + std::sqrt(z + r12 - r02));
      }
    }
    return(ret);
  }

  double forestStand::subsampleCircle(const double &px, const double &py, const double &r, const bool &makeBorderCorrection) {
    double sum=0.;
    for(auto&& i : trees) {
      double dist = sqrt(pow(i.x - px,2) + pow(i.y - py,2));
      double aoverlap = circleCircleIntersectionArea(dist, i.influence0, r);
	if(aoverlap > 0.) {
	  double borderCor = 1.;
	  if(makeBorderCorrection) {
	    borderCor /= poly.shareInside(px, py, sqrt(pow(i.x-px,2) + pow(i.y-py,2)));}
	  sum += i.impact * aoverlap * borderCor;
	}
    }
    sum /= pow(r,2) * M_PI;
    return(sum);
  }

  std::valarray<double> forestStand::subsamplePointTree(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree), const bool &makeBorderCorrection) {
    std::valarray<double> ret(0., trees.size());
    double sum=0.;
    for(size_t j=0; j<trees.size(); ++j) {
      auto&& i = trees[j];
      double wgt = f(px, py, i);
      if(wgt > 0.) {
	ret[j] = 1.;
	double borderCor = 1.;
	if(makeBorderCorrection) {
	  borderCor /= poly.shareInside(px, py, sqrt(pow(i.x-px,2) + pow(i.y-py,2)));}
	sum += i.impact * wgt * borderCor;
      }
    }
    ret *= sum;
    return(ret);
  }

  std::valarray<double> forestStand::influencePoint(const double &px, const double &py, double f(const double &px, const double &py, const tree &tree), const bool &makeBorderCorrection) {
    std::valarray<double> ret(0., trees.size());
    for(size_t j=0; j<trees.size(); ++j) {
      auto&& i = trees[j];
      double wgt = f(px, py, i);
      if(wgt > 0.) {
	double borderCor = 1.;
	if(makeBorderCorrection) {
	  borderCor /= poly.shareInside(px, py, sqrt(pow(i.x-px,2) + pow(i.y-py,2)));}
	ret[j] += i.impact * wgt * borderCor;
      }
    }
    return(ret);
  }

  std::valarray<double> forestStand::rasterize(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const double &dx, const double &dy, const double &c1, const double &c2, double f0(const tree &tree, const double &dist2, const double &maxDist2, const double &c1, const double &c2)) {
    std::valarray<double> ret(0., trees.size());
    for(auto&& i : trees) {i.maxDist2 = pow(i.influence0,2);}
    std::stack<tree *> treeSortedX; //Trees sorted in x by there influence 
    {
      std::vector<tree *> tmp;
      for(std::vector<tree>::iterator it = trees.begin() ; it != trees.end(); ++it) {
	tmp.push_back(&*it);
      }
      std::sort(tmp.begin(), tmp.end(), [](tree *a, tree *b) {
	  return (a->x - a->influence0) > (b->x - b->influence0);});
      for(auto i : tmp) {treeSortedX.push(i);}
    }
    std::forward_list<tree *> treeWithInfluenceY;
    for(double x = xmin; x <= xmax; x += dx) {
      { //Remove trees which lost influence
	std::forward_list<tree *>::iterator it0 = treeWithInfluenceY.before_begin();
	std::forward_list<tree *>::iterator it1 = treeWithInfluenceY.begin();
	for(;it1 != treeWithInfluenceY.end();) { 
	  if(((*it1)->x + (*it1)->influence0) < x) {
	    ++it1;
	    treeWithInfluenceY.erase_after(it0);
	  } else {it0 = it1++;}
	}
      }
      { //Insert trees which get influence
	std::forward_list<tree *> tmp;
	while(!treeSortedX.empty() && (treeSortedX.top()->x - treeSortedX.top()->influence0) < x) {
	  tmp.push_front(treeSortedX.top());
	  treeSortedX.pop();
	}
	tmp.sort([](tree *a, tree *b) {return (a->y - a->influence0) < (b->y - b->influence0);});
	treeWithInfluenceY.merge(tmp, [](tree *a, tree *b) {return (a->y - a->influence0) < (b->y - b->influence0);});
      }
      std::forward_list<tree *> treeWithInfluence(treeWithInfluenceY);
      for(double y = ymin; y <= ymax; y += dy) {
	if(poly.pointInPolyCN(x, y)) {
	  double inflSum = 0.;
	  std::forward_list<tree *>::iterator it0 = treeWithInfluence.before_begin();
	  std::forward_list<tree *>::iterator it1 = treeWithInfluence.begin();
	  for(;it1 != treeWithInfluence.end() && ((*it1)->y - (*it1)->influence0) < y;) {
	    
	    tree &ctree = **it1;
	    double distance2 = pow(ctree.x - x, 2) + pow(ctree.y - y, 2);
	    if(distance2 < ctree.maxDist2) {
	      ctree.pointInfl = f0(ctree, distance2, ctree.maxDist2, c1, c2);
	      inflSum += ctree.pointInfl;
	    } else if(((*it1)->y + (*it1)->influence0) < y) { //Lost influence
	      it1 = it0;
	      treeWithInfluenceY.erase_after(it0);
	    } else {ctree.pointInfl = 0.;}
	    it0 = it1++;
	  }
	  if(inflSum > 0) {
	    for(it0 = treeWithInfluence.begin(); it0 != it1; ++it0) {
	      ret[*it0 - &*trees.begin()] += (*it0)->pointInfl / inflSum;
	    }
	  }
	}
      }
    }
    ret *= dx*dy;
    return(ret);
  }


    std::valarray<double> forestStand::rasterizeWta(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const double &dx, const double &dy, const double &c1, const double &c2, double f0(const tree &tree, const double &dist2, const double &maxDist2, const double &c1, const double &c2)) {
    std::valarray<double> ret(0., trees.size());
    for(auto&& i : trees) {i.maxDist2 = pow(i.influence0,2);}
    std::stack<tree *> treeSortedX; //Trees sorted in x by there influence 
    {
      std::vector<tree *> tmp;
      for(std::vector<tree>::iterator it = trees.begin() ; it != trees.end(); ++it) {
	tmp.push_back(&*it);
      }
      std::sort(tmp.begin(), tmp.end(), [](tree *a, tree *b) {
	  return (a->x - a->influence0) > (b->x - b->influence0);});
      for(auto i : tmp) {treeSortedX.push(i);}
    }
    std::forward_list<tree *> treeWithInfluenceY;
    for(double x = xmin; x <= xmax; x += dx) {
      { //Remove trees which lost influence
	std::forward_list<tree *>::iterator it0 = treeWithInfluenceY.before_begin();
	std::forward_list<tree *>::iterator it1 = treeWithInfluenceY.begin();
	for(;it1 != treeWithInfluenceY.end();) { 
	  if(((*it1)->x + (*it1)->influence0) < x) {
	    ++it1;
	    treeWithInfluenceY.erase_after(it0);
	  } else {it0 = it1++;}
	}
      }
      { //Insert trees which get influence
	std::forward_list<tree *> tmp;
	while(!treeSortedX.empty() && (treeSortedX.top()->x - treeSortedX.top()->influence0) < x) {
	  tmp.push_front(treeSortedX.top());
	  treeSortedX.pop();
	}
	tmp.sort([](tree *a, tree *b) {return (a->y - a->influence0) < (b->y - b->influence0);});
	treeWithInfluenceY.merge(tmp, [](tree *a, tree *b) {return (a->y - a->influence0) < (b->y - b->influence0);});
      }
      std::forward_list<tree *> treeWithInfluence(treeWithInfluenceY);
      for(double y = ymin; y <= ymax; y += dy) {
	if(poly.pointInPolyCN(x, y)) {
	  std::forward_list<tree *>::iterator it0 = treeWithInfluence.before_begin();
	  std::forward_list<tree *>::iterator it1 = treeWithInfluence.begin();
	  std::forward_list<tree *>::iterator winner = treeWithInfluence.end();
	  double winnersInfluece = INFINITY;
	  for(;it1 != treeWithInfluence.end() && ((*it1)->y - (*it1)->influence0) < y;) {
	    tree &ctree = **it1;
	    double distance2 = pow(ctree.x - x, 2) + pow(ctree.y - y, 2);
	    if(distance2 < ctree.maxDist2) {
	      double influence = f0(ctree, distance2, ctree.maxDist2, c1, c2);
	      if(winnersInfluece > influence) {
		winnersInfluece = influence;
		winner = it1;
	      }
	    } else if(((*it1)->y + (*it1)->influence0) < y) { //Lost influence
	      it1 = it0;
	      treeWithInfluenceY.erase_after(it0);
	    }
	    it0 = it1++;
	  }
	  if(winner != treeWithInfluence.end()) {
	    ++ret[*winner - &*trees.begin()];
	  }
	}
      }
    }
    ret *= dx*dy;
    return(ret);
  }

  
}

