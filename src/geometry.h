#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

struct Point;
class Vector;
class Shape;
class Ellipse;
class Circle;
class Polygon;
class Rectangle;
class Square;
class Triangle;

const double DELTA = 1e-6;
bool isClose(double a, double b);
double determinant(double a11, double a12, double a21, double a22);

struct Point {
    double x, y;
    Point& operator+=(const Vector& a);
    Point& operator-=(const Vector& a);
    Point() : x(0), y(0){};
    Point(double x, double y) : x(x), y(y) {}
    void rotate(const Point& center, double phi);
};

double distance(const Point& a, const Point& b);
bool operator==(const Point& a, const Point& b);
bool operator!=(const Point& a, const Point& b);
bool isClose(const Point& a, const Point& b);
Point midPoint(const Point& a, const Point& b);

class Vector {
  public:
    Vector();
    explicit Vector(const Point& a);
    Vector(double x, double y);
    Vector(const Point& begin, const Point& end);

    Vector& operator+=(const Vector& a);
    Vector& operator-=(const Vector& a);
    Vector& operator*=(double lambda);
    Vector& operator/=(double lambda);
    bool operator==(const Vector& a) const;
    bool operator!=(const Vector& a) const;
    Vector rotate(double phi);
    double norm() const;
    Vector orthogonal() const;
    Vector normalize();

  private:
    double x, y;

    friend double operator*(const Vector& a, const Vector& b);
    friend struct Point;
};

Vector operator+(Vector a, const Vector& b);
Vector operator-(Vector a, const Vector& b);
Vector operator*(Vector a, double lambda);
Vector operator/(Vector a, double lambda);
double operator*(const Vector& a, const Vector& b);
Point operator+(Point a, const Vector& v);
Point operator-(Point a, const Vector& v);
bool isPositivelyOriented(const Vector& a, const Vector& b);
bool orientationSign(const Point& a, const Point& b, const Point& c);
double cosineBetweenVectors(const Vector& a, const Vector& b);

class Line {
  public:
    Line();
    Line(const Point& a, const Point& b);
    Line(const Point& a, double phi);
    Line(const Point& a, const Vector& v);
    Line(double phi, double move);

    double signedDistance(const Point& a) const;
    double distance(const Point& a) const;
    bool onLine(const Point& a) const;
    bool operator==(const Line& another) const;
    bool operator!=(const Line& another) const;
    Vector directingVector() const;

  private:
    double A, B, C;
    inline static const double radius = 100;
    inline static const double zero = 0;
    friend Point commonPoint(const Line& a, const Line& b);
    Line normalize();
};

Line midPerpendicular(const Point& a, const Point& b);
Point commonPoint(const Line& a, const Line& b);
Line midLine(const Point& a, const Point& b, const Point& c);
Line pointPerpendicular(const Point& a, const Line& line);

class Shape {
  public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;

    bool operator==(const Shape& another) const;
    bool operator!=(const Shape& another) const;
    bool isCongruentTo(const Shape& another) const;
    bool isSimilarTo(const Shape& another) const;
    virtual bool containsPoint(const Point& point) const = 0;
    virtual void rotate(const Point& center, double phi) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
    virtual ~Shape() = default;
};

class Polygon : public Shape {
  public:
    Polygon() {}
    explicit Polygon(std::vector<Point>& v) : pointsv(v) {}

    template <typename... Args>
    explicit Polygon(Args... args) {
        (pointsv.push_back(args), ...);
    }

    size_t verticesCount() const;
    const std::vector<Point>& getVertices() const;
    bool isConvex() const;
    Point massCenter() const;
    bool operator==(const Polygon& another) const;
    bool operator!=(const Polygon& another) const;
    bool isCongruentTo(const Polygon& another) const;
    bool isSimilarAsPolygon(const Polygon& another) const;
    double perimeter() const override;
    double area() const override;
    bool containsPoint(const Point& point) const override;
    void rotate(const Point& center, double phi) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;

  protected:
    std::vector<Point> pointsv;
    double subArea(size_t i, size_t j) const;

  private:
    static bool arePointsPermuted(const std::vector<Point>& v1,
                                  const std::vector<Point>& v2, size_t step);

    bool arePolygonsSimilarWithCoefficient(const Polygon& another,
                                           double coefficient) const;

    bool arePolygonsPermuted(const std::vector<Point>& v, size_t start_point,
                             double coefficient) const;
    bool arePolygonsDirectedOrReversed(const std::vector<Point>& v,
                                       size_t start_point,
                                       int64_t multiplicationConst,
                                       double coefficient) const;
};

class Ellipse : public Shape {
  public:
    Ellipse();
    Ellipse(const Point& F1, const Point& F2, double radius);

    Point center() const;
    std::pair<double, double> halfAxes() const;
    double eccentricity() const;
    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    bool operator==(const Ellipse& another) const;
    bool operator!=(const Ellipse& another) const;
    bool isCongruentTo(const Ellipse& another) const;
    bool isSimilarAsEllipse(const Ellipse& another) const;
    double perimeter() const override;
    double area() const override;
    bool containsPoint(const Point& point) const override;
    void rotate(const Point& center, double phi) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;

  protected:
    Point F1, F2;
    double radius;
};

class Circle : public Ellipse {
  public:
    Circle(const Point& center, double radius);
    Circle(const Point& a, const Point& b, const Point& c);
    double radius() const;
};

class Rectangle : public Polygon {
  public:
    Rectangle(const Point& a, const Point& b, double proportion);
    Point center() const;
    std::pair<Line, Line> diagonals() const;
    void Print() const;

  private:
    static std::vector<Point> rectangleVerticesByProportion(const Point& a,
                                                            const Point& b,
                                                            double proportion);
};

class Square : public Rectangle {
  public:
    Square(const Point& a, const Point& b);
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

class Triangle : public Polygon {
  public:
    Triangle(const Point& a, const Point& b, const Point& c);
    Point centroid() const;
    Point orthocenter() const;
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Circle ninePointsCircle() const;
    Line EulerLine() const;
};

#endif  //GEOMETRY_GEOMETRY_H