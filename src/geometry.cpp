#include "geometry.h"
//
// Point
//

bool isClose(double a, double b) {
    return std::abs(a - b) <= DELTA;
}

double determinant(double a11, double a12, double a21, double a22) {
    return a11 * a22 - a12 * a21;
}

double distance(const Point& a, const Point& b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

bool operator==(const Point& a, const Point& b) {
    return isClose(a.x, b.x) && isClose(a.y, b.y);
}

bool operator!=(const Point& a, const Point& b) {
    return !(a == b);
}

Point& Point::operator+=(const Vector& a) {
    x += a.x;
    y += a.y;
    return *this;
}
Point& Point::operator-=(const Vector& a) {
    x -= a.x;
    y -= a.y;
    return *this;
}

Point midPoint(const Point& a, const Point& b) {
    return {(a.x + b.x) / 2, (a.y + b.y) / 2};
}

bool isClose(const Point& a, const Point& b) {
    return a == b;
}

//
// Vector
//

Vector::Vector() : Vector(0, 0) {}

Vector::Vector(double x, double y) : x(x), y(y) {}

Vector::Vector(const Point& a) : x(a.x), y(a.y){};

Vector::Vector(const Point& begin, const Point& end)
    : x(end.x - begin.x), y(end.y - begin.y){};

Vector& Vector::operator+=(const Vector& a) {
    x += a.x;
    y += a.y;
    return *this;
}

Vector& Vector::operator-=(const Vector& a) {
    x -= a.x;
    y -= a.y;
    return *this;
}

Vector& Vector::operator*=(double lambda) {
    x *= lambda;
    y *= lambda;
    return *this;
}
Vector& Vector::operator/=(double lambda) {
    x /= lambda;
    y /= lambda;
    return *this;
}

bool Vector::operator==(const Vector& a) const {
    return isClose(x, a.x) && isClose(y, a.y);
}

bool Vector::operator!=(const Vector& a) const {
    return !(*this == a);
}

Vector Vector::rotate(double phi) {
    *this = Vector(x * cos(phi) - y * sin(phi), x * sin(phi) + y * cos(phi));
    return *this;
}

Vector operator+(Vector a, const Vector& b) {
    return a += b;
}

Vector operator-(Vector a, const Vector& b) {
    return a -= b;
}

Vector operator*(Vector a, double lambda) {
    return a *= lambda;
}
Vector operator/(Vector a, double lambda) {
    return a /= lambda;
}

double operator*(const Vector& a, const Vector& b) {
    return a.x * b.x + a.y * b.y;
}

double Vector::norm() const {
    return sqrt(pow(x, 2) + pow(y, 2));
}

Vector Vector::orthogonal() const {
    Vector copy = *this;
    return copy.rotate(std::numbers::pi / 2);
}

Vector Vector::normalize() {
    double modulus = norm();
    x /= modulus;
    y /= modulus;
    return *this;
}

bool isPositivelyOriented(const Vector& a, const Vector& b) {
    return a.orthogonal() * b > 0 - DELTA;
}

bool orientationSign(const Point& a, const Point& b, const Point& c) {
    Vector v1(b, a), v2(c, b);
    return isPositivelyOriented(v1, v2);
}

Point operator+(Point a, const Vector& v) {
    return (a += v);
}

Point operator-(Point a, const Vector& v) {
    return (a -= v);
}

double cosineBetweenVectors(const Vector& a, const Vector& b) {
    return a * b / (a.norm() * b.norm());
}

//
// Line
//

bool Line::operator==(const Line& another) const {
    return isClose(A, another.A) && isClose(B, another.B) &&
           isClose(C, another.C);
}

bool Line::operator!=(const Line& another) const {
    return !(*this == another);
}

Line::Line() : Line(Point(0, 0), Point(1, 0)) {}

Line::Line(const Point& a, const Point& b)
    : A(a.y - b.y), B(b.x - a.x), C(a.x * (b.y - a.y) + a.y * (a.x - b.x)) {
    normalize();
}

Line::Line(const Point& a, double phi)
    : Line(a, Point(a.x + radius, a.y + radius * phi)) {}

Line::Line(const Point& a, const Vector& v) : Line(a, a + v) {}

Line::Line(double phi, double move)
    : Line(Point(0, move), Point(radius, move + radius * phi)) {}

Line Line::normalize() {
    double modulus = sqrt(pow(A, 2) + pow(B, 2));
    A /= modulus;
    B /= modulus;
    C /= modulus;
    if (A < 0) {
        A = -A;
        B = -B;
        C = -C;
    }
    return *this;
}

double Line::signedDistance(const Point& a) const {
    return (A * a.x + B * a.y) / sqrt(pow(A, 2) + pow(B, 2));
}

double Line::distance(const Point& a) const {
    return std::abs(signedDistance(a));
}

bool Line::onLine(const Point& a) const {
    return isClose(distance(a), zero);
}

Vector Line::directingVector() const {
    return Vector(B, -A);
}

Point commonPoint(const Line& a, const Line& b) {
    double detMain = determinant(a.A, a.B, b.A, b.B);
    if (detMain == 0) {
        return Point();
    }
    double detX = determinant(-a.C, a.B, -b.C, b.B);
    double detY = determinant(a.A, -a.C, b.A, -b.C);
    return {detX / detMain, detY / detMain};
}

Line midPerpendicular(const Point& a, const Point& b) {
    return {midPoint(a, b), Vector(a, b).orthogonal()};
}

Line midLine(const Point& a, const Point& b, const Point& c) {
    Vector v1(b, a), v2(b, c);
    v1.normalize();
    v2.normalize();
    Vector v3 = v1 + v2;
    return Line(b, v3);
}

Line pointPerpendicular(const Point& a, const Line& line) {
    return Line(a, line.directingVector().orthogonal());
}

void Point::rotate(const Point& center, double phi) {
    Vector v(center, *this);
    v.rotate(phi);
    *this = center + v;
}

//
// Shape
//

bool Shape::operator==(const Shape& another) const {
    std::string typeName1 = typeid(*this).name();
    typeName1 = typeName1.substr(1);
    std::string typeName2 = typeid(another).name();
    typeName2 = typeName2.substr(1);
    bool isEllipse1 = (typeName1 == "Ellipse" || typeName1 == "Circle");
    bool isEllipse2 = (typeName2 == "Ellipse" || typeName2 == "Circle");
    if (isEllipse1 == isEllipse2) {
        return (isEllipse1 ? *dynamic_cast<const Ellipse*>(this) ==
                                 *dynamic_cast<const Ellipse*>(&another)
                           : *dynamic_cast<const Polygon*>(this) ==
                                 *dynamic_cast<const Polygon*>(&another));
    }
    return false;
}

bool Shape::operator!=(const Shape& another) const {
    return !(*this == another);
}

bool Shape::isCongruentTo(const Shape& another) const {
    std::string typeName1 = typeid(*this).name();
    typeName1 = typeName1.substr(1);
    std::string typeName2 = typeid(another).name();
    typeName2 = typeName2.substr(1);
    bool isEllipse1 = (typeName1 == "Ellipse" || typeName1 == "Circle");
    bool isEllipse2 = (typeName2 == "Ellipse" || typeName2 == "Circle");
    if (isEllipse1 == isEllipse2) {
        return (isEllipse1 ? dynamic_cast<const Ellipse*>(this)->isCongruentTo(
                                 *dynamic_cast<const Ellipse*>(&another))
                           : dynamic_cast<const Polygon*>(this)->isCongruentTo(
                                 *dynamic_cast<const Polygon*>(&another)));
    }
    return false;
}

bool Shape::isSimilarTo(const Shape& another) const {
    std::string typeName1 = typeid(*this).name();
    std::string typeName2 = typeid(another).name();
    bool isEllipse1 = (typeName1 == typeid(Ellipse).name() ||
                       typeName1 == typeid(Circle).name());
    bool isEllipse2 = (typeName2 == typeid(Ellipse).name() ||
                       typeName2 == typeid(Circle).name());
    if (isEllipse1 == isEllipse2) {
        return (isEllipse1
                    ? dynamic_cast<const Ellipse*>(this)->isSimilarAsEllipse(
                          *dynamic_cast<const Ellipse*>(&another))
                    : dynamic_cast<const Polygon*>(this)->isSimilarAsPolygon(
                          *dynamic_cast<const Polygon*>(&another)));
    }
    return false;
}

//
//Polygon
//

size_t Polygon::verticesCount() const {
    return pointsv.size();
}

const std::vector<Point>& Polygon::getVertices() const {
    return pointsv;
}

bool Polygon::isConvex() const {
    size_t size = pointsv.size();
    bool orientation_sign =
        orientationSign(pointsv[size - 1], pointsv[0], pointsv[1]);
    for (size_t i = 1; i < size; ++i) {
        bool tmp_sign = orientationSign(pointsv[i - 1], pointsv[i],
                                        pointsv[(i + 1) % size]);
        if (tmp_sign != orientation_sign) {
            return false;
        }
    }
    return true;
}

Point Polygon::massCenter() const {
    size_t size = pointsv.size();
    double x = 0, y = 0;
    for (size_t i = 0; i < size; ++i) {
        x += pointsv[i].x;
        y += pointsv[i].y;
    }
    x /= size;
    y /= size;
    return {x, y};
}

bool Polygon::operator==(const Polygon& another) const {
    if (verticesCount() != another.verticesCount()) {
        return false;
    }
    return arePointsPermuted(pointsv, another.pointsv, 1);
}

bool Polygon::operator!=(const Polygon& another) const {
    return !(*this == another);
}

bool Polygon::isCongruentTo(const Polygon& another) const {
    return arePolygonsSimilarWithCoefficient(another, 1);
}

bool Polygon::isSimilarAsPolygon(const Polygon& another) const {
    return arePolygonsSimilarWithCoefficient(another,
                                             perimeter() / another.perimeter());
}

double Polygon::perimeter() const {
    size_t size = pointsv.size();
    double ans = 0;
    for (size_t i = 0; i < size; ++i) {
        ans += distance(pointsv[i], pointsv[(i + 1) % size]);
    }
    return ans;
}

double Polygon::area() const {
    size_t size = pointsv.size();
    double ans = 0;
    for (size_t i = 0; i < size; ++i) {
        ans += subArea(i, (i + 1) % size);
    }
    return std::abs(ans);
}

double Polygon::subArea(size_t i, size_t j) const {
    return (pointsv[i].x - pointsv[j].x) * (pointsv[i].y + pointsv[j].y) / 2;
}

bool Polygon::containsPoint(const Point& point) const {
    size_t size = pointsv.size();
    bool sign = isPositivelyOriented(Vector(pointsv[0], pointsv[1]),
                                     Vector(pointsv[0], point));
    for (size_t i = 0; i < size; ++i) {
        Vector v(pointsv[i], pointsv[(i + 1) % size]);
        if (Line(pointsv[i], v).onLine(point) &&
            v * Vector(pointsv[i], point) > 0 &&
            distance(pointsv[i], point) <
                distance(pointsv[i], pointsv[(i + 1) % size]) + DELTA) {
            return true;
        }
        if (sign != isPositivelyOriented(v, Vector(pointsv[i], point)) &&
            !Line(pointsv[i], v).onLine(point)) {
            return false;
        }
    }
    return true;
}

void Polygon::rotate(const Point& center, double phi) {
    for (size_t i = 0; i < pointsv.size(); ++i) {
        pointsv[i].rotate(center, phi);
    }
}

void Polygon::reflect(const Point& center) {
    for (size_t i = 0; i < pointsv.size(); ++i) {
        pointsv[i] = center + Vector(pointsv[i], center);
    }
}

void Polygon::reflect(const Line& axis) {
    Point dot;
    size_t i = 0;
    while (i < pointsv.size() && axis.onLine(pointsv[i])) {
        ++i;
    }
    if (i != pointsv.size()) {
        dot = pointsv[i];
        Vector v = axis.directingVector().orthogonal().normalize();
        double dist = axis.signedDistance(dot);
        if (!axis.onLine(dot + v * dist)) {
            v *= -1;
        }
        for (size_t i = 0; i < pointsv.size(); ++i) {
            pointsv[i] += v * 2 * axis.signedDistance(pointsv[i]);
        }
    }
}

void Polygon::scale(const Point& center, double coefficient) {
    for (size_t i = 0; i < pointsv.size(); ++i) {
        Vector v(center, pointsv[i]);
        pointsv[i] = center + v * coefficient;
    }
}

bool Polygon::arePolygonsSimilarWithCoefficient(const Polygon& another,
                                                double coefficient) const {
    if (verticesCount() != another.verticesCount()) {
        return false;
    }
    size_t size = verticesCount();
    double cosanotherBegin = cosineBetweenVectors(
        Vector(another.pointsv[0], another.pointsv[1]),
        Vector(another.pointsv[0], another.pointsv.back()));
    for (size_t i = 0; i < size; ++i) {
        if (isClose(cosanotherBegin,
                    cosineBetweenVectors(
                        Vector(pointsv[i], pointsv[(i + 1) % size]),
                        Vector(pointsv[i], pointsv[(i + size - 1) % size]))) &&
            arePolygonsPermuted(another.pointsv, i, coefficient)) {
            return true;
        }
    }
    return false;
}

bool Polygon::arePointsPermuted(const std::vector<Point>& v1,
                                const std::vector<Point>& v2, size_t step) {
    if (v1.size() != v2.size()) {
        return false;
    }
    size_t start_index = 0, size = v1.size();
    while (start_index < size) {
        while (start_index < size && v1[start_index] != v2[0]) {
            start_index += step;
        }
        bool directPerm = true, reversePerm = true;
        for (size_t i = 0; i < size && directPerm; ++i) {
            directPerm = isClose(v1[(start_index + i) % size], v2[i]);
        }
        for (size_t i = 0; i < size && !directPerm && reversePerm; ++i) {
            reversePerm =
                isClose(v1[(start_index + i) % size], v2[(size - i) % size]);
        }
        if (directPerm || reversePerm) {
            return true;
        }
    }
    return false;
}

bool Polygon::arePolygonsPermuted(const std::vector<Point>& v,
                                  size_t start_point,
                                  double coefficient) const {
    return arePolygonsDirectedOrReversed(v, start_point, 1, coefficient) ||
           arePolygonsDirectedOrReversed(v, start_point, -1, coefficient);
}

bool Polygon::arePolygonsDirectedOrReversed(const std::vector<Point>& v,
                                            size_t start_point,
                                            int64_t multiplicationConst,
                                            double coefficient) const {
    int64_t size = pointsv.size();
    bool isEqual = true;
    for (int64_t i = 0; i < size && isEqual; ++i) {
        Vector v1(pointsv[(start_point + i) % size],
                  pointsv[(start_point + i + 1) % size]);
        Vector v2(pointsv[(start_point + i) % size],
                  pointsv[(start_point + i + size - 1) % size]);
        double cosA = cosineBetweenVectors(v1, v2);
        Vector anotherV1(v[(size + i * multiplicationConst) % size],
                         v[(size + i * multiplicationConst + 1) % size]);
        Vector anotherV2(v[(size + i * multiplicationConst) % size],
                         v[(size + i * multiplicationConst - 1) % size]);
        double cosB = cosineBetweenVectors(anotherV1, anotherV2);
        isEqual =
            isClose(cosA, cosB) &&
            isClose(
                distance(pointsv[(start_point + i) % size],
                         pointsv[(start_point + i + 1) % size]),
                coefficient *
                    distance(v[(size + i * multiplicationConst) % size],
                             v[(size + i * multiplicationConst - 1) % size]));
    }
    return isEqual;
}

//
// Ellipse
//

Ellipse::Ellipse() : F1(Point()), F2(Point()), radius(0) {}

Ellipse::Ellipse(const Point& F1, const Point& F2, double radius)
    : F1(F1), F2(F2), radius(radius) {}

Point Ellipse::center() const {
    return {(F1.x + F2.x) / 2, (F1.y + F2.y) / 2};
}

std::pair<double, double> Ellipse::halfAxes() const {
    Point center = this->center();
    double center_focus_distance = distance(center, F1);
    double b = sqrt(pow(radius, 2) - pow(2 * center_focus_distance, 2)) / 2;
    double a = radius / 2;
    return std::make_pair(a, b);
}

double Ellipse::eccentricity() const {
    std::pair<double, double> halfAxes = this->halfAxes();
    double a = halfAxes.first;
    double b = halfAxes.second;
    return sqrt(1 - pow(b / a, 2));
}

std::pair<Point, Point> Ellipse::focuses() const {
    return std::make_pair(F1, F2);
}

std::pair<Line, Line> Ellipse::directrices() const {
    Point center = this->center(), d1, d2;
    Vector F1F2(F1, F2);
    F1F2.normalize();
    std::pair<double, double> halfAxes = this->halfAxes();
    double eccentricity = this->eccentricity();
    double a = halfAxes.first;
    double delta = a / eccentricity;
    F1F2 *= delta;
    d1 = center + F1F2;
    d2 = center - F1F2;
    F1F2.orthogonal();
    return std::make_pair(Line(d1, F1F2), Line(d2, F1F2));
}

bool Ellipse::operator==(const Ellipse& another) const {
    return ((F1 == another.F1 && F2 == another.F2) ||
            (F1 == another.F2 && F2 == another.F1)) &&
           isClose(radius, another.radius);
}

bool Ellipse::operator!=(const Ellipse& another) const {
    return !(*this == another);
}

bool Ellipse::isCongruentTo(const Ellipse& another) const {
    return isClose(distance(another.F1, another.F2), distance(F1, F2)) &&
           isClose(radius, another.radius);
}

bool Ellipse::isSimilarAsEllipse(const Ellipse& another) const {
    Ellipse copy = another;
    copy.scale(another.center(), area() / another.area());
    return isCongruentTo(copy);
}

double Ellipse::perimeter() const {
    std::pair<double, double> halfAxes = this->halfAxes();
    double& a = halfAxes.first;
    double& b = halfAxes.second;
    return std::numbers::pi * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)));
}

double Ellipse::area() const {
    std::pair<double, double> halfAxes = this->halfAxes();
    return std::numbers::pi * halfAxes.first * halfAxes.second;
}

bool Ellipse::containsPoint(const Point& point) const {
    return distance(point, F1) + distance(point, F2) <= radius + DELTA;
}

void Ellipse::rotate(const Point& center, double phi) {
    F1.rotate(center, phi);
    F2.rotate(center, phi);
}

void Ellipse::reflect(const Point& center) {
    Vector v1(F1, center), v2(F2, center);
    F1 = center + v1;
    F2 = center + v2;
}

void Ellipse::reflect(const Line& axis) {
    Vector v = axis.directingVector().orthogonal().normalize();
    double dist1 = axis.signedDistance(F1);
    if (!axis.onLine(F1 + v * dist1)) {
        v *= -1;
    }
    F1 += v * 2 * axis.signedDistance(F1);
    F2 += v * 2 * axis.signedDistance(F2);
}

void Ellipse::scale(const Point& center, double coefficient) {
    Vector v1(center, F1);
    F1 = center + v1 * coefficient;
    Vector v2(center, F2);
    F2 = center + v2 * coefficient;
    radius *= coefficient;
}

//
// Circle
//

Circle::Circle(const Point& center, double radius)
    : Ellipse(center, center, 2 * radius){};

Circle::Circle(const Point& a, const Point& b, const Point& c) {
    *this = Triangle(a, b, c).circumscribedCircle();
}

double Circle::radius() const {
    return Ellipse::radius / 2;
}

//
// Rectangle
//

std::vector<Point> Rectangle::rectangleVerticesByProportion(const Point& a,
                                                            const Point& b,
                                                            double proportion) {
    proportion = (proportion > 1 ? 1 / proportion : proportion);
    double hypotenuse = distance(a, b);
    double sin = 1 / sqrt(pow(proportion, 2) + 1);
    double cos = sqrt(1 - pow(sin, 2));
    double long_side = hypotenuse * cos;
    double short_side = hypotenuse * sin;
    Vector v(a, b), v1(a, b);
    v.normalize();
    v1.normalize();
    v.rotate(-std::numbers::pi / 2 + acos(cos));
    v1.rotate(acos(cos));
    v *= short_side;
    v1 *= long_side;
    Point c = a + v;
    Point d = a + v1;
    return {a, d, b, c};
}

Rectangle::Rectangle(const Point& a, const Point& b, double proportion) {
    pointsv = rectangleVerticesByProportion(a, b, proportion);
}

Point Rectangle::center() const {
    return massCenter();
}

std::pair<Line, Line> Rectangle::diagonals() const {
    return std::make_pair(Line(pointsv[0], pointsv[2]),
                          Line(pointsv[1], pointsv[3]));
}

//
// Square
//

Square::Square(const Point& a, const Point& b) : Rectangle(a, b, 1){};

Circle Square::circumscribedCircle() const {
    return Circle(center(), distance(pointsv[0], pointsv[2]) / 2);
}

Circle Square::inscribedCircle() const {
    return Circle(center(), distance(pointsv[0], pointsv[1]) / 2);
}

//
// Triangle
//

Triangle::Triangle(const Point& a, const Point& b, const Point& c)
    : Polygon(a, b, c){};

Point Triangle::centroid() const {
    return massCenter();
}

Point Triangle::orthocenter() const {
    const Point &a = pointsv[0], &b = pointsv[1], &c = pointsv[2];
    Line orthogonalLine1 = pointPerpendicular(a, Line(b, c));
    Line orthogonalLine2 = pointPerpendicular(b, Line(a, c));
    return commonPoint(orthogonalLine1, orthogonalLine2);
}

Circle Triangle::circumscribedCircle() const {
    const Point &a = pointsv[0], &b = pointsv[1], &c = pointsv[2];
    Point center = commonPoint(midPerpendicular(a, b), midPerpendicular(a, c));
    return {center, distance(center, a)};
}

Circle Triangle::inscribedCircle() const {
    const Point &a = pointsv[0], &b = pointsv[1], &c = pointsv[2];
    Line midLine1 = midLine(a, b, c), midLine2 = midLine(b, c, a), line(a, b);
    Point center = commonPoint(midLine1, midLine2);
    return Circle(center, line.distance(center));
}

Circle Triangle::ninePointsCircle() const {
    Circle circle = this->circumscribedCircle();
    return {midPoint(orthocenter(), circle.center()), circle.radius() / 2};
}

Line Triangle::EulerLine() const {
    return Line(orthocenter(), centroid());
}