#include "euclidean.h"

#include <stdexcept>

namespace TCERC {

Point::Point(int D)
{
	this->D = D;
	this->coordinate = std::vector<iRRAM::REAL>(D, 0);
}


Point::Point(const std::vector<iRRAM::REAL> &coordinate)
{
	this->D = coordinate.size();
	this->coordinate = coordinate;
}

Point::~Point()
{

}

int Point::ambient_dimension() const
{
	return this->D;
}

Point Point::operator+(const Vector &v) const
{
	if (this->D != v.ambient_dimension())
		throw std::invalid_argument("std::vector addition of different size");

	std::vector<iRRAM::REAL> result;
	for (int i = 0; i < this->D; i++) {
		result.push_back(this->coordinate[i] + v[i]);
	}
	return Point(result);
}

const iRRAM::REAL &Point::operator[](const int i) const
{
	return this->coordinate[i];
}

Point &Point::operator+=(const Vector &v)
{
	if (this->D != v.ambient_dimension())
		throw std::invalid_argument("std::vector addition of different size");

	for (int i = 0; i < this->D; i++) {
		this->coordinate[i] += v[i];
	}
	return *this;
}

Vector Point::operator-(const Point& x) const
{
	if (this->D != x.D)
		throw std::invalid_argument("point subtraction of different size");

	std::vector<iRRAM::REAL> result;
	for (int i = 0; i < this->D; i++) {
		result.push_back(this->coordinate[i] - x.coordinate[i]);
	}
	return Vector(result);
}

Point Point::operator-(const Vector &v) const
{
	if (this->D != v.ambient_dimension())
		throw std::invalid_argument("vector subtraction of different size");

	std::vector<iRRAM::REAL> result;
	for (int i = 0; i < this->D; i++) {
		result.push_back(this->coordinate[i] - v[i]);
	}
	return Point(result);
}

Point &Point::operator-=(const Vector &v)
{
	if (this->D != v.ambient_dimension())
		throw std::invalid_argument("vector subtraction of different size");

	for (int i = 0; i < this->D; i++) {
		this->coordinate[i] -= v[i];
	}
	return *this;
}

Vector::Vector(int D)
{
	this->D = D;
	this->coordinate = std::vector<iRRAM::REAL>(D, 0);
}


Vector::Vector(const std::vector<iRRAM::REAL> &coordinate)
{
	this->D = coordinate.size();
	this->coordinate = coordinate;
}

Vector::~Vector()
{

}

Vector Vector::operator+(const Vector &v) const
{
	if (this->D != v.D)
		throw std::invalid_argument("vector addition of different size");

	std::vector<iRRAM::REAL> result;
	for (int i = 0; i < this->D; i++) {
		result.push_back(this->coordinate[i] + v.coordinate[i]);
	}
	return Vector(result);
}

Vector &Vector::operator+=(const Vector &v)
{
	if (this->D != v.D)
		throw std::invalid_argument("vector addition of different size");

	for (int i = 0; i < this->D; i++) {
		this->coordinate[i] += v.coordinate[i];
	}
	return *this;
}

Vector Vector::operator-(const Vector &v) const
{
	if (this->D != v.D)
		throw std::invalid_argument("vector subtraction of different size");

	std::vector<iRRAM::REAL> result;
	for (int i = 0; i < this->D; i++) {
		result.push_back(this->coordinate[i] - v.coordinate[i]);
	}
	return Vector(result);
}

Vector &Vector::operator-=(const Vector &v)
{
	if (this->D != v.D)
		throw std::invalid_argument("vector subtraction of different size");

	for (int i = 0; i < this->D; i++) {
		this->coordinate[i] -= v.coordinate[i];
	}
	return *this;
}

iRRAM::REAL Vector::operator*(const Vector& v) const
{
	if (this->D != v.D)
		throw std::invalid_argument("dot product of different size");

	iRRAM::REAL result(0);
	for (int i = 0; i < this->D; i++){
		result += this->coordinate[i] * v.coordinate[i];
	}
	return result;
}

Vector Vector::operator*(const iRRAM::REAL &r) const
{
	std::vector<iRRAM::REAL> result;
	for (int i = 0; i < this->D; i++){
		result.push_back(this->coordinate[i] * r);
	}
	return Vector(result);
}

Vector &Vector::operator*=(const iRRAM::REAL &r)
{
	for (int i = 0; i < this->D; i++) {
		this->coordinate[i] *= r;
	}
	return *this;
}

Vector Vector::operator/(const iRRAM::REAL& r) const
{
	std::vector<iRRAM::REAL> result;
	for (int i = 0; i < this->D; i++){
		result.push_back(this->coordinate[i] / r);
	}
	return Vector(result);
}

Vector &Vector::operator/=(const iRRAM::REAL &r)
{
	for (int i = 0; i < this->D; i++) {
		this->coordinate[i] /= r;
	}
	return *this;
}

const iRRAM::REAL &Vector::operator[](const int i) const
{
	return this->coordinate[i];
}

iRRAM::REAL &Vector::operator[](const int i)
{
	return this->coordinate[i];
}

int Vector::ambient_dimension() const
{
	return this->D;
}

std::vector<iRRAM::REAL> Vector::get_coordinate() const
{
	return this->coordinate;
}

}
