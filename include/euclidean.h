/**
 * @file euclidean.h
 * @author Joonhyung Shin (tonyshin@kaist.ac.kr)
 * @author Namjo Ahn
 * @author Seungwoo Lee
 * @date May, 2018
 * @brief This is the header file of Point and Vector.
 */

#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H

#include <iRRAM/lib.h>
#include <vector>

namespace TCERC {

class Vector;
class Point;

/**
 * @brief An affine point in a Euclidean space.
 */
class Point {
private:
	int D;
	std::vector<iRRAM::REAL> coordinate;
public:
	/* Constructors. */

	/**
	 * @brief Construct a point located at the origin.
	 *
	 * @param D the dimension of a space the point lies in
	 */
	Point(int D);
	Point(const std::vector<iRRAM::REAL> &coordinate);

	/* Destructor. */
	~Point();

	/* Basic arithmetics. */
	Point operator+(const Vector &v) const;
	Point &operator+=(const Vector &v);
	Vector operator-(const Point &x) const;
	Point operator-(const Vector &v) const;
	Point &operator-=(const Vector &v);

	/* Comparisons. */
	friend inline iRRAM::LAZY_BOOLEAN operator==(const Point &x, const Point &y);
	friend inline iRRAM::LAZY_BOOLEAN operator!=(const Point &x, const Point &y);

	/* Basic functionalities. */

	/**
	 * @brief The dimension of the ambient space.
	 *
	 * @return the dimension of the ambient space
	 */
	int ambient_dimension() const;

	const iRRAM::REAL &operator[](const int i) const;
	/* Do we really need this for affine object? */
	//std::vector<iRRAM::REAL> get_coordinate() const;
};

/**
 * @brief A vector lying in a Euclidean (vector) space.
 */
class Vector {
private:
	int D;
	std::vector<iRRAM::REAL> coordinate;
public:
	/* Constructors. */

	/**
	 * @brief Construct a zero vector.
	 *
	 * @param the dimension of the vector.
	 */
	Vector(int D);

	Vector(const std::vector<iRRAM::REAL> &coordinate);

	/* Destructor. */
	~Vector();

	/* Basic arithmetics. */
	Vector operator+(const Vector &v) const;
	Vector &operator+=(const Vector &v);
	Vector operator-(const Vector &v) const;
	Vector &operator-=(const Vector &v);
	iRRAM::REAL operator*(const Vector &v) const;	//Dot Product
	Vector operator*(const iRRAM::REAL &r) const;	//Scalar Multiplication
	Vector &operator*=(const iRRAM::REAL &r);
	Vector operator/(const iRRAM::REAL &r) const;	//Scalar Division - Can't check 0
	Vector &operator/=(const iRRAM::REAL &r);

	/* Comparisons. */
	friend inline iRRAM::LAZY_BOOLEAN operator==(const Vector &x, const Vector &y);
	friend inline iRRAM::LAZY_BOOLEAN operator!=(const Vector &x, const Vector &y);

	/* Basic functionalities. */
	const iRRAM::REAL &operator[](const int i) const;
	iRRAM::REAL &operator[](const int i);

	/**
	 * @brief The dimension of the vector.
	 *
	 * @return the dimension of the vector
	 */
	int ambient_dimension() const;

	/**
	 * @brief The standard coordinates of the vector.
	 *
	 * @return the standard coordinates of the vector.
	 */
	std::vector<iRRAM::REAL> get_coordinate() const;
};

inline iRRAM::LAZY_BOOLEAN operator==(const Point &x, const Point &y)
{
	return (x - y) * (x - y) == 0;
}

inline iRRAM::LAZY_BOOLEAN operator!=(const Point &x, const Point &y)
{
	return (x - y) * (x - y) != 0;
}

inline iRRAM::LAZY_BOOLEAN operator==(const Vector &x, const Vector &y)
{
	return (x - y) * (x - y) == 0;
}

inline iRRAM::LAZY_BOOLEAN operator!=(const Vector &x, const Vector &y)
{
	return (x - y) * (x - y) != 0;
}

}

#endif /* euclidean.h */
