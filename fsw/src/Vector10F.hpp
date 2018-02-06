#ifndef Vector10F_HPP
#define Vector10F_HPP

#include "cfe.h"

namespace math
{

class Vector10F
{
private:
	static const int SIZE = 10;
	float data[SIZE];
	float nan;

public:
	Vector10F();
	Vector10F(float A, float B, float C, float D, float E, float F, float G, float H, float I, float J);
	~Vector10F();
	float& operator [] (uint32 i);
	float operator [] (uint32 i) const;
	float Length(void) const;
	void Zero(void);
	Vector10F operator+(const Vector10F &vecIn);
	Vector10F operator+(const Vector10F &vecIn) const;
	Vector10F operator-(const Vector10F &vecIn);
	Vector10F operator-(const Vector10F &vecIn) const;
	//Vector10F operator%(const Vector10F &vecIn);
	Vector10F operator*(const float scalar);
	Vector10F operator/(const float scalar);
	float operator*(const Vector10F &vecIn);
	Vector10F operator-() const;
	Vector10F EMult(const Vector10F &vecIn);
	void Constrain(uint32 i, float min, float max);
	Vector10F Normalized(void);
	void Normalize(void);
	const Vector10F EDivide(const Vector10F &v) const;

protected:

};

}

#endif
