#pragma once

#include <cmath>

template <typename T>
struct Vector
{
    T       x;
    T       y;

    Vector(T x_, T y_) : x(x_), y(y_) {}
    Vector() : x(0), y(0) {}
    Vector(const Vector& v) : x(v.x), y(v.y) {}

    Vector& operator=(const Vector& v)
    {
        if (this != &v) {
            x = v.x;
            y = v.y;
        }
        return *this;
    }

    Vector operator+(const Vector& v) const
    {
        return Vector(x + v.x, y + v.y);
    }

    Vector operator-(const Vector& v) const
    {
        return Vector(x - v.x, y - v.y);
    }

    Vector operator*(const T& s) const
    {
        return Vector(x * s, y * s);
    }

    Vector operator/(const T& s) const
    {
        return Vector(x / s, y / s);
    }

    T operator*(const Vector& v) const
    {
        return x * v.x + y * v.y;
    }

    Vector operator&(const Vector& v) const
    {
        return Vector(x * v.x, y * v.y);
    }

    T operator^(const Vector& v) const
    {
        return x * v.y - y * v.x;
    }

    T length() const
    {
        return std::sqrt(x * x + y * y);
    }

    Vector norm() const
    {
        return *this / length();
    }

    Vector operator-() const {
        return Vector(-x, -y);
    }

    Vector operator+=(const Vector& v)
    {
        x += v.x;
        y += v.y;
        return *this;
    }

    Vector operator-=(const Vector& v)
    {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    Vector operator*=(const T& s)
    {
        x *= s;
        y *= s;
        return *this;
    }

    Vector operator/=(const T& s)
    {
        x /= s;
        y /= s;
        return *this;
    }


    T& operator[](int i)
    {
        switch (i)
        {
        case 0:
            return x;
            break;
        case 1:
            return y;
            break;
        default:
            assert(false);
            break;
        }
    }

    const T& operator[](int i) const
    {
        switch (i)
        {
        case 0:
            return x;
            break;
        case 1:
            return y;
            break;
        default:
            assert(false);
            break;
        }
    }

};

template <typename T>
Vector<T> operator*(const T& s, const Vector<T>& v)
{
    return v * s;
}