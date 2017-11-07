#ifndef FOCSS_VECTOR_H_
#define FOCSS_VECTOR_H_

namespace focss {
template <typename Scalar>
class Vector {
    int size_;
    Scalar* data_;
    bool proxy_instance_;

  public:
    Vector();
    explicit Vector(const int& size);
    explicit Vector(const int& size, const Scalar& fill_value);

    Vector(const Vector& other);
    Vector(Vector&& other);

    Vector& operator=(const Vector& other);
    Vector& operator=(Vector&& other);

    ~Vector();

  public:
    bool is_proxy() const;
    static Vector proxy(const int& size, Scalar* data);

  public:
    Vector& operator*=(const Scalar& multiplier);
    Vector& operator*=(const Vector& multipliers);
    Vector operator*(const Scalar& multiplier) const;
    Vector operator*(const Vector& multipliers) const;

    Vector& operator+=(const Scalar& summand);
    Vector& operator+=(const Vector& summands);
    Vector operator+(const Scalar& summand) const;
    Vector operator+(const Vector& summands) const;

    Vector& operator-=(const Scalar& subtrahend);
    Vector& operator-=(const Vector& subtrahends);
    Vector operator-(const Scalar& subtrahend) const;
    Vector operator-(const Vector& subtrahends) const;

  public:
    int size() const;
    Vector sub(const int& begin, const int& end) const;
    Vector chomp(const int& at_begin, const int& at_end) const;

    Scalar* raw();
    Scalar& operator[](const int& index);
    Scalar operator[](const int& index) const;
};
}  // namespace focss

#endif  // FOCSS_VECTOR_H_
