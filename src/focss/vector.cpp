#include "vector.h"
#include <cassert>

namespace focss {
// ----------------------------------------------------------------------
// -------------------- Vector<Scalar> constructors and destructor
// ----------------------------------------------------------------------
template <typename Scalar>
Vector<Scalar>::Vector() : size_(0), data_(nullptr), proxy_instance_(false) {}

template <typename Scalar>
Vector<Scalar>::Vector(const int& size) : Vector(size, 0) {}

template <typename Scalar>
Vector<Scalar>::Vector(const int& size, const Scalar& fill_value) {
    assert(size > 0);

    size_ = size;
    data_ = new Scalar[size_];
    proxy_instance_ = false;
    for (int i = 0; i < size_; ++i)
        data_[i] = fill_value;
}

template <typename Scalar>
Vector<Scalar>::Vector(const Vector<Scalar>& other) {
    size_ = other.size_;
    data_ = new Scalar[size_];
    proxy_instance_ = false;

    for (int i = 0; i < size_; ++i)
        data_[i] = other.data_[i];
}

template <typename Scalar>
Vector<Scalar>::Vector(Vector&& other) {
    if (other.proxy_instance_) {
        size_ = other.size_;
        data_ = new Scalar[size_];
        proxy_instance_ = false;

        for (int i = 0; i < size_; ++i)
            data_[i] = other.data_[i];
    } else {
        size_ = other.size_;
        data_ = other.data_;
        proxy_instance_ = false;

        other.size_ = 0;
        other.data_ = nullptr;
    }
}

template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator=(const Vector<Scalar>& other) {
    if (this != &other) {
        if (proxy_instance_) assert(size_ == other.size_);

        if (size_ != other.size__) {
            delete[] data_;
            size_ = other.size_;
            data_ = new Scalar[size_];
        }

        for (int i = 0; i < size_; ++i)
            data_[i] = other.data_[i];
    }

    return *this;
}

template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator=(Vector&& other) {
    if (this != &other) {
        if (proxy_instance_) {
            assert(size_ == other.size_);

            for (int i = 0; i < size_; ++i)
                data_[i] = other.data_[i];
        } else if (!proxy_instance_ && other.proxy_instance_) {
            if (size_ != other.size_) {
                delete[] data_;
                size_ = other.size_;
                data_ = new Scalar[size_];
            }

            for (int i = 0; i < size_; ++i)
                data_[i] = other.data_[i];
        } else {
            delete[] data_;
            size_ = other.size_;
            data_ = other.data_;

            other.size_ = 0;
            other.data_ = nullptr;
        }
    }

    return *this;
}

template <typename Scalar>
Vector<Scalar>::~Vector() {
    if (!proxy_instance_) delete[] data_;
}

// ----------------------------------------------------------------------
// -------------------- Vector<Scalar> proxy handling
// ----------------------------------------------------------------------
template <typename Scalar>
bool Vector<Scalar>::is_proxy() const {
    return proxy_instance_;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::proxy(const int& size, Scalar* data) {
    Vector<Scalar> proxy;
    proxy.size_ = size;
    proxy.data_ = data;
    proxy.proxy_instance_ = true;

    return proxy;
}

// ----------------------------------------------------------------------
// -------------------- Vector<Scalar> multiplication
// ----------------------------------------------------------------------
template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator*=(const Scalar& multiplier) {
    for (int i = 0; i < size_; ++i)
        data_[i] *= multiplier;

    return *this;
}

template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator*=(const Vector<Scalar>& multipliers) {
    assert(size_ == multipliers.size_);

    for (int i = 0; i < size_; ++i)
        data_[i] *= multipliers.data_[i];

    return *this;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::operator*(const Scalar& multiplier) const {
    return Vector(*this) *= multiplier;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::operator*(const Vector<Scalar>& multipliers) const {
    return Vector(*this) *= multipliers;
}

// ----------------------------------------------------------------------
// -------------------- Vector<Scalar> addition
// ----------------------------------------------------------------------
template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator+=(const Scalar& summand) {
    for (int i = 0; i < size_; ++i)
        data_[i] += summand;
    http://www.cpo
    return *this;
}

template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator+=(const Vector<Scalar>& summands) {
    assert(size_ == summands.size_);

    for (int i = 0; i < size_; ++i)
        data_[i] += summands.data_[i];

    return *this;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::operator+(const Scalar& summand) const {
    return Vector(*this) += summand;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::operator+(const Vector<Scalar>& summands) const {
    return Vector(*this) += summands;
}

// ----------------------------------------------------------------------
// -------------------- Vector<Scalar> substraction
// ----------------------------------------------------------------------
template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator-=(const Scalar& subtrahend) {
    for (int i = 0; i < size_; ++i)
        data_[i] -= subtrahend;

    return *this;
}

template <typename Scalar>
Vector<Scalar>& Vector<Scalar>::operator-=(const Vector<Scalar>& subtrahends) {
    assert(size_ == subtrahends.size_);

    for (int i = 0; i < size_; ++i)
        data_[i] -= subtrahends.data_[i];

    return *this;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::operator-(const Scalar& subtrahend) const {
    return Vector(*this) -= subtrahend;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::operator-(const Vector<Scalar>& subtrahends) const {
    return Vector(*this) -= subtrahends;
}

// ----------------------------------------------------------------------
// -------------------- Vector<Scalar> general methods
// ----------------------------------------------------------------------
template <typename Scalar>
int Vector<Scalar>::size() const {
    return size_;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::sub(const int& begin, const int& end) const {
    assert(0 <= begin);
    assert(begin < end);
    assert(end <= size_);

    Vector<Scalar> subvector(end - begin);
    for (int i = 0; i < subvector.size_; ++i)
        subvector.data_[i] = data_[i + begin];

    return subvector;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::chomp(const int& at_begin, const int& at_end) const {
    return sub(at_begin, size_ - at_end);
}

template <typename Scalar>
Scalar* Vector<Scalar>::raw() {
    return data_;
}

template <typename Scalar>
Scalar& Vector<Scalar>::operator[](const int& index) {
    assert((0 <= index && index < size_));
    return data_[index];
}

template <typename Scalar>
Scalar Vector<Scalar>::operator[](const int& index) const {
    assert((0 <= index && index < size_));
    return data_[index];
}
}  // namespace focss
