#ifndef FOCSS_VECTOR_H_
#define FOCSS_VECTOR_H_

#include <cassert>

namespace focss {
template <typename Scalar>
class Vector {
    int size_;
    Scalar* data_;
    bool proxy_instance_;

  public:
    // ------------------------------------------------------------
    // Vector<Scalar> constructors and destructor
    // ------------------------------------------------------------
    inline Vector() : size_(0), data_(nullptr), proxy_instance_(false) {}

    inline explicit Vector(const int& size) {
        assert(size >= 0);

        size_ = size;
        data_ = new Scalar[size_];
        proxy_instance_ = false;
        for (int i = 0; i < size_; ++i)
            data_[i] = Scalar();
    }

    inline explicit Vector(const int& size, const Scalar& fill_value) {
        assert(size >= 0);

        size_ = size;
        data_ = new Scalar[size_];
        proxy_instance_ = false;
        for (int i = 0; i < size_; ++i)
            data_[i] = fill_value;
    }

    inline Vector(const Vector& other) {
        size_ = other.size_;
        data_ = new Scalar[size_];
        proxy_instance_ = false;
        for (int i = 0; i < size_; ++i)
            data_[i] = other.data_[i];
    }

    inline Vector(Vector&& other) {
        size_ = other.size_;
        data_ = other.data_;
        proxy_instance_ = other.proxy_instance_;

        other.size_ = 0;
        other.data_ = nullptr;
    }

    inline Vector& operator=(const Vector& other) {
        if (this != &other) {
            if (size_ != other.size_) {
                assert(proxy_instance_);

                delete[] data_;
                size_ = other.size_;
                data_ = new Scalar[size_];
            }

            for (int i = 0; i < size_; ++i)
                data_[i] = other.data_[i];
        }

        return *this;
    }

    inline Vector& operator=(Vector&& other) {
        if (this != &other) {
            if (proxy_instance_) {
                assert(other.proxy_instance_);
            } else {
                delete[] data_;
            }

            size_ = other.size_;
            data_ = other.data_;
            proxy_instance_ = other.proxy_instance_;

            other.size_ = 0;
            other.data_ = nullptr;
        }

        return *this;
    }

    inline virtual ~Vector() {
        if (!proxy_instance_) delete[] data_;
    }

  public:
    // ------------------------------------------------------------
    // Vector<Scalar> proxy handling
    // ------------------------------------------------------------
    inline bool is_proxy() const { return proxy_instance_; }

    inline static Vector proxy(const int& size, Scalar* data) {
        Vector proxy;
        proxy.size_ = size;
        proxy.data_ = data;
        proxy.proxy_instance_ = true;

        return proxy;
    }

  public:
    // ------------------------------------------------------------
    // Vector<Scalar> multiplication
    // ------------------------------------------------------------
    inline Vector& operator*=(const Scalar& multiplier) {
        for (int i = 0; i < size_; ++i)
            data_[i] *= multiplier;

        return *this;
    }

    inline Vector& operator*=(const Vector& multipliers) {
        assert(size_ == multipliers.size_);

        for (int i = 0; i < size_; ++i)
            data_[i] *= multipliers.data_[i];

        return *this;
    }

    inline Vector operator*(const Scalar& multiplier) const {
        return Vector(*this) *= multiplier;
    }

    inline Vector operator*(const Vector& multipliers) const {
        return Vector(*this) *= multipliers;
    }

    // ------------------------------------------------------------
    // Vector<Scalar> division
    // ------------------------------------------------------------
    inline Vector& operator/=(const Scalar& devisor) {
        for (int i = 0; i < size_; ++i)
            data_[i] /= devisor;

        return *this;
    }

    inline Vector& operator/=(const Vector& devisors) {
        assert(size_ == devisors.size_);

        for (int i = 0; i < size_; ++i)
            data_[i] /= devisors.data_[i];

        return *this;
    }

    inline Vector operator/(const Scalar& devisor) const {
        return Vector(*this) /= devisor;
    }

    inline Vector operator/(const Vector& devisors) const {
        return Vector(*this) /= devisors;
    }

    // ------------------------------------------------------------
    // Vector<Scalar> addition
    // ------------------------------------------------------------
    inline Vector& operator+=(const Scalar& summand) {
        for (int i = 0; i < size_; ++i)
            data_[i] += summand;

        return *this;
    }

    inline Vector& operator+=(const Vector& summands) {
        assert(size_ == summands.size_);

        for (int i = 0; i < size_; ++i)
            data_[i] += summands.data_[i];

        return *this;
    }

    inline Vector operator+(const Scalar& summand) const {
        return Vector(*this) += summand;
    }

    inline Vector operator+(const Vector& summands) const {
        return Vector(*this) += summands;
    }

    // ------------------------------------------------------------
    // Vector<Scalar> substraction
    // ------------------------------------------------------------
    inline Vector& operator-=(const Scalar& subtrahend) {
        for (int i = 0; i < size_; ++i)
            data_[i] -= subtrahend;

        return *this;
    }

    inline Vector& operator-=(const Vector& subtrahends) {
        assert(size_ == subtrahends.size_);

        for (int i = 0; i < size_; ++i)
            data_[i] -= subtrahends.data_[i];

        return *this;
    }

    inline Vector operator-(const Scalar& subtrahend) const {
        return Vector(*this) -= subtrahend;
    }

    inline Vector operator-(const Vector& subtrahends) const {
        return Vector(*this) -= subtrahends;
    }

  public:
    // ------------------------------------------------------------
    // Vector<Scalar> general methods
    // ------------------------------------------------------------
    inline int size() const { return size_; }

    inline Vector sub(const int& begin, const int& end) const {
        assert(0 <= begin);
        assert(begin < end);
        assert(end <= size_);

        Vector<Scalar> subvector(end - begin);
        for (int i = 0; i < subvector.size_; ++i)
            subvector.data_[i] = data_[i + begin];

        return subvector;
    }

    inline Vector chomp(const int& at_begin, const int& at_end) const {
        return sub(at_begin, size_ - at_end);
    }

    inline Scalar* raw() { return data_; }

    inline Scalar& operator[](const int& index) {
        assert(0 <= index && index < size_);
        return data_[index];
    }

    inline Scalar operator[](const int& index) const {
        assert(0 <= index && index < size_);
        return data_[index];
    }
};
}  // namespace focss

#endif  // FOCSS_VECTOR_H_
