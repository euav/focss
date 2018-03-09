#ifndef FOCSS_VECTOR_H_
#define FOCSS_VECTOR_H_

#include <cassert>
#include <initializer_list>

namespace focss {
template <typename Scalar>
class Vector {
    int size_;
    Scalar* data_;
    bool owner_;

  public:
    // ------------------------------------------------------------
    // Vector<Scalar> constructors and destructor
    // ------------------------------------------------------------
    inline Vector() : size_(0), data_(nullptr), owner_(true) {}

    inline explicit Vector(const int& size) {
        assert(size >= 0);

        size_ = size;
        data_ = new Scalar[size_];
        owner_ = true;
        for (int i = 0; i < size_; ++i)
            data_[i] = Scalar();
    }

    inline explicit Vector(const int& size, const Scalar& fill_value) {
        assert(size >= 0);

        size_ = size;
        data_ = new Scalar[size_];
        owner_ = true;
        for (int i = 0; i < size_; ++i)
            data_[i] = fill_value;
    }

    inline Vector(std::initializer_list<Scalar> data_list) {
        size_ = data_list.size();
        data_ = new Scalar[size_];
        owner_ = true;

        auto it = data_list.begin();
        for (int i = 0; i < size_; ++i)
            data_[i] = *it++;
    }

    inline Vector(const Vector& other) {
        size_ = other.size_;
        data_ = new Scalar[size_];
        owner_ = true;
        for (int i = 0; i < size_; ++i)
            data_[i] = other.data_[i];
    }

    inline Vector(Vector&& other) {
        if (other.owner_) {
            size_ = other.size_;
            data_ = other.data_;
            owner_ = true;

            other.size_ = 0;
            other.data_ = nullptr;
        } else {
            size_ = other.size_;
            data_ = new Scalar[size_];
            owner_ = true;

            for (int i = 0; i < size_; ++i)
                data_[i] = other.data_[i];
        }
    }

    inline Vector& operator=(const Vector& other) {
        if (this != &other) {
            if (size_ != other.size_) {
                assert(owner_);

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
            if (owner_ && other.owner_) {
                delete[] data_;

                size_ = other.size_;
                data_ = other.data_;
                owner_ = other.owner_;

                other.size_ = 0;
                other.data_ = nullptr;
            } else {
                if (size_ != other.size_) {
                    assert(owner_);

                    delete[] data_;
                    size_ = other.size_;
                    data_ = new Scalar[size_];
                }

                for (int i = 0; i < size_; ++i)
                    data_[i] = other.data_[i];
            }
        }

        return *this;
    }

    inline virtual ~Vector() {
        if (owner_) delete[] data_;
    }

  public:
    // ------------------------------------------------------------
    // Vector<Scalar> weak ownership handling
    // ------------------------------------------------------------
    inline bool is_owner() const { return owner_; }

    inline static Vector proxy(const int& size, Scalar* data) {
        Vector proxy;
        proxy.size_ = size;
        proxy.data_ = data;
        proxy.owner_ = false;

        return proxy;
    }

  public:
    // ------------------------------------------------------------
    // Vector<Scalar> access methods
    // ------------------------------------------------------------
    inline int size() const { return size_; }

    inline Scalar* raw() { return data_; }

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

    inline Scalar& operator[](const int& index) {
        assert(0 <= index && index < size_);
        return data_[index];
    }

    inline Scalar operator[](const int& index) const {
        assert(0 <= index && index < size_);
        return data_[index];
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
};
}  // namespace focss

#endif  // FOCSS_VECTOR_H_
