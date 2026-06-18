#ifndef JGAP_VALUEPTR_HPP
#define JGAP_VALUEPTR_HPP

#include <memory>
#include <utility>
#include <concepts>

namespace jgap {

    template <typename T>
    concept Cloneable = requires(const T& obj) {
            // Requires obj.clone() to be valid, and its return type
            // must be implicitly convertible to std::unique_ptr<T>
            { obj.clone() } -> std::convertible_to<std::unique_ptr<T>>;
    };

    template <Cloneable T>
    class ValuePtr {
    public:
        using element_type = T;

        // Constructors
        constexpr ValuePtr() noexcept = default;
        constexpr ValuePtr(std::nullptr_t) noexcept {}

        // Implicit constructor from an object of type T (or derived from T)
        template <typename U>
        requires std::derived_from<std::decay_t<U>, T> &&
                 (!std::same_as<std::decay_t<U>, ValuePtr<T>>) // Prevent hijacking the copy/move constructors
        ValuePtr(U&& obj) noexcept(std::is_nothrow_constructible_v<std::decay_t<U>, U>)
            : ptr(std::make_unique<std::decay_t<U>>(std::forward<U>(obj))) {}


        explicit ValuePtr(std::unique_ptr<T> p) noexcept
            : ptr(std::move(p)) {}

        // Deep-copy constructor
        ValuePtr(const ValuePtr& other)
            : ptr(other.ptr ? other.ptr->clone() : nullptr) {}

        // Move operations
        ValuePtr(ValuePtr&&) noexcept = default;
        ValuePtr& operator=(ValuePtr&&) noexcept = default;

        // Copy assignment (copy-and-swap)
        ValuePtr& operator=(const ValuePtr& other)
        {
            if (this != &other) {
                ValuePtr tmp(other);
                swap(tmp);
            }
            return *this;
        }

        ~ValuePtr() = default;

        // Observers
        T* get() noexcept { return ptr.get(); }
        const T* get() const noexcept { return ptr.get(); }

        T& operator*() noexcept { return *ptr; }
        const T& operator*() const noexcept { return *ptr; }

        T* operator->() noexcept { return ptr.get(); }
        const T* operator->() const noexcept { return ptr.get(); }

        explicit operator bool() const noexcept
        {
            return static_cast<bool>(ptr);
        }

        template <typename U>
        U* as() const {
            return dynamic_cast<U*>(ptr.get());
        }

        // Modifiers
        void reset(std::unique_ptr<T> p = nullptr) noexcept
        {
            ptr = std::move(p);
        }

        [[nodiscard]]
        std::unique_ptr<T> release() noexcept
        {
            return std::move(ptr);
        }

        void swap(ValuePtr& other) noexcept
        {
            ptr.swap(other.ptr);
        }

    private:
        std::unique_ptr<T> ptr;
    };

    template <typename T>
    void swap(ValuePtr<T>& a, ValuePtr<T>& b) noexcept
    {
        a.swap(b);
    }
}

#endif