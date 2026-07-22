#ifndef JGAP_VALUEPTR_HPP
#define JGAP_VALUEPTR_HPP

#include <memory>
#include <utility>
#include <concepts>
#include <type_traits>

namespace jgap {

    /// Ensures one can deep-copy an object.
    /// \note Using a raw pointer is unsafe,
    /// but it avoids the boilerplate when dealing with multi-inheritance
    /// (e.g. with \ref EamPairFunction).
    template <typename T>
    concept Cloneable = requires(const T& obj) {
            { obj.clone() } -> std::convertible_to<T*>;
    } && std::has_virtual_destructor_v<T>;

    /// @brief Nullable r-value-like smart pointer.
    ///
    /// A smart pointer that makes a deep copy unless ownership is transferred.
    ///
    /// @note Intended to simplify the construction of objects that are relatively small in memory,
    /// can benefit from template deduction
    /// (which would be prevented by \ref std::unique_ptr),
    /// but need to be cast to a virtual base class,
    /// all while maintaining memory safety (unlike using "new").
    /// So e.g.:
    /// @code
    /// std::unique_ptr<TVirtual> A = std::make_unique<Derived<T1, T2, T3..>>(T1 a, T2 b...);
    /// @endcode
    /// could be simplified:
    /// @code
    /// ValuePtr<TVirtual> A = Derived(T1 a, T2 b...);
    /// @endcode
    template <Cloneable T>
    class ValuePtr {
    public:
        using element_type = T;

        constexpr ValuePtr() noexcept = default;
        constexpr ValuePtr(std::nullptr_t) noexcept {}

        /// Universal constructor from underlying type U (creates a new allocation)
        template <typename U>
        requires std::derived_from<std::remove_cvref_t<U>, T> &&
                 (!std::same_as<std::remove_cvref_t<U>, ValuePtr<T>>)
        ValuePtr(U&& obj)
            : ptr(std::make_unique<std::remove_cvref_t<U>>(std::forward<U>(obj))) {}

        /// Cross-type Move Constructor (steals the pointer, no allocation)
        template <typename U>
        requires std::derived_from<U, T> && (!std::same_as<U, T>)
        ValuePtr(ValuePtr<U>&& obj) noexcept
            : ptr(obj.release()) {}

        /// 3. Cross-type Copy Constructor (clones the underlying object)
        template <typename U>
        requires std::derived_from<U, T> && (!std::same_as<U, T>)
        ValuePtr(const ValuePtr<U>& obj)
            : ptr(obj ? std::unique_ptr<T>(obj->clone()) : nullptr) {}

        /// Explicit constructor from unique_ptr
        explicit ValuePtr(std::unique_ptr<T> p) noexcept : ptr(std::move(p)) {}

        /// Standard Copy Constructor
        ValuePtr(const ValuePtr& other)
            : ptr(other.ptr ? std::unique_ptr<T>(other.ptr->clone()) : nullptr) {}

        /// Standard Move Constructor
        ValuePtr(ValuePtr&&) noexcept = default;
        ValuePtr& operator=(ValuePtr&&) noexcept = default;

        // Copy assignment via copy-and-swap idiom
        ValuePtr& operator=(const ValuePtr& other)
        {
            if (this != &other) {
                ValuePtr tmp(other);
                swap(tmp);
            }
            return *this;
        }

        ~ValuePtr() = default;

        T* get() noexcept { return ptr.get(); }
        const T* get() const noexcept { return ptr.get(); }

        T& operator*() noexcept { return *ptr; }
        const T& operator*() const noexcept { return *ptr; }

        T* operator->() noexcept { return ptr.get(); }
        const T* operator->() const noexcept { return ptr.get(); }

        /// Is nullptr?
        explicit operator bool() const noexcept
        {
            return static_cast<bool>(ptr);
        }

        /// Tries `dynamic_cast<U*>` on the underlying object's pointer.
        /// @returns !WARN! original pointer still owned by this ValuePtr.
        template <typename U>
        U* as() const {
            return dynamic_cast<U*>(ptr.get());
        }

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
        template <Cloneable U> friend class ValuePtr;

        std::unique_ptr<T> ptr;
    };

    template <typename T>
    void swap(ValuePtr<T>& a, ValuePtr<T>& b) noexcept
    {
        a.swap(b);
    }
}

#endif