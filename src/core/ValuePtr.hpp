#ifndef JGAP_VALUEPTR_HPP
#define JGAP_VALUEPTR_HPP

#include <memory>
#include <utility>
#include <concepts>
#include <type_traits>

namespace jgap {

    /**
     * Ensures one can deep-copy an object.
     *
     * @note Using raw pointers since smart pointers cannot "fall back"
     * into base class pointers upon rewriting,
     * which would either prevent non-base class {@link ValuePtr}
     * or would require too big of a boilerplate.
     */
    template <typename T>
    concept Cloneable = requires(const T& obj) {
            { obj.clone() } -> std::convertible_to<T*>;
    };

    /**
     * @brief Nullable implementation \approx {@see std::polymorphic from C++26}.
     * 
     * A smart pointer that makes a deep copy unless ownership is transferred.
     *
     * @note Intended to simplify the construction of objects that are relatively small in memory,
     * can benefit from template deduction
     * (which would be prevented by {@link std::unique_ptr}),
     * but need to be cast to a virtual base class,
     * all while maintaining memory safety (unlike using "new").
     * So e.g.:
     * {@snippet :
     * std::unique_ptr<Base> A = std::make_unique<Derived<T1, T2, T3..>>(T1 a, T2 b...);
     * }
     * would become
     * {@snippet :
     * ValuePtr<Base> A = Derived(T1 a, T2 b...);
     * }.
     */
    template <Cloneable T>
    class ValuePtr {
    public:
        using element_type = T;

        constexpr ValuePtr() noexcept = default;
        constexpr ValuePtr(std::nullptr_t) noexcept {}

        // Implicit constructor from T, as opposed to a ValuePtr<T> reference
        template <typename U>
        requires std::derived_from<
            std::decay_t<U>, // Strips away references (&, &&),
            T
        > && (!std::same_as<std::decay_t<U>, ValuePtr<T>>) // Avoid this constructor to be called on ValuePtr<U>
        ValuePtr(U&& obj) noexcept(std::is_nothrow_constructible_v<std::decay_t<U>, U>)
            : ptr(std::make_unique<std::decay_t<U>>(std::forward<U>(obj))) {}

        // Dangerous to allow implicit.
        explicit ValuePtr(std::unique_ptr<T> p) noexcept
            : ptr(std::move(p)) {}

        ValuePtr(const ValuePtr& other)
            : ptr(other.ptr ? std::unique_ptr<T>(other.ptr->clone()) : nullptr) {}

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

        // Is nullptr?
        explicit operator bool() const noexcept
        {
            return static_cast<bool>(ptr);
        }

        /**
         * Tries dynamic_cast<U*> on the underlying object's pointer.
         * @returns !WARN! original poineter still owned by this ValuePtr.
         */
        template <typename U>
        U* as() const {
            return dynamic_cast<U*>(ptr.get());
        }

        template <typename U>
        const U* as() const {
            return dynamic_cast<const U*>(ptr.get());
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
        std::unique_ptr<T> ptr;
    };

    template <typename T>
    void swap(ValuePtr<T>& a, ValuePtr<T>& b) noexcept
    {
        a.swap(b);
    }
}

#endif