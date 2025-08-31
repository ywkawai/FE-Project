#pragma once

#include <memory>
#include <stdexcept>

namespace scale{ namespace FElib{ namespace cbind{
    struct Finalizer {
        using Fn = void(*)(void*);
        Fn fn = nullptr;
        void operator()(void* p) const noexcept { if (p && fn) fn(p); }
    };
    
    using Handle = std::unique_ptr<void, Finalizer>;
    template <class InitFn, class... Args>
    inline Handle init_handle(InitFn init, Finalizer::Fn fin, const char* err_msg, Args&&... args) {
        void* raw = init(std::forward<Args>(args)...);
        if (!raw) throw std::runtime_error(err_msg ? err_msg : "Failed to initialize a handle");
        return Handle(raw, Finalizer{fin});
    }

    template <class T, class Fn>    
    inline T get_value(const Handle& handle, Fn fn) noexcept {
        T v{};
        fn(handle.get(), &v);
        return v;
    }    

    template <class Fn>
    inline std::vector<double> fetch_vector(const Handle& handle, Fn fn, int n) {
        std::vector<double> buf(static_cast<size_t>(n));
        fn(handle.get(), buf.data(), n);
        return buf;
    }
    template <class Fn>
    inline void fill_vector(const Handle& handle, Fn fn, double* out, int n) {
        fn(handle.get(), out, n);
    }

    template <class Fn>
    inline std::vector<double> fetch_matrix(const Handle& handle, Fn fn, int nx, int ny) {
        std::vector<double> buf(static_cast<size_t>(nx) * static_cast<size_t>(ny));
        fn(handle.get(), buf.data(), nx, ny);
        return buf;
    }
    template <class Fn>
    inline void fill_matrix(const Handle& handle, Fn fn, double* out, int nx, int ny) {
        fn(handle.get(), out, nx, ny);
    }
    
}}}
