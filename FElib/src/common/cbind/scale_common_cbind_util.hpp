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
}}}